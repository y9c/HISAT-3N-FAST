/*
 * Process-parallel alignment for HISAT-3N.
 *
 * A single invocation can be stretched across N worker OS processes by
 * splitting the read input at record boundaries, aligning each chunk in a
 * separate worker (with a fraction of the threads), and merging the per-chunk
 * SAM output in chunk order.  Because every read is a self-contained
 * alignment (no cross-read state is carried), the merged result is equivalent
 * to a single-process run on the same input (same alignment set; SAM header taken
 * from the first worker).
 *
 * This is opt-in via --proc-parallel <N>; without it main.cpp never routes
 * here and the single-process path runs unchanged.
 */

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <spawn.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

extern char **environ;

using namespace std;

extern "C" {
	int hisat2(int argc, const char **argv);
}

static string find_short_arg(int argc, const char **argv, char opt, int *idx) {
	for(int i = 1; i + 1 < argc; i++) {
		if(strlen(argv[i]) == 2 && argv[i][0] == '-' && argv[i][1] == opt) {
			if(idx != NULL) *idx = i;
			return argv[i + 1];
		}
	}
	return string();
}

static string find_long_arg(int argc, const char **argv, const char* name, int *idx) {
	int len = (int)strlen(name);
	for(int i = 1; i + 1 < argc; i++) {
		if(strcmp(argv[i], name) == 0) {
			if(idx != NULL) *idx = i;
			return argv[i + 1];
		}
		if(strncmp(argv[i], name, len) == 0 && argv[i][len] == '=') {
			if(idx != NULL) *idx = i;
			return argv[i] + len + 1;
		}
	}
	return string();
}

static int find_threads(int argc, const char **argv, int fallback) {
	int idx = -1;
	string v = find_short_arg(argc, argv, 'p', &idx);
	if(v.empty()) return fallback;
	int t = atoi(v.c_str());
	return t > 0 ? t : fallback;
}

static bool is_fastq(const string& in) {
	FILE* f = fopen(in.c_str(), "rb");
	if(f == NULL) return false;
	char c;
	bool fastq = false;
	if(fread(&c, 1, 1, f) == 1 && c == '@') fastq = true;
	fclose(f);
	return fastq;
}

static bool split_reads(const string& in, int nproc, vector<string>& chunks, bool fastq) {
	chunks.clear();
	const char* cdir = getenv("HISAT3N_TMPDIR");
	if(!cdir || !*cdir) cdir = "/tmp";
	char tmpl[128];
	for(int k = 0; k < nproc; k++) {
		snprintf(tmpl, sizeof(tmpl), "%s/hisat3n_chunk_%d_%d.fa", cdir, (int)getpid(), k);
		chunks.push_back(tmpl);
	}
	// mmap the whole input read-only; scanning and chunk writes address it
	// directly so we never copy the file into heap memory.
	int fd = open(in.c_str(), O_RDONLY);
	if(fd < 0) return false;
	struct stat sb;
	if(fstat(fd, &sb) != 0 || sb.st_size <= 0) { close(fd); return false; }
	const char* base = (const char*)mmap(NULL, (size_t)sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	close(fd);
	if(base == MAP_FAILED) return false;
	size_t n = (size_t)sb.st_size;
	vector<pair<size_t,size_t> > recs;
	if(!fastq) {
		// Vectorized FASTA scan: find '>' at line starts via memchr.
		size_t st = 0;
		const char* cur = base;
		const char* eof = base + n;
		while(cur < eof) {
			const char* gt = (const char*)memchr(cur, '>', (size_t)(eof - cur));
			if(gt == NULL) break;
			if(gt == base || *(gt - 1) == '\n') {
				if(st < (size_t)(gt - base)) recs.push_back(make_pair(st, (size_t)(gt - base)));
				st = (size_t)(gt - base);
			}
			cur = gt + 1;
		}
		if(st < n) recs.push_back(make_pair(st, n));
	} else {
		size_t line = 0, st = 0, i = 0;
		while(i < n) {
			const char* nl = (const char*)memchr(base + i, '\n', n - i);
			size_t off = nl == NULL ? n : (size_t)(nl - base);
			i = off + 1;
			line++;
			if(line % 4 == 0 && i < n) {
				recs.push_back(make_pair(st, i));
				st = i;
			}
		}
		if(st < n) recs.push_back(make_pair(st, n));
	}
	size_t nrec = recs.size();
	if(nrec == 0) { munmap((void*)base, n); return false; }
	size_t chunk_size = (nrec + nproc - 1) / nproc;
	if(chunk_size == 0) chunk_size = 1;
	size_t cur = 0;
	size_t chunk_cnt = 0;
	size_t cur_begin = recs[0].first;
	for(size_t r = 0; r < nrec; r++) {
		if(r > 0 && chunk_cnt >= chunk_size && cur + 1 < nproc) {
			size_t klen = recs[r-1].second - cur_begin;
			if(klen > 0) {
				FILE* o = fopen(chunks[cur].c_str(), "wb");
				if(o == NULL) { munmap((void*)base, n); return false; }
				fwrite(base + cur_begin, 1, klen, o);
				fclose(o);
			}
			cur++;
			chunk_cnt = 0;
			cur_begin = recs[r].first;
		}
		chunk_cnt++;
	}
	size_t klen = recs[nrec-1].second - cur_begin;
	if(klen > 0) {
		FILE* o = fopen(chunks[cur].c_str(), "wb");
		if(o == NULL) { munmap((void*)base, n); return false; }
		fwrite(base + cur_begin, 1, klen, o);
		fclose(o);
	}
	munmap((void*)base, n);
	return true;
}

static void merge_sams(const vector<string>& sam_files, const string& out_path, int nproc) {
	FILE* fo = NULL;
	if(out_path.empty()) { fo = stdout; } else { fo = fopen(out_path.c_str(), "wb"); }
	if(fo == NULL) return;
	// The first file carries the SAM header; all later files drop their own
	// header lines. Bulk-copy by mmap so we avoid 2M line-by-line reads.
	bool past_header = false;
	for(int k = 0; k < nproc; k++) {
		int fd = open(sam_files[k].c_str(), O_RDONLY);
		if(fd < 0) continue;
		struct stat sb;
		if(fstat(fd, &sb) != 0 || sb.st_size <= 0) { close(fd); continue; }
		size_t sz = (size_t)sb.st_size;
		const char* base = (const char*)mmap(NULL, sz, PROT_READ, MAP_PRIVATE, fd, 0);
		close(fd);
		if(base == MAP_FAILED) continue;
		if(!past_header) {
			// Walk only until the first non-header line of the first file.
			size_t pos = 0;
			while(pos < sz) {
				const char* nl = (const char*)memchr(base + pos, '\n', sz - pos);
				size_t line_end = nl == NULL ? sz : (size_t)(nl - base) + 1;
				if(base[pos] != '@') { past_header = true; break; }
				pos = line_end;
			}
			fwrite(base, 1, pos, fo);
			fwrite(base + pos, 1, sz - pos, fo);
		} else {
			// Later files: skip leading '@' header lines, bulk-copy the rest.
			size_t pos = 0;
			while(pos < sz && base[pos] == '@') {
				const char* nl = (const char*)memchr(base + pos, '\n', sz - pos);
				if(nl == NULL) { pos = sz; break; }
				pos = (size_t)(nl - base) + 1;
			}
			fwrite(base + pos, 1, sz - pos, fo);
		}
		munmap((void*)base, sz);
	}
	if(fo != stdout) fclose(fo);
}

int proc_align(int argc, const char **argv) {
	int nproc = 0;
	for(int i = 1; i < argc; i++) {
		if(strncmp(argv[i], "--proc-parallel=", 16) == 0) { nproc = atoi(argv[i] + 16); break; }
		if(strcmp(argv[i], "--proc-parallel") == 0 && i + 1 < argc) { nproc = atoi(argv[i + 1]); break; }
	}
	if(nproc < 2) {
		cerr << "proc-parallel needs N >= 2" << endl;
		return 1;
	}
	int uidx = -1, sidx = -1;
	string in   = find_short_arg(argc, argv, 'U', &uidx);
	string out  = find_short_arg(argc, argv, 'S', &sidx);
	if(in.empty()) {
		return hisat2(argc, argv);
	}
	int threads = find_threads(argc, argv, 1);
	bool fastq = is_fastq(in);
	vector<string> chunks;
	string presplit = find_long_arg(argc, argv, "--presplit-dir", NULL);
	bool presplit_mode = !presplit.empty();
	if(presplit_mode) {
		for(int k = 0; k < nproc; k++) {
			string fn = presplit + "/rec" + std::to_string(k) + ".fa";
			chunks.push_back(fn);
		}
	} else if(!split_reads(in, nproc, chunks, fastq)) {
		return hisat2(argc, argv);
	}
	// split_reads may write fewer chunk files than nproc when the input has
	// fewer records than workers.  Count the files that actually exist so we
	// never spawn a worker onto a missing chunk (which would make it error out
	// with "Could not open read file ..." and abort the whole run).
	{
		int nactive = 0;
		for(int k = 0; k < nproc; k++) {
			if(access(chunks[k].c_str(), F_OK) == 0) nactive++;
		}
		if(nactive < 1) return hisat2(argc, argv);
		nproc = nactive;
	}
	vector<vector<string> > store(nproc);
	vector<const char**> worker_argv(nproc);
	vector<int> worker_argc(nproc);
	vector<string> sam_files(nproc);
	vector<pid_t> pids(nproc);
	// Route per-worker intermediate SAM files to a RAM-backed dir (tmpfs) if
	// available, so both the worker writes and the final merge never touch
	// spinning/network disk.
	const char* ramdir = getenv("HISAT3N_SAM_TMPDIR");
	if(!ramdir || !*ramdir) ramdir = "/dev/shm";
	struct stat st;
	bool ram_ok = (stat(ramdir, &st) == 0 && S_ISDIR(st.st_mode));
	if(!ram_ok) ramdir = "/tmp";
	for(int k = 0; k < nproc; k++) {
		char tmp[160];
		snprintf(tmp, sizeof(tmp), "%s/hisat3n_sam_%d_%d.sam", ramdir, (int)getpid(), k);
		sam_files[k] = tmp;
	}
	for(int k = 0; k < nproc; k++) {
		vector<string>& args = store[k];
		args.push_back(argv[0]);
		for(int i = 1; i < argc; i++) {
			if(strcmp(argv[i], "--proc-parallel") == 0) { i++; continue; }
			if(strncmp(argv[i], "--proc-parallel=", 16) == 0) continue;
			if(strcmp(argv[i], "--presplit-dir") == 0) { i++; continue; }
			if(strncmp(argv[i], "--presplit-dir=", 16) == 0) continue;
			if(i == uidx && uidx >= 0) { i++; continue; }
			if(i == sidx && sidx >= 0) { i++; continue; }
			if(strcmp(argv[i], "-p") == 0) { i++; continue; }
			args.push_back(argv[i]);
		}
		// Load the index with --mm (memory-mapped, MAP_SHARED): every worker
		// mmaps the same index files, so the OS shares the physical index
		// pages across workers instead of each process allocating its own copy
		// in heap. On a big genome index this removes a large per-worker
		// memory multiplication for free, with no change to alignment output.
		args.push_back("--mm");
		// Split the requested thread budget evenly across workers so the total
		// thread count matches a single-process run (per-worker ~= threads/nproc),
		// keeping the parallel path an apples-to-apples comparison.
		args.push_back("-U"); args.push_back(chunks[k]);
		int per_worker = threads / nproc;
		int rem = threads % nproc;
		if(per_worker < 1) per_worker = 1;
		args.push_back("-p"); args.push_back(to_string(per_worker + (k < rem ? 1 : 0)));
		args.push_back("-S"); args.push_back(sam_files[k]);
		const char **a = (const char**)malloc(sizeof(char*) * (args.size() + 1));
		for(size_t j = 0; j < args.size(); j++) a[j] = args[j].c_str();
		a[args.size()] = NULL;
		worker_argv[k] = a;
		worker_argc[k] = (int)args.size();
	}
	int rc = 0;
	for(int k = 0; k < nproc; k++) {
		posix_spawn_file_actions_t fa;
		posix_spawn_file_actions_init(&fa);
		int r = posix_spawn(&pids[k], worker_argv[k][0], &fa, NULL,
			               (char* const*)worker_argv[k], environ);
		posix_spawn_file_actions_destroy(&fa);
		if(r != 0) { rc = 127; break; }
	}
	for(int k = 0; k < nproc; k++) {
		int st;
		waitpid(pids[k], &st, 0);
		if(WIFEXITED(st) && WEXITSTATUS(st) != 0) rc = WEXITSTATUS(st);
	}
	merge_sams(sam_files, out, nproc);
	for(int k = 0; k < nproc; k++) {
		if(!presplit_mode) remove(chunks[k].c_str());
		remove(sam_files[k].c_str());
		free(worker_argv[k]);
	}
	return rc;
}
