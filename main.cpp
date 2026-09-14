/*
 * hisat3n unified dispatcher.
 *
 * Single binary serving both align (hisat2) and index-build (hisat2_build)
 * routes.  Routing is selected either by an explicit subcommand
 *     hisat3n align ...     (alias: hisat3n-build? no - use argv0 below)
 *     hisat3n build <index> ...
 * or, for backward compatibility with the historical multi-binary layout,
 * by the program name used to invoke it:
 *     hisat-3n  -> align   (was hisat2-align-s / hisat2-align-l)
 *     hisat-3n-build -> build (was hisat2-build-s / hisat2-build-l)
 *
 * The -A <file> batch mode (one arg set per line) is preserved here so it
 * applies uniformly to whichever route is selected.
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include "tokenize.h"
#include "ds.h"
#include "mem_ids.h"

using namespace std;

extern "C" {
	int hisat2(int argc, const char **argv);
	int hisat2_build(int argc, const char **argv);
}

int proc_align(int argc, const char **argv);

enum struct Route { ALIGN, BUILD, UNKNOWN };

static Route detect_route(int argc, char **argv) {
	// 1) explicit subcommand
	if(argc >= 2) {
		if(strcmp(argv[1], "align") == 0) return Route::ALIGN;
		if(strcmp(argv[1], "build") == 0) return Route::BUILD;
	}
	// 2) argv[0] shim for backward compatibility
	const char *base = strrchr(argv[0], '/');
	base = base ? base + 1 : argv[0];
	if(strstr(base, "build") != NULL)  return Route::BUILD;
	if(strstr(base, "align") != NULL)  return Route::ALIGN;
	if(strstr(base, "3n") != NULL &&
	   strstr(base, "build") == NULL)  return Route::ALIGN; // hisat-3n / hisat3n
	return Route::UNKNOWN;
}

typedef int (*route_fn)(int, const char **);

static int dispatch(Route route, int argc, const char **argv) {
	route_fn fn = NULL;
	switch(route) {
		case Route::ALIGN: fn = hisat2;       break;
		case Route::BUILD: fn = hisat2_build; break;
		default:                               return 2;
	}

	// -A <file>: one arg set per line, dispatched to the selected route.
	// argv[1] is the -A flag literally; after subcommand stripping a real
	// user arg in position 1 is never "align"/"build", so this always fires
	// for the batch spelling regardless of entry style.
	if(argc > 2 && strcmp(argv[1], "-A") == 0) {
		const char *file = argv[2];
		ifstream in;
		in.open(file);
		char buf[4096];
		int lastret = -1;
		while(in.getline(buf, 4095)) {
			EList<string> args(MISC_CAT);
			args.push_back(string(argv[0]));
			tokenize(buf, " \t", args);
			const char **myargs = (const char**)malloc(sizeof(char*)*args.size());
			for(size_t i = 0; i < args.size(); i++) {
				myargs[i] = args[i].c_str();
			}
			if(args.size() == 1) continue;
			lastret = fn((int)args.size(), myargs);
			free(myargs);
		}
		if(lastret == -1) {
			cerr << "Warning: No arg strings parsed from " << file << endl;
			return 0;
		}
		return lastret;
	}
	return fn(argc, argv);
}

static void print_usage(const char *prog) {
	cerr << "Usage: " << prog << " <align|build> [options...]" << endl
	     << "  align   map nucleotide-conversion reads (was hisat-3n)" << endl
	     << "  build   build a 3N index from a reference (was hisat-3n-build)" << endl;
}

static bool has_proc_parallel(int argc, char **argv, int *n) {
	for(int i = 1; i < argc; i++) {
		if(strncmp(argv[i], "--proc-parallel=", 16) == 0) {
			*n = atoi(argv[i] + 16); return true;
		}
		if(strcmp(argv[i], "--proc-parallel") == 0 && i + 1 < argc) {
			*n = atoi(argv[i + 1]); return true;
		}
	}
	return false;
}

// Strip --proc-parallel[=]N from argv when it cannot be honored (N <= 1 or
// otherwise non-positive), so the default single-process path tolerates it
// being present but harmless.
static int strip_proc_parallel(int argc, const char **argv) {
	vector<const char*> kept;
	kept.push_back(argv[0]);
	for(int i = 1; i < argc; i++) {
		if(strncmp(argv[i], "--proc-parallel=", 16) == 0) continue;
		if(strcmp(argv[i], "--proc-parallel") == 0) { i++; continue; }
		kept.push_back(argv[i]);
	}
	for(size_t i = 0; i < kept.size(); i++) {
		const_cast<char**>(argv)[i] = const_cast<char*>(kept[i]);
	}
	return (int)kept.size();
}

int main(int argc, char **argv) {
	Route route = detect_route(argc, argv);
	if(route == Route::ALIGN) {
		int nproc = 0;
		if(has_proc_parallel(argc, argv, &nproc) && nproc > 1) {
			return proc_align(argc, (const char**)argv);
		}
		if(has_proc_parallel(argc, argv, &nproc)) {
			argc = strip_proc_parallel(argc, (const char**)argv);
		}
	}
	if(route == Route::UNKNOWN) {
		// No subcommand, no known argv[0]: if argv[1] starts with '-' it is a
		// bare align invocation (e.g. "hisat3n -x idx -U r.fq"), so default to
		// align; otherwise print usage.
		if(argc > 1 && argv[1][0] == '-') {
			route = Route::ALIGN;
		} else {
			print_usage(argv[0]);
			return 1;
		}
	} else if(argc >= 2 &&
	          (strcmp(argv[1], "align") == 0 || strcmp(argv[1], "build") == 0)) {
		// strip the explicit subcommand token before forwarding
		argv++; argc--;
	}
	return dispatch(route, argc, (const char**)argv);
}