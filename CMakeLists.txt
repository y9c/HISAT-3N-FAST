cmake_minimum_required(VERSION 3.16)
project(hisat3n VERSION 3.0.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(BOWTIE_MM         "Enable memory-mapped index"              ON)
option(BOWTIE_SHARED_MEM "Enable shared-memory index"              OFF)
option(USE_SRA           "Enable NCBI SRA Toolkit input"           OFF)
option(USE_64BIT_INDEX   "Use 64-bit (large) index offsets (-l)"   ON)

if(MINGW OR CYGWIN)
    set(BOWTIE_MM OFF)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

file(STRINGS HISAT2_VERSION HISAT2_VERSION)
string(TIMESTAMP BUILD_DATE "%Y-%m-%d")

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
# Linked into both align and build (shared core)
set(CORE_CPPS
    alphabet.cpp
    ccnt_lut.cpp
    ds.cpp
    edit.cpp
    gfm.cpp
    limit.cpp
    multikey_qsort.cpp
    random_source.cpp
    reference.cpp
    ref_read.cpp
    shmem.cpp
    tinythread.cpp
)

# Align route only
set(ALIGN_CPPS
    aligner_bt.cpp
    aligner_cache.cpp
    aligner_driver.cpp
    aligner_result.cpp
    aligner_seed.cpp
    aligner_seed2.cpp
    aligner_seed_policy.cpp
    aligner_sw.cpp
    aligner_sw_driver.cpp
    aligner_swsse.cpp
    aligner_swsse_ee_i16.cpp
    aligner_swsse_ee_u8.cpp
    aligner_swsse_loc_i16.cpp
    aligner_swsse_loc_u8.cpp
    aln_sink.cpp
    dp_framer.cpp
    outq.cpp
    pat.cpp
    pe.cpp
    presets.cpp
    qual.cpp
    random_util.cpp
    read_qseq.cpp
    ref_coord.cpp
    mask.cpp
    scoring.cpp
    simple_func.cpp
    splice_site.cpp
    sse_util.cpp
    unique.cpp
    hisat2.cpp
    procmain.cpp
)

# Build route only
set(BUILD_CPPS
    diff_sample.cpp
    hisat2_build.cpp
)

# 3N extension (shared by align + build)
set(THREE_N_CPPS
    alignment_3n.cpp
    position_3n.cpp
    utility_3n.cpp
)

# Library backend (exported via hisat2lib for JNI/Python bindings)
set(HT2LIB_CPPS
    hisat2lib/ht2_init.cpp
    hisat2lib/ht2_repeat.cpp
    hisat2lib/ht2_index.cpp
)

set(ALL_SOURCES ${CORE_CPPS} ${ALIGN_CPPS} ${BUILD_CPPS} ${THREE_N_CPPS} ${HT2LIB_CPPS})

# ---------------------------------------------------------------------------
# Find dependencies
# ---------------------------------------------------------------------------
find_package(Threads REQUIRED)
find_package(OpenMP)

# ---------------------------------------------------------------------------
# OBJECT library: compiled once (PIC), shared by both libraries + the binary
# ---------------------------------------------------------------------------
add_library(hisat3n_objs OBJECT ${ALL_SOURCES})
set_target_properties(hisat3n_objs PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(hisat3n_objs
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}          # flat in-tree headers
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/hisat2lib
)

target_compile_options(hisat3n_objs PRIVATE
    $<$<CONFIG:Release>:-O3 -funroll-loops -march=native -faligned-new>
    $<$<BOOL:${OpenMP_FOUND}>:-fopenmp>
)

target_compile_definitions(hisat3n_objs PRIVATE
    BOWTIE2
    HISAT2_VERSION="${HISAT2_VERSION}"
    BUILD_TIME="${BUILD_DATE}"
    BUILD_HOST="${CMAKE_HOST_SYSTEM}"
    COMPILER_VERSION="${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}"
    COMPILER_OPTIONS="${CMAKE_CXX_FLAGS}"
    _LARGEFILE_SOURCE
    _FILE_OFFSET_BITS=64
    _GNU_SOURCE
    POPCNT_CAPABILITY
    $<$<BOOL:${BOWTIE_MM}>:BOWTIE_MM>
    $<$<BOOL:${BOWTIE_SHARED_MEM}>:BOWTIE_SHARED_MEM>
    $<$<BOOL:${USE_SRA}>:USE_SRA>
    $<$<BOOL:${USE_64BIT_INDEX}>:BOWTIE_64BIT_INDEX>
)

# ---------------------------------------------------------------------------
# Static + shared libraries (export hisat2()/hisat2_build()/ht2_* for bindings)
# ---------------------------------------------------------------------------
add_library(hisat3n_static STATIC $<TARGET_OBJECTS:hisat3n_objs>)
add_library(hisat3n_shared SHARED $<TARGET_OBJECTS:hisat3n_objs>)

foreach(target hisat3n_static hisat3n_shared)
    target_link_libraries(${target} PUBLIC Threads::Threads)
    if(OpenMP_FOUND)
        target_link_libraries(${target} PUBLIC OpenMP::OpenMP_CXX)
    endif()
    if(USE_SRA)
        if(NOT NCBI_NGS_DIR)
            message(FATAL_ERROR "USE_SRA requires -DNSRA_DIR=/path/to/ngs (via -DNCBI_NGS_DIR)")
        endif()
        target_include_directories(${target} PRIVATE "${NCBI_NGS_DIR}/include")
        target_link_libraries(${target} PUBLIC
            ${NCBI_NGS_DIR}/lib64/libncbi-ngs-c++-static.a
            ${NCBI_NGS_DIR}/lib64/libngs-c++-static.a
            ${NCBI_NGS_DIR}/lib64/libncbi-vdb-static.a
            ${CMAKE_DL_LIBS})
    endif()
endforeach()

set_target_properties(hisat3n_static PROPERTIES OUTPUT_NAME hisat3n)
set_target_properties(hisat3n_shared PROPERTIES OUTPUT_NAME hisat3n)

# Link-time optimization: cross-TU inlining in this template-heavy codebase
# (the dual-width index_t instantiations) gives a measurable runtime win.
include(CheckIPOSupported)
check_ipo_supported(RESULT lto_supported OUTPUT lto_error)
if(lto_supported)
    set_target_properties(hisat3n_objs hisat3n_static hisat3n_shared PROPERTIES
        INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ---------------------------------------------------------------------------
# Single binary: dispatches by subcommand (align | build) or argv[0] shim.
# main.cpp is the ONLY executable-source file; the library holds no main().
# ---------------------------------------------------------------------------
add_executable(hisat3n main.cpp)
target_link_libraries(hisat3n PRIVATE hisat3n_static)

# ---------------------------------------------------------------------------
# Command-name shims: historical official-style entry points.
#   hisat-3n       -> align route
#   hisat-3n-build -> build route
# main.cpp's argv[0] inspection dispatches to the matching route, so these
# are plain symlinks to the single binary (no separate executables).
# ---------------------------------------------------------------------------
add_custom_command(TARGET hisat3n POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_NAME:hisat3n>
            $<TARGET_FILE_DIR:hisat3n>/hisat-3n
    COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_NAME:hisat3n>
            $<TARGET_FILE_DIR:hisat3n>/hisat-3n-build
    COMMENT "Creating hisat-3n / hisat-3n-build symlinks")

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS hisat3n hisat3n_static hisat3n_shared
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(CODE "
    execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
        hisat3n \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/hisat-3n\")
    execute_process(COMMAND \"${CMAKE_COMMAND}\" -E create_symlink
        hisat3n \"${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/hisat-3n-build\")
")