cmake_minimum_required(VERSION 2.0)
project(spaced_seed)
include_directories(
    "/home/ming/software/gtest-1.6.0"
    "/home/ming/software/gtest-1.6.0/include"
    "${CMAKE_SOURCE_DIR}/src" 
    "${CMAKE_SOURCE_DIR}/test"
)
add_library(gtest
    /home/ming/software/gtest-1.6.0/src/gtest-all.cc
)
add_library(gtest_main
    /home/ming/software/gtest-1.6.0/src/gtest_main.cc
)
add_executable(
    src/spaced_seed 
    src/spaced_seed.cpp 
    src/common.h 
    src/seq_aligner.h 
    src/dna_seq.h
)
add_executable(
    src/visual_align 
    src/visual_align.cpp 
    src/common.h 
    src/seq_aligner.h 
    src/dna_seq.h
)
add_executable(
    src/binary_test
    src/binary_test.cpp 
    src/common.h 
    src/dna_seq.h
)
add_executable(
    src/locator
    src/locator.cpp 
    src/common.h 
    src/seq_aligner.h 
    src/dna_seq.h
)
add_executable(
    test/dna_test 
    test/dna_test.cpp
)
add_executable(
    test/aligner_test 
    test/aligner_test.cpp
)
add_executable(
    test/ref_test 
    test/ref_test.cpp
)
target_link_libraries(test/dna_test gtest gtest_main pthread)
target_link_libraries(test/aligner_test gtest gtest_main pthread)
target_link_libraries(test/ref_test gtest gtest_main pthread)
enable_testing()
add_test(
    NAME dna_test
    COMMAND test/dna_test 
)
add_test(
    NAME aligner_test
    COMMAND test/aligner_test
)
add_test(
    NAME ref_test
    COMMAND test/ref_test
)
