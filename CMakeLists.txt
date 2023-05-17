cmake_minimum_required(VERSION 3.2)

set(CMAKE_C_COMPILER icc)
set(CMAKE_CXX_COMPILER icpc)
set(CMAKE_Fortran_COMPILER ifort)

project(QE CXX Fortran)

# Options
set(ASSERT on CACHE BOOL "Use assertion (on/off)")
set(INTEGER_64 on CACHE BOOL "Use 64-bit integer (on/off)")
set(PARALLEL "SEQ" CACHE STRING "Parallel mode (SEQ/OPENMP)")

if(NOT ${PARALLEL} MATCHES "SEQ")
    set(USE_PARALLEL on)
else()
    set(USE_PARALLEL off)
endif()

# Compiler Flags
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -wd3058")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -march=native -diag-disable 7713")
add_compile_options(-m64)

# Include
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/external)

# Assert
if(ASSERT)
    add_definitions(-DCHECK_MKL_SPARSE_CALL)
    add_definitions(-DCHECK_CVODE_CALL)
    add_definitions(-DUSE_ASSERT)
    function(define_file_basename_for_sources targetname)
        get_target_property(source_files "${targetname}" SOURCES)
        foreach(sourcefile ${source_files})
            # Add the FILE_BASENAME=filename compile definition to the list.
            get_filename_component(basename "${sourcefile}" NAME)
            # Set the updated compile definitions on the source file.
            set_property(
                SOURCE "${sourcefile}" APPEND
                PROPERTY COMPILE_DEFINITIONS "__FILE_NAME__=\"${basename}\"")
        endforeach()
    endfunction()
endif(ASSERT)

# 64 int
if(INTEGER_64)
    add_definitions(-DUSE_64_INT)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
endif(INTEGER_64)

# MKL
add_definitions(-DUSE_MKL)
add_definitions(-DMKL_Complex8=std::complex<float> -DMKL_Complex16=std::complex<double>)

include_directories($ENV{MKLROOT}/include)
link_directories($ENV{MKLROOT}/lib/intel64)

if(INTEGER_64)
    add_definitions(-DMKL_ILP64)
    set(MKL_64_LIB mkl_intel_ilp64)
else()
    set(MKL_64_LIB mkl_intel_lp64)
endif(INTEGER_64)

if (${PARALLEL} MATCHES "SEQ")
    set(MKL_THREAD_LIB mkl_sequential)
# elseif(${PARALLEL} MATCHES "TBB")
#     set(MKL_THREAD_LIB mkl_tbb_thread)
elseif(${PARALLEL} MATCHES "OPENMP")
    set(MKL_THREAD_LIB mkl_intel_thread)
    set(MKL_THREAD_LIB2 iomp5)
else()
    message(FATAL_ERROR "Unknown parallel mode: ${PARALLEL}")
endif()

set(mkl ${MKL_64_LIB} ${MKL_THREAD_LIB} mkl_core ${MKL_THREAD_LIB2} pthread m dl)

# Blaze
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/external/blaze-cmake)
include(Blaze_Import)

if(INTEGER_64)
    set(BLA_VENDOR_BIT "ilp")
else()
    set(BLA_VENDOR_BIT "lp")
endif(INTEGER_64)

if (${PARALLEL} MATCHES "SEQ")
    set(BLA_VENDOR_PAR "_seq")
elseif(${PARALLEL} MATCHES "OPENMP")
    set(BLA_VENDOR_PAR "")
else()
    message(FATAL_ERROR "Unknown parallel mode: ${PARALLEL}")
endif()

set(BLA_VENDOR "Intel10_64${BLA_VENDOR_BIT}${BLA_VENDOR_PAR}")
message(STATUS "BLA_VENDOR: ${BLA_VENDOR}")

Blaze_Import(
   BLAS on
   BLAS_64BIT ${INTEGER_64}
   BLAS_PARALLEL ${USE_PARALLEL}
   BLAS_MV on
   BLAS_MM on
   BLAS_INCLUDE \"base/blas.h\"
   LAPACK on
   THREADING off
   CACHE_SIZE auto
   VECTORIZATION on
)

# ZVODE
file(GLOB_RECURSE ZVODE_SRC ${CMAKE_SOURCE_DIR}/external/zvode/*.f)
add_library(zvode STATIC ${ZVODE_SRC})

# QE LIB
file(GLOB_RECURSE SRC ${CMAKE_SOURCE_DIR}/src/*.cpp)

add_library(${PROJECT_NAME} STATIC ${SRC})
target_link_libraries(${PROJECT_NAME} Blaze)
if(ASSERT)
    define_file_basename_for_sources(${PROJECT_NAME})
endif(ASSERT)

# EXP
include_directories(${CMAKE_SOURCE_DIR}/exp/include)
file(GLOB_RECURSE EXP_SRC ${CMAKE_SOURCE_DIR}/exp/src/*.cpp)

add_library(exp STATIC ${EXP_SRC})
target_link_libraries(exp Blaze)

function(add_exp exp_name)
    file(GLOB_RECURSE ${exp_name}_SRC ${CMAKE_SOURCE_DIR}/exp/${exp_name}/*.cpp)
    add_executable(${exp_name} ${${exp_name}_SRC})
    # target_compile_options(${exp_name} PUBLIC -g)
    target_link_libraries(${exp_name} exp ${PROJECT_NAME} Blaze zvode ${mkl})
    if(ASSERT)
        define_file_basename_for_sources(${exp_name})
    endif(ASSERT)
endfunction(add_exp)

add_exp(surface)
add_exp(speed)
add_exp(single)

# RUN
file(GLOB_RECURSE RUN_SRC ${CMAKE_SOURCE_DIR}/exp/run/*.cpp)
add_executable(run ${RUN_SRC})