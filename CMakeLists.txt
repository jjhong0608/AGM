cmake_minimum_required(VERSION 3.17)
project(AGM)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 17)

include_directories(Header /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/include Header/Eigen)

add_executable(AGM main.cpp Header/ReadFile.cpp Header/ReadFile.h Header/Greenfunction.cpp Header/Greenfunction.h Header/util.cpp Header/util.h Header/Point.cpp Header/Point.h Header/matrixrow.cpp Header/matrixrow.h Header/function.cpp Header/function.h Header/AxialLine.cpp Header/AxialLine.h Header/ReadFile.cpp Header/ReadFile.h Header/matrix.cpp Header/matrix.h Header/WriteFile.cpp Header/WriteFile.h Header/solver.cpp Header/solver.h Header/coordinate.cpp Header/coordinate.h Header/value.cpp Header/value.h Header/GreenfunctionLinear.cpp Header/GreenfunctionLinear.h Header/PointHeat.cpp Header/PointHeat.h Header/PointConvectionDiffusion.cpp Header/PointConvectionDiffusion.h Header/matrixOthers.cpp Header/matrixOthers.h Header/GreenfunctionConvectionDiffusion.cpp Header/GreenfunctionConvectionDiffusion.h Header/GreenfunctionConvectionDiffusionLinear.cpp Header/GreenfunctionConvectionDiffusionLinear.h Header/PointConvectionDiffusionT.cpp Header/PointConvectionDiffusionT.h Header/function2D.cpp Header/function2D.h Header/Eigen Header/matrixEIGEN.cpp Header/matrixEIGEN.h Header/GreenfunctionNeumann.cpp Header/GreenfunctionNeumann.h Header/GreenfunctionReactionDiffusion.cpp Header/GreenfunctionReactionDiffusion.h Header/GreenfunctionReactionDiffusionLinear.cpp Header/GreenfunctionReactionDiffusionLinear.h)

find_library(
        ILP64_MODULE
        NAMES libmkl_intel_lp64.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        THREAD_MODULE
        NAMES libmkl_intel_thread.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        CORE_MODULE
        NAMES libmkl_core.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        SEQUENTIAL_MODULE
        NAMES libmkl_sequential.so
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/mkl/lib/intel64
        REQUIRED
)

find_library(
        IOMP_MODULE
        NAMES libiomp5.a
        HINTS /usr/local/intel/compilers_and_libraries_2018.5.274/linux/compiler/lib/intel64
        REQUIRED
)

find_library(
        LPTHREAD
        NAMES pthread
)

find_library(
        LM
        NAMES m
)

find_library(
        LDL
        NAMES dl
)

target_link_libraries(AGM ${ILP64_MODULE} ${THREAD_MODULE} ${CORE_MODULE} ${SEQUENTIAL_MODULE} ${IOMP_MODULE} ${LPTHREAD} ${LM} ${LDL})

target_compile_options(AGM PUBLIC -xCOMMON-AVX512)
