PROJECT(SWE)

#--------------------------------------------------------------------------
# Instructions
#
# Notes
# Only use unix syle / path separators, even on windows
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# cmake old/new compatibility
#--------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8)

if (COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0017 NEW)
endif (COMMAND cmake_policy)

#-----------------------------------------------------------------------------
# Setup output Directories
#-----------------------------------------------------------------------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Executables."
)
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all Libraries"
)
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/bin CACHE PATH "Single Directory for all static libraries."
)

#-----------------------------------------------------------------------------
# Build not shared on the Cray
#-----------------------------------------------------------------------------
if(NOT WIN32)
    OPTION(BUILD_NOT_SHARED "Build the program not shared" ON)

    if(BUILD_NOT_SHARED)
        SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)
        SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
        SET(CMAKE_SKIP_RPATH ON)
    endif(BUILD_NOT_SHARED)
endif(NOT WIN32)

#-----------------------------------------------------------------------------
# Hide the advanced variables
#-----------------------------------------------------------------------------
MARK_AS_ADVANCED(
  CMAKE_RUNTIME_OUTPUT_DIRECTORY 
  CMAKE_LIBRARY_OUTPUT_DIRECTORY 
  CMAKE_ARCHIVE_OUTPUT_DIRECTORY
  BUILD_NOT_SHARED
)
#--------------------------------------------------------------------------
# provide a wave model option for user selection 
#--------------------------------------------------------------------------
set(SWE_WaveModel "AugRie" CACHE STRING "WaveModel : one of 'Hybrid FWave AugRie AugRieGeoClaw Rusanov FWavevec'")
set_property(CACHE SWE_WaveModel PROPERTY STRINGS Hybrid FWave AugRie AugRieGeoClaw Rusanov FWavevec)
string(TOLOWER "${SWE_WaveModel}" WaveModel_lower )

if (${WaveModel_lower} STREQUAL "hybrid")
  add_definitions(-D"WAVE_PROPAGATION_SOLVER=0")
elseif (${WaveModel_lower} STREQUAL "fwave")
  add_definitions(-D"WAVE_PROPAGATION_SOLVER=1")
elseif (${WaveModel_lower} STREQUAL "augrie")
  add_definitions(-D"WAVE_PROPAGATION_SOLVER=2")
elseif (${WaveModel_lower} STREQUAL "augriegeoclaw")
  add_definitions(-D"WAVE_PROPAGATION_SOLVER=3")
else()
  message(STATUS "${WaveModel_lower} not yet implemented, using hybrid")
  add_definitions(-D"WAVE_PROPAGATION_SOLVER=0")
endif (${WaveModel_lower} STREQUAL "hybrid")

#--------------------------------------------------------------------------
# Select the output writer
#--------------------------------------------------------------------------
OPTION(ENABLE_NETCDF "Enable NetCDF results writer" OFF)
if(ENABLE_NETCDF)
	add_definitions(-D"WRITENETCDF")
endif(ENABLE_NETCDF)

#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
# provide parallelization options for user selection 
#--------------------------------------------------------------------------
option(SWE_ENABLE_MPI  "Enable MPI Parallel version" OFF)
option(SWE_ENABLE_CUDA "Enable CUDA Parallel version" OFF)

if (SWE_ENABLE_MPI)
 find_package(MPI REQUIRED)
endif (SWE_ENABLE_MPI)

if (SWE_ENABLE_CUDA)
 find_package(CUDA REQUIRED)
endif (SWE_ENABLE_CUDA)

#--------------------------------------------------------------------------
# include directories
#
# any includes listed here apply to this project and 
# any projects in sub dirs
#--------------------------------------------------------------------------
include_directories(
  src
  submodules/swe_solvers/src
)

set(SWE_SRCS 
  src/examples/swe_wavepropagation.cpp
  src/SWE_Block.cpp
  src/SWE_RusanovBlock.cpp
  src/SWE_WavePropagationBlock.cpp
)
if(ENABLE_NETCDF)
add_executable(swe_simple
  src/blocks/SWE_Block.cpp
  src/blocks/SWE_WavePropagationBlock.cpp
  src/writer/NetCdfWriter.cpp
  src/tools/Logger.cpp
  src/examples/swe_simple.cpp
)
else(ENABLE_NETCDF)

#--------------------------------------------------------------------------
# simple non mpi example
#--------------------------------------------------------------------------

add_executable(swe_simple
  src/blocks/SWE_Block.cpp
  src/blocks/SWE_WavePropagationBlock.cpp
  src/writer/VtkWriter.cpp
  src/tools/Logger.cpp
  src/examples/swe_simple.cpp
)
endif(ENABLE_NETCDF)

#--------------------------------------------------------------------------
# mpi  example
#--------------------------------------------------------------------------
if (SWE_ENABLE_MPI)
  include_directories(${MPI_INCLUDE_PATH})
  
  add_executable(swe_mpi
    src/blocks/SWE_Block.cpp
    src/blocks/SWE_WavePropagationBlock.cpp
    src/writer/VtkWriter.cpp
    src/tools/Logger.cpp
    src/examples/swe_mpi.cpp
  )

  target_link_libraries(swe_mpi ${MPI_LIBRARY})

endif (SWE_ENABLE_MPI)
  
#--------------------------------------------------------------------------
# cuda example
#--------------------------------------------------------------------------
if (SWE_ENABLE_CUDA AND SWE_ENABLE_MPI)

  option(CUDA_VERBOSE_PTXAS
    "On to enable verbose output from the PTXAS assembler."
    OFF
  )

  set(GPGPU_BUILD_TYPE STATIC)
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  set(BUILD_SHARED_LIBS OFF)
  set(CUDA_SEPARABLE_COMPILATION ON)
  list(APPEND CUDA_NVCC_FLAGS "-arch=sm_21")
#  list(APPEND CUDA_NVCC_FLAGS "rdc=false")

  set(SWE_CUDA_SRCS
    src/blocks/cuda/SWE_BlockCUDA.cu
    src/blocks/cuda/SWE_BlockCUDA_kernels.cu
#    src/blocks/SWE_RusanovBlockCUDA.cu
#    src/blocks/SWE_RusanovBlockCUDA_kernels.cu
    src/blocks/cuda/SWE_WavePropagationBlockCuda.cu
    src/blocks/cuda/SWE_WavePropagationBlockCuda_kernels.cu
#    src/SWE_WavePropagationBlockCuda.hh
#    src/SWE_WavePropagationBlockCuda_kernels.hh
  )


  set(GENCODE_SM10 -gencode=arch=compute_10,code=sm_10 -gencode=arch=compute_10,code=compute_10)
  set(GENCODE_SM13 -gencode=arch=compute_13,code=sm_13 -gencode=arch=compute_13,code=compute_13)
  set(GENCODE_SM20 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20)


  if (CUDA_VERBOSE_PTXAS)
    set(VERBOSE_PTXAS --ptxas-options=-v)
  endif (CUDA_VERBOSE_PTXAS)

  cuda_add_library(gpuacceleration ${LIB_TYPE}
    ${CCFILES}
    ${HFILES}
    ${CUHFILES}
    ${HFILES_PUBLIC}
    ${SWE_CUDA_SRCS}
    #OPTIONS ${GENCODE_SM20} ${VERBOSE_PTXAS}
    OPTIONS ${GENCODE_SM20} ${VERBOSE_PTXAS}
  )

#  set(LIB_NAME "gpuacceleration")
#  cuda_add_library(${LIB_NAME}
#    ${SWE_CUDA_SRCS} 
 #   OPTIONS -DSTUFF="blah blah"
 #   RELEASE -DNDEBUG
 #   DEBUG -g -DDEBUG
 # )

  source_group("CUDA Source Files" FILES ${SWE_CUDA_SRCS})

  cuda_compile(Module_CUDA_C_SRCS ${SWE_CUDA_SRCS} ${GPGPU_BUILD_TYPE})

  add_executable(swe_mpi_cuda
    src/blocks/SWE_Block.cpp
    src/blocks/SWE_WavePropagationBlock.cpp
    src/writer/VtkWriter.cpp
    src/tools/Logger.cpp
    src/examples/swe_mpi.cpp
  )

  set_property(TARGET swe_mpi_cuda APPEND PROPERTY COMPILE_DEFINITIONS CUDA)

  target_link_libraries(swe_mpi_cuda 

    gpuacceleration 
    ${MPI_LIBRARY}
    ${CUDA_LIBRARIES}
  )




endif (SWE_ENABLE_CUDA AND SWE_ENABLE_MPI)
  
#set(SWE_EXAMPLE_SIMPLE
#  src/examples/swe_opengl.cpp
#)

set(SWE_HDRS
  src/scenarios/SWE_AsagiScenario.hpp
  src/scenarios/SWE_Scenario.h
  src/scenarios/SWE_simple_scenarios.h
  src/solvers/augrie.hpp
  src/solvers/augriegeoclaw.hpp
  src/solvers/fwave.hpp
  src/solvers/FWaveCuda.h
  src/solvers/hybrid.hpp
  src/solvers/wavepropagation.hpp
  src/tools/Logger.hpp
)

#      <CodeGeneration Condition="'$(Configuration)|$(Platform)'=='CUDA|Win32'">compute_20,sm_20</CodeGeneration>

set(SWE_CUDA_HDRS
  src/swe_block.hh
  src/swe_blockcuda.hh
  src/SWE_BlockCUDA_kernels.hh
  src/swe_rusanovblock.hh
  src/SWE_RusanovBlockCUDA.hh
  src/SWE_RusanovBlockCUDA_kernels.hh
  src/swe_wavepropagationblock.hh
)

#add_executable(SWE 
#  ${SWE_SRCS}
#)

