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
MARK_AS_ADVANCED(
  CMAKE_RUNTIME_OUTPUT_DIRECTORY 
  CMAKE_LIBRARY_OUTPUT_DIRECTORY 
  CMAKE_ARCHIVE_OUTPUT_DIRECTORY
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
if (0 AND SWE_ENABLE_CUDA)
  set(GPGPU_BUILD_TYPE SHARED)
  set(SWE_CUDA_SRCS
    src/SWE_BlockCUDA.cu
    src/SWE_BlockCUDA_kernels.cu
    src/SWE_RusanovBlockCUDA.cu
    src/SWE_RusanovBlockCUDA_kernels.cu
    src/SWE_WavePropagationBlockCuda.cu
#    src/SWE_WavePropagationBlockCuda.hh
    src/SWE_WavePropagationBlockCuda_kernels.cu
#    src/SWE_WavePropagationBlockCuda_kernels.hh
  )

  cuda_compile(Module_CUDA_C_SRCS ${SWE_CUDA_SRCS} ${GPGPU_BUILD_TYPE})
endif (0 AND SWE_ENABLE_CUDA)
  
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

