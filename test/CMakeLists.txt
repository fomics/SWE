#--------------------------------------------------------------------------
# Simple test to get some practice
#--------------------------------------------------------------------------
add_executable(test_simple test.c)
add_test(
  NAME SimpleTest
  COMMAND $<TARGET_FILE:test_simple> 
)

#--------------------------------------------------------------------------
# complex test to try netCDF check
#--------------------------------------------------------------------------
add_executable(test_netcdf netcdf_test.cpp)
target_link_libraries(test_netcdf ${NETCDF_LIBRARIES})

set(X_SIZE 64)
set(Y_SIZE 64)

add_test(NAME NetCDF_IO_TEST
    COMMAND ${CMAKE_COMMAND}
            -DCMD1=$<TARGET_FILE:swe_simple>
            -DARG1=--grid-size-x=${X_SIZE} 
            -DARG2=--grid-size-y=${Y_SIZE}
            -DARG3=--output-basepath=${CMAKE_CURRENT_BINARY_DIR}/temp

            -DCMD2=$<TARGET_FILE:test_netcdf>
            -DARG4=${CMAKE_CURRENT_BINARY_DIR}/temp_00.nc
            -DARG5=${X_SIZE} 
            -DARG6=${Y_SIZE} 
    -P ${CMAKE_CURRENT_SOURCE_DIR}/run_test.cmake
)


