include_directories("${brr_INCLUDE_DIRS}")
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

# Now simply link against gtest or gtest_main as needed. 

add_executable(tests
    optionstest.cpp
    ppbayestest.cpp
    bayeswtest.cpp
    analysisrunnertest.cpp
    markersubsettest.cpp
)

target_compile_definitions(tests PRIVATE TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../test/data/")
target_compile_definitions(tests PRIVATE GROUPS_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../test.groups/")
target_compile_definitions(tests PRIVATE GAUSS_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/../test.gauss/")
target_compile_definitions(tests PRIVATE TEST_RESULTS="${CMAKE_CURRENT_BINARY_DIR}/results/")
target_compile_definitions(tests PRIVATE WORKING_DIRECTORY="${CMAKE_CURRENT_BINARY_DIR}/working-directory/")

set_property(TARGET tests PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET tests PROPERTY CXX_STANDARD 17)

target_link_libraries(tests ${ZLIB_LIBRARIES})
target_link_libraries(tests gtest_main bayes ${TBB_IMPORTED_TARGETS} )

add_test(NAME example_test COMMAND tests)
