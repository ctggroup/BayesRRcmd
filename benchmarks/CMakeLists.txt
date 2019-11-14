include_directories(${EIGEN3_INCLUDE_DIR})

if (ENABLE_MPI)
    add_executable(dataTransfer data_transfer.cpp)

    set_property(TARGET dataTransfer PROPERTY CXX_STANDARD_REQUIRED ON)
    set_property(TARGET dataTransfer PROPERTY CXX_STANDARD 17)

    target_link_libraries(dataTransfer ${MPI_LIBRARIES})
    set_target_properties(dataTransfer PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(dataTransfer PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
