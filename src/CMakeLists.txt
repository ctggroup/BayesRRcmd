project(brr)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(bayes
    analysis.cpp
    analysisgraph.cpp
    analysisrunner.cpp
    BayesRRm.cpp
    BayesW_arms.cpp
    bayesrkernel.cpp
    bayeswbase.cpp
    bayeswkernel.cpp
    data.cpp
    densebayesrkernel.cpp
    densebayeswkernel.cpp
    distributions_boost.cpp
    eigenbayesrkernel.cpp
    gadgets.cpp
    kernel.cpp
    options.cpp
    samplewriter.cpp
    compression.cpp
    limitsequencegraph.cpp
    parallelgraph.cpp
    BayesRBase.cpp
    DenseBayesRRmz.cpp
    densebayesw.cpp
    SparseBayesRRG.cpp
    sequential.cpp
    sparsebayesrkernel.cpp
    sparsebayesw.cpp
    marker.cpp
    markercache.cpp
    densemarker.cpp
    sparsemarker.cpp
    raggedbayesrkernel.cpp
    raggedbayeswkernel.cpp
    raggedsparsemarker.cpp
    eigensparsemarker.cpp
    common.cpp
    markerbuilder.cpp
    densemarkerbuilder.cpp
    raggedsparsemarkerbuilder.cpp
    preprocessgraph.cpp
    eigensparsemarkerbuilder.cpp
    writer.cpp
    logwriter.cpp
    colwriter.cpp
)

set_property(TARGET bayes PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET bayes PROPERTY CXX_STANDARD 17)

target_link_libraries(bayes ${ZLIB_LIBRARIES} ${TBB_IMPORTED_TARGETS}
    ${CMAKE_THREAD_LIBS_INIT} stdc++fs)

add_executable(${PROJECT_NAME}
    main.cpp
)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

target_link_libraries(${PROJECT_NAME} bayes)

set(${PROJECT_NAME}_INCLUDE_DIRS ${PROJECT_SOURCE_DIR} ${PROJECT_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR} CACHE INTERNAL "${PROJECT_NAME}: Include Directories"
    FORCE)
