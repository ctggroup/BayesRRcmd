cmake_minimum_required(VERSION 2.8.12)

option(BRR_TRACE_FLOW_GRAPH "Enable flow graph tracing" OFF)

if(BRR_TRACE_FLOW_GRAPH)
    add_definitions(-DTBB_USE_THREADING_TOOLS)
    add_definitions(-DTBB_PREVIEW_FLOW_GRAPH_TRACE)
    add_definitions(-DTBB_PREVIEW_FLOW_GRAPH_FEATURES)
    add_definitions(-DTBB_PREVIEW_ALGORITHM_TRACE)
endif()

add_subdirectory(src)

if(TEST)
	add_subdirectory(test)
endif()
