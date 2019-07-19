cmake_minimum_required(VERSION 2.6)

set (CMAKE_BUILD_TYPE Release)

INCLUDE_DIRECTORIES(${MSCVISUS_SOURCE_DIR}/include)

ADD_DEFINITIONS(-DPROBABILITY_CUTOFF_VALUE=0.05)

file(GLOB_RECURSE Incs "${MSCVISUS_SOURCE_DIR}/include/*.h")

ADD_EXECUTABLE(test_block_streaming 
		${Incs}		
		
		main_block_streaming.cpp)

