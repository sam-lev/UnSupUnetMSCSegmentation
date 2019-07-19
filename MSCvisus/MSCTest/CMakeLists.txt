cmake_minimum_required(VERSION 2.6)

set (CMAKE_BUILD_TYPE Release)

INCLUDE_DIRECTORIES(${MSCVISUS_SOURCE_DIR}/include)

ADD_DEFINITIONS(-DPROBABILITY_CUTOFF_VALUE=0.05)

ADD_EXECUTABLE(test 
		${MSCVISUS_SOURCE_DIR}/include/mscIndexTypes.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicIterator.h		
		${MSCVISUS_SOURCE_DIR}/include/mscBasicArray.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicDataHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicMeshHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicMeshFunction.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicGradientField.h
		${MSCVISUS_SOURCE_DIR}/include/mscDumbGradientField.h
		${MSCVISUS_SOURCE_DIR}/include/mscRegularArray.h		
		
		main.cpp)

