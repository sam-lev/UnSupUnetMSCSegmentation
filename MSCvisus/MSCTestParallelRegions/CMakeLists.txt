cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${MSCVISUS_SOURCE_DIR}/include)

ADD_EXECUTABLE(test_parallel_regions 
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

