cmake_minimum_required(VERSION 2.6)

INCLUDE_DIRECTORIES(${MSCVISUS_SOURCE_DIR}/include)

file(GLOB_RECURSE Incs "${MSCVISUS_SOURCE_DIR}/include/*.h")

ADD_EXECUTABLE(test_indep_set
		${Incs}	
		main_indep_set.cpp)

