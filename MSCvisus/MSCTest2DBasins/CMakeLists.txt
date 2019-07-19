cmake_minimum_required(VERSION 2.6)

set (CMAKE_BUILD_TYPE Release)

INCLUDE_DIRECTORIES(${MSCVISUS_SOURCE_DIR}/include
   ${MSCVISUS_SOURCE_DIR}/MSCTest2DBasins/include
)

ADD_DEFINITIONS(-DPROBABILITY_CUTOFF_VALUE=0.05)

ADD_EXECUTABLE(test_2d_basins

		${MSCVISUS_SOURCE_DIR}/include/mscIndexTypes.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicIterator.h		
		${MSCVISUS_SOURCE_DIR}/include/mscBasicArray.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicDataHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicMeshHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicMeshFunction.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicGradientField.h
		${MSCVISUS_SOURCE_DIR}/include/mscDumbGradientField.h
		${MSCVISUS_SOURCE_DIR}/include/mscRegularArray.h		
		${MSCVISUS_SOURCE_DIR}/include/mscArrayFactory.h
		${MSCVISUS_SOURCE_DIR}/include/mscRegularRawDataHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscRegularGrid3DImplicitMeshHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscRegularGrid2DImplicitMeshHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscDumbStoringMeshFunction.h
		${MSCVISUS_SOURCE_DIR}/include/mscDumbStoringMinFunction.h
		${MSCVISUS_SOURCE_DIR}/include/mscSimpleGradientBuilder.h
		${MSCVISUS_SOURCE_DIR}/include/mscSimpleRandomGradientBuilder.h
		${MSCVISUS_SOURCE_DIR}/include/mscSimpleGradientUsingAlgorithms.h
		${MSCVISUS_SOURCE_DIR}/include/mscRegularGrid3DGradientField.h
		${MSCVISUS_SOURCE_DIR}/include/mscRegularGrid3DMeshFunction.h
		${MSCVISUS_SOURCE_DIR}/include/mscTwoWay3DGradientBuilder.h
		${MSCVISUS_SOURCE_DIR}/include/mscConvergentGradientBuilder.h
		${MSCVISUS_SOURCE_DIR}/include/mscNegatingMeshFunction.h
		${MSCVISUS_SOURCE_DIR}/include/mscComplementMeshHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscModifiedBoundaryMeshHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscCombinationGradientField.h
		${MSCVISUS_SOURCE_DIR}/include/mscBasicMSC.h
		${MSCVISUS_SOURCE_DIR}/include/mscNegatingDataHandler.h
		${MSCVISUS_SOURCE_DIR}/include/mscTopoArray.h
		${MSCVISUS_SOURCE_DIR}/include/mscPreClassifier.h
		${MSCVISUS_SOURCE_DIR}/include/mscSimpleConstrainedGradientBuilder.h
		${MSCVISUS_SOURCE_DIR}/include/mscAssistedGradientBuilder.h

		main_2d_basins.cpp)

