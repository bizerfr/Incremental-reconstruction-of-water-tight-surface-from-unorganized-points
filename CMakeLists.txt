cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)

project (SurfaceReconstruction)

find_package(VTK REQUIRED)

include(${VTK_USE_FILE})

add_executable(SurfaceReconstruction MACOSX_BUNDLE 
  DataArray.h
  DataStructure.h 
  GeometricTraits.h 
  NumericComputation.h 
  PointLocator.h 
  SurfaceDelaunay.h
  SurfaceDelaunayUtils.h
  targetver.h
  TypeList.h
  Utility.h
  stdafx.cpp 
  SurfaceReconstruction.cpp)

  
if(VTK_LIBRARIES)
  target_link_libraries(SurfaceReconstruction  ${VTK_LIBRARIES})
else()
  target_link_libraries(SurfaceReconstruction  vtkHybrid vtkWidgets)
endif()


