cmake_minimum_required (VERSION 2.6)
project (DelphesAnalysis)

set(ENV{DelphesAnaPath} "${CMAKE_SOURCE_DIR}")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
find_package(CLHEPM)
find_package(GSL)

if (CLHEP_FOUND)
  include_directories(${CLHEP_INCLUDE_DIRS})
  link_directories(${CLHEP_LIBRARIES})
endif (CLHEP_FOUND)

if (GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  link_directories(${GSL_LIBRARIES})
endif (GSL_FOUND)

#include_directories(${CMAKE_SOURCE_DIR}/production/delphesCode $ENV{CLHEP_PARAM_PATH}/include ${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_SOURCE_DIR}/production/delphesCode ${CMAKE_SOURCE_DIR}/include)
link_directories(${CMAKE_SOURCE_DIR}/production/delphesCode ${CMAKE_SOURCE_DIR}/lib)


list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS GenVector Gui Core RIO Net Hist Graf Graf3d Gpad Tree Rint Postscript Matrix Physics MathCore Thread PATHS $ENV{ROOTSYS})
#find_package(ROOT)
include(${ROOT_USE_FILE})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


add_subdirectory(src "${CMAKE_SOURCE_DIR}/lib")
add_subdirectory(test)
