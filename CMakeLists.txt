cmake_minimum_required (VERSION 2.6)
project (DelphesAnalysis)

include_directories("${CMAKE_SOURCE_DIR}/production/delphesCode" include)
link_directories("${CMAKE_SOURCE_DIR}/production/delphesCode" ${CMAKE_SOURCE_DIR}/lib)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net PATHS $ENV{ROOTSYS})
include(${ROOT_USE_FILE})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


add_subdirectory(src "${CMAKE_SOURCE_DIR}/lib")
add_subdirectory(test "${CMAKE_SOURCE_DIR}/test")
