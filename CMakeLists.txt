cmake_minimum_required (VERSION 2.6)
project (DelphesAnalysis)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/production/delphesCode")
link_directories("${CMAKE_CURRENT_SOURCE_DIR}/production/delphesCode")

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net PATHS $ENV{ROOTSYS})
include(${ROOT_USE_FILE})

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")


add_executable(JetTestAna JetTestAna.cpp)
add_executable(PrintParm  PrintParm.cpp)
add_library(DecayChannel  DecayChannel.cpp)
target_link_libraries(JetTestAna Delphes ${ROOT_LIBRARIES})
target_link_libraries(PrintParm Delphes ${ROOT_LIBRARIES} DecayChannel)
#target_link_libraries(decayChannel Delphes ${ROOT_LIBRARIES})
