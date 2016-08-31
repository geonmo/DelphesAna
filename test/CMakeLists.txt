cmake_minimum_required (VERSION 2.6)
project (DelphesAnalysis)

add_executable(JetTestAna JetTestAna.cpp)
add_executable(PrintParm  PrintParm.cpp)
target_link_libraries(JetTestAna Delphes ${ROOT_LIBRARIES})
target_link_libraries(PrintParm Delphes ${ROOT_LIBRARIES} DecayChannel)
#target_link_libraries(decayChannel Delphes ${ROOT_LIBRARIES})
