cmake_minimum_required (VERSION 2.6)
project (DelphesAnalysis)

add_library(DecayChannel SHARED DecayChannel.cpp)

file(GLOB KINFILES *.cc)


add_library(CMSKin SHARED ${KINFILES})
 
