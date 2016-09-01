# - Try to find the CLHEP library.
#
# The following are set after configuration is done: 
#  CLHEP_FOUND
#  CLHEP_INCLUDE_DIRS
#  CLHEP_LIBRARY_DIRS
#  CLHEP_LIBRARIES

find_path(CLHEP_INCLUDE_DIR NAMES CLHEP/Random/RandomEngine.h PATHS /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/clhep/2.2.0.4-giojec/include)
find_library(CLHEP_LIBRARY NAMES CLHEP PATHS /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/clhep/2.2.0.4-giojec/lib)

message("CLHEP include dir = ${CLHEP_INCLUDE_DIR}")
message("CLHEP lib = ${CLHEP_LIBRARY}")

set(CLHEP_LIBRARIES ${CLHEP_LIBRARY})
set(CLHEP_INCLUDE_DIRS ${CLHEP_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set the CLHEP_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CLHEP DEFAULT_MSG
                                  CLHEP_LIBRARY CLHEP_INCLUDE_DIR)

mark_as_advanced(CLHEP_INCLUDE_DIR CLHEP_LIBRARY)
