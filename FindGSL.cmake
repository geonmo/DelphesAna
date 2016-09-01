# - Try to find the GSL library.
#
# The following are set after configuration is done: 
#  GSL_FOUND
#  GSL_INCLUDE_DIRS
#  GSL_LIBRARY_DIRS
#  GSL_LIBRARIES

#find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_multimin.h PATHS $ENV{ROOT_INCLUDE_PATH})
#find_library(GSL_LIBRARY NAMES gsl PATHS $ENV{ROOT_INCLUDE_PATH})
find_path(GSL_INCLUDE_DIR NAMES gsl/gsl_multimin.h PATHS /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/1.16/include)
find_library(GSL_LIBRARY NAMES gsl PATHS /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/gsl/1.16/lib)

message("GSL include dir = ${GSL_INCLUDE_DIR}")
message("GSL lib = ${GSL_LIBRARY}")

set(GSL_LIBRARIES ${GSL_LIBRARY})
set(GSL_INCLUDE_DIRS ${GSL_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
# Handle the QUIETLY and REQUIRED arguments and set the GSL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(GSL DEFAULT_MSG
                                  GSL_LIBRARY GSL_INCLUDE_DIR)

mark_as_advanced(GSL_INCLUDE_DIR GSL_LIBRARY)
