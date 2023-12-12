if (COMP_NAME STREQUAL elm)
  # See Land NaNs in conditionals: https://github.com/E3SM-Project/E3SM/issues/4996
  string(APPEND CMAKE_Fortran_FLAGS " -hfp0")
endif()
# Disable ipa and zero initialization are for other NaN isues:
# https://github.com/E3SM-Project/E3SM/pull/5208
string(APPEND CMAKE_Fortran_FLAGS " -hipa0 -hzero")

string(APPEND CMAKE_C_FLAGS_RELEASE " -O1")
string(APPEND CMAKE_Fortran_FLAGS_RELEASE " -O1")
if (COMP_NAME STREQUAL eam)
  string(APPEND CMAKE_Fortran_FLAGS " -vector0")
endif()
if (compile_threaded)
  string(APPEND CMAKE_CXX_FLAGS " -fopenmp")
endif()
string(APPEND CMAKE_CXX_FLAGS_DEBUG " -g -Wall")
string(APPEND CMAKE_CXX_FLAGS_RELEASE " -O1")

set(E3SM_LINK_WITH_FORTRAN "TRUE")
