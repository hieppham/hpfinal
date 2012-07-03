project(hpfinal CXX C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.4)

# Debug mode switch
SET(DEBUG_MODE 1 CACHE BOOL "Turn off to supress debug outputs")
  SET(DEBUG_FLAG "-g")
IF(DEBUG_MODE)
ENDIF(DEBUG_MODE)
IF(NOT DEBUG_MODE)
	ADD_DEFINITIONS(-DNDEBUG)
ENDIF(NOT DEBUG_MODE)

# Check GAlib availability
find_library(GALIB ga ${GALIB_DIR}/ga)
IF(NOT EXISTS ${GALIB})
  MESSAGE(FATAL_ERROR "Cannot find GAlib library")
ENDIF(NOT EXISTS ${GALIB})

# Check for C++ compiler flags
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG(-std=c++98 HAVE_STD)
IF(HAVE_STD)
	SET(STD_FLAG "-std=c++98")
ENDIF(HAVE_STD)
CHECK_CXX_COMPILER_FLAG("-pedantic -Wall" HAVE_PEDANTIC)
IF(HAVE_PEDANTIC)
	SET(PEDANTIC_FLAG "-pedantic -Wall")
ENDIF(HAVE_PEDANTIC)
CHECK_CXX_COMPILER_FLAG(-O3 HAVE_O3)

# Set C++ compiler flags
SET(CMAKE_CXX_FLAGS "${STD_FLAG} ${PEDANTIC_FLAG} ${DEBUG_FLAG} -I${GALIB_DIR}" CACHE STRING "C++ compiler flags" FORCE)

# Executable binary rrv-visualize
ADD_EXECUTABLE(hpfinal
  hpfinal.cc Customer.cc RouteInfo.cc HGAGenome.cc init.cc BaseOp.cc HPGradProjectGA.cc)
TARGET_LINK_LIBRARIES(hpfinal ${GALIB})

#TARGET_LINK_LIBRARIES(rrv-visualize rrv)
# make install
#INSTALL(TARGETS rrv-compute rrv-visualize DESTINATION bin)
