#BOOST
set(Boost_NO_SYSTEM_PATHS TRUE)
if (Boost_NO_SYSTEM_PATHS)
  set(BOOST_ROOT "/XF/Libraries/Boost/1.62.0/Installed-gcc5.4.0")
  set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
  set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
endif (Boost_NO_SYSTEM_PATHS)
FIND_PACKAGE(Boost)
IF(Boost_FOUND)
	INCLUDE_DIRECTORIES ( ${Boost_INCLUDE_DIR} )
	MESSAGE(STATUS "...Boost found : " ${Boost_INCLUDE_DIR})
ELSE(Boost_FOUND)
	MESSAGE(FATAL_ERROR "Cannot build Boost project without Boost. Please set Boost_DIR.")
ENDIF(Boost_FOUND)
