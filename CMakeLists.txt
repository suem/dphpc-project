cmake_minimum_required(VERSION 2.8.12)

project(MaximumMatching)

set(CMAKE_BUILD_TYPE Release)

# Boost
set(Boost_USE_MULTITHREADED ON)

if (WIN32)
	set (BOOST_ROOT "C:/Program Files/boost/boost_1_62_0")
endif()

# if we ever need to use the precompiled libraries, precompile them, put into the lib folder and use:
#if (WIN32)
#	set (BOOST_LIBARYDIR "${BOOST_ROOT}/lib")
#endif()

find_package(Boost 1.62 REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# OpenMP
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
	message("found OpenMP")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
	message("OpenMP NOT FOUND")
endif()

if(WIN32)
	include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
	if(COMPILER_SUPPORTS_CXX11)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
	elseif(COMPILER_SUPPORTS_CXX0X)
	    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
	else()
	        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
	endif()
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -std=c++11")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall")
endif()

if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
  add_definitions (/D "_CRT_SECURE_NO_WARNINGS")
endif()

file(GLOB_RECURSE DPHCP_Files
	"${PROJECT_SOURCE_DIR}/include/*.h"
	"${PROJECT_SOURCE_DIR}/src/*.cpp"
)

set (CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lpthread -ldl")

add_executable(MaximumMatching ${DPHCP_Files})
TARGET_LINK_LIBRARIES(MaximumMatching ${Boost_LIBRARIES} )
