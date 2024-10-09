#####################################################################################################################
#
# CMake and System Settings
#
#####################################################################################################################

cmake_minimum_required(VERSION 2.8)

if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum OS X deployment version")
endif()

if (UNIX AND NOT CMAKE_C_COMPILER)
	set(CMAKE_C_COMPILER gcc)
	set(CMAKE_CXX_COMPILER g++)
endif()
set(CMAKE_CXX_STANDARD 11)


#####################################################################################################################
#
# Project Settings 
#
#####################################################################################################################

Project(solarpilot_codelibs)

if(NOT BUILD_COPILOT_API)
	add_subdirectory(splinter)
endif()
add_subdirectory(shared)
add_subdirectory(solarpilot)
add_subdirectory(nlopt)
add_subdirectory(lpsolve)
