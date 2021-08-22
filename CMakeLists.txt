#
#OpenCSLSolver: An open-source, parallel, pressure based CFD solver.
#Copyright (C) 2021  Bertan ÖZKAN
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
# any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

cmake_minimum_required(VERSION 3.5)

project(OpenCSLSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-O3")

add_compile_definitions(COMPILE_WITH_OPENMP)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++ -static-libgcc")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")

add_executable(
        OpenCSLSolver
	cconfiguration.cpp
	cflow_solver.cpp
	cflow_solver_cell_face_loops.cpp
	cflow_solver_supply.cpp
	cfluid_properties.cpp
	cfvm.cpp
	cgeometry.cpp
	cheat_solver.cpp
	cmatrix_solver.cpp
	cpost_processing.cpp
	cpost_processing_monitoring.cpp
	croot.cpp
	main.cpp
)

include_directories(.)


