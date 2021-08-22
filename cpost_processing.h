/*
OpenCSLSolver: An open-source, parallel, pressure based CFD solver.
Copyright (C) 2021  Bertan ÖZKAN

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
 any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CPOST_PROCESSING_H
#define CPOST_PROCESSING_H

#include "vector"
#include "math.h"
#include "string"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

class Cpost_processing
{
public:
    Cpost_processing();

    std::fstream fstream_VTK;

    void write_VTK_data();
    void write_boundary_VTK_data();

   // void write_CSV_file(std::string filename);
    //void write_boundary_CSV_file(std::string filename);

    void write_VTK_file(std::string filename);
    void write_boundary_VTK_file(std::string filename);

    void monitoring();

};

#include "croot.h"

#endif // CPOST_PROCESSING_H
