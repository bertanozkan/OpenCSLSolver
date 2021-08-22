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

#ifdef COMPILE_WITH_OPENMP
    #include <omp.h>
#endif

#include "croot.h"
#include "tests.h"


int main(int argc, char** argv)
{
    //  Reading config file:
    if (argc == 2) {
        root.run_config(argv[1]);
    } else {
        std::cout << "Error. There are no input arguments." << std::endl;
        std::terminate();
    }

    return 0;
}
