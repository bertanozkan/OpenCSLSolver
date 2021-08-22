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

#ifndef CFLUID_PROPERTIES_SUPPLY_H
#define CFLUID_PROPERTIES_SUPPLY_H

#include "croot.h"
#include "tensor.h"

static inline Cfluid_properties get_AIR_properties()
{
    Cfluid_properties air;
    air.molar_mass = 28.9651159/1000.0;
    /*air.a_thermo[0] = 1.009950160e+04; air.a_thermo[1] = -1.968275610e+02;
    air.a_thermo[2] = 5.009155110e+00; air.a_thermo[3] = -5.761013730e-03;
    air.a_thermo[4] = 1.066859930e-05; air.a_thermo[5] = -7.940297970e-09;
    air.a_thermo[6] = 2.185231910e-12; air.a_thermo[7] = -1.767967310e+02;
    air.a_thermo[8] = -3.921504225e+00;*/


    air.Tref = 298.15;
    air.H_formation = 0.0;

    air.Sutherland_coefs[0] = 1.458e-6;
    air.Sutherland_coefs[1] = 110.4;
    air.Prandlt_number = 0.713;

    air.calc_R_specific();
    air.calc_CP_fixed();

    return air;
}

namespace unit_tests {
static inline void test_AIR_calculations()
{
    Cfluid_properties air;
    air = get_AIR_properties();
    Svector vel;
    std::cout << "300 K air rho = " << air.get_Density(101325.0,300.0) << std::endl;
    std::cout << "500 K air rho = " << air.get_Density(101325.0,500.0) << std::endl;

    std::cout << "300 K air cp = " << air.get_CP(300.0) << std::endl;
    std::cout << "500 K air cp = " << air.get_CP(500.0) << std::endl;

    std::cout << "300 K air gamma = " << air.get_gamma(300.0) << std::endl;
    std::cout << "500 K air gamma = " << air.get_gamma(500.0) << std::endl;

    std::cout << "300 K air Dp/Drho = " << air.get_Pressure_Rho_derivative(300.0) << std::endl;
    std::cout << "500 K air Dp/Drho = " << air.get_Pressure_Rho_derivative(500.0) << std::endl;

    std::cout << "300 K air Speed of Sound = " << air.get_speed_of_sound(300.0) << std::endl;
    std::cout << "500 K air Speed of Sound = " << air.get_speed_of_sound(500.0) << std::endl;

    std::cout << "300 K air Enthalpy = " << air.get_Enthalpy(300.0,vel) << std::endl;
    std::cout << "500 K air Enthalpy = " << air.get_Enthalpy(500.0,vel) << std::endl;

    std::cout << "300 K air Viscosity = " << air.get_viscosity(300.0) << std::endl;
    std::cout << "500 K air Viscosity = " << air.get_viscosity(500.0) << std::endl;
}
}

#endif // CFLUID_PROPERTIES_SUPPLY_H
