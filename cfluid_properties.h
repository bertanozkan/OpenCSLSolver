/*
 * OpenCSLSolver: An open-source, parallel, pressure based CFD solver.
 * Copyright (C) 2021  Bertan ÖZKAN
 *
 * Email: ozkanbertan@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CFLUID_PROPERTIES_H
#define CFLUID_PROPERTIES_H

#include "tensor.h"

/*!
 * \brief Universal Gas constant [J/(K*mol)]
 */
const double R_universal = 8.31446261815324;


class Cfluid_properties
{
public:
    Cfluid_properties();

    //double a_thermo[9];
    double cp = 0.0;
    double gamma = 1.4;
    /*!
     * \brief Molar Mass [Kg/mol]
     */
    double molar_mass = 28.9651159/1000.0;
    // [K]
    double Tref = 298.15;
    //[J/(Kg)]
    double H_formation = 0.0;

    double Sutherland_coefs[2] = {1.458e-6, 110.4};
    double Prandlt_number = 0.713;

    /*!
     * \brief Specific Gas constant [J/(Kg*K)]
     */
    double R_specific;
    void calc_CP_fixed();

    void calc_R_specific();
    double get_Density(const double Pressure, const double Temperature);
    double get_CP( const double Temperature);
    double get_gamma( const double Temperature);
    double get_Pressure_Rho_derivative(const double Temperature);
    double get_speed_of_sound(const double Temperature);
    double get_Enthalpy(const double Temperature , const Svector velocity);
    double get_Temperature(const double init_Temperature , const double Enthalpy, const Svector velocity);
    double get_viscosity(const double Temperature);

};



#endif // CFLUID_PROPERTIES_H
