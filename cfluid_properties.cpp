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

#include "cfluid_properties.h"

#include "croot.h"

Cfluid_properties::Cfluid_properties()
{

}

void Cfluid_properties::calc_R_specific()
{
    this->R_specific = R_universal / this->molar_mass;
}

void Cfluid_properties::calc_CP_fixed()
{
    this->cp = (this->gamma * this->R_specific)/(this->gamma - 1.0);
}

/*!
 * \brief Function to get density [Kg/M^3]
 * \param Pressure [Pa]
 * \param Temperature [K]
 * \return Density [Kg/M^3]
 */
double Cfluid_properties::get_Density(const double Pressure, const double Temperature)
{
    double rho;
    rho = Pressure/this->R_specific/Temperature;
    return rho;
}

//[J/(Kg*K)]
double Cfluid_properties::get_CP(const double Temperature)
{
    /*double cp;
    cp = this->a_thermo[0] * pow(Temperature,-2.0);
    cp += this->a_thermo[1] * pow(Temperature,-1.0);
    cp += this->a_thermo[2];
    cp += this->a_thermo[3] * Temperature;
    cp += this->a_thermo[4] * pow(Temperature,2.0);
    cp += this->a_thermo[5] * pow(Temperature,3.0);
    cp += this->a_thermo[6] * pow(Temperature,4.0);
    cp = cp * this->R_specific;
    return cp;*/

    return this->cp;
}

double Cfluid_properties::get_gamma( const double Temperature)
{
    /*double cp = this->get_CP(Temperature);
    double gamma = cp/(cp - this->R_specific);
    return gamma;*/

    return this->gamma;
}

double Cfluid_properties::get_Pressure_Rho_derivative(const double Temperature)
{
    double result;
    result = this->R_specific * Temperature;
    double gamma = this->get_gamma(Temperature);
    result = result * gamma;
    return 1.0/result;
}

//[m/s]
double Cfluid_properties::get_speed_of_sound(const double Temperature)
{
    double result;
    result = this->R_specific * Temperature;
    double gamma = this->get_gamma(Temperature);
    result = sqrt(result * gamma);
    return result;
}

//[J/(Kg)]
double Cfluid_properties::get_Enthalpy(const double Temperature, const Svector velocity)
{
    double velmag = mag(velocity);
    double result = H_formation;
    result += (Temperature - Tref) * get_CP(Temperature);
    result += velmag*velmag*0.5;
    return result;
}

// [K]
double Cfluid_properties::get_Temperature(const double init_Temperature , const double Enthalpy, const Svector velocity)
{
    double velmag = mag(velocity);
    double result = init_Temperature;
    /*for (int i = 0 ; i < 7 ; i++) {
        result = Tref + (Enthalpy - H_formation - velmag*velmag*0.5)/get_CP(result);
    }*/
    result = Tref + (Enthalpy - H_formation - velmag*velmag*0.5)/get_CP(result);

    if(result <= 200.0) result = 200.0;
    if(result >= 6000.0) result = 6000.0;
    return result;
}

double Cfluid_properties::get_viscosity(const double Temperature)
{
    double result = Sutherland_coefs[0] * pow(Temperature,3.0/2.0);
    result = result/(Temperature + Sutherland_coefs[1]);
    return result;
}




