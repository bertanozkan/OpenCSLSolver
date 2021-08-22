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

#ifndef TESTS_H
#define TESTS_H

#include "cconfiguration.h"

namespace unit_tests {

void run_all_unit_tests()
{
    test_AIR_calculations();
    test_geo_calculations();
    test_matrix_solver_cg();
    test_matrix_solver_ILU_BICGSTAB();
}

}

namespace regression_tests {


Cconfiguration oneD_conduction_test()
{
    Cconfiguration config;

    config.openmp_number_of_threads = 1;
    config.eproblem_type = enums::Eproblem_type::Heat_conduction;
    config.ematrix_solver_type = enums::Ematrix_solver_type::ILU_BICSTAB;
    config.emesh_source = enums::Emesh_source::generate_1D;
    config.geo_parameters_1D[0] = 1.0;
    config.geo_parameters_1D[1] = 0.1;
    config.geo_parameters_1D[2] = 0.1;
    config.geo_parameters_1D[3] = 50;

    config.global_delta_t = 200.0;
    config.thermal_diffusivity = 1.0e-3;
    config.inital_condition->Temperature = 300.0;
    config.set_BC_number(3);
    config.boundary_conditions[0].ethermal_BC_type = enums::Ethermal_BC_type::constant_temperature;
    config.boundary_conditions[0].Temperature = 200.0;
    config.boundary_conditions[1].ethermal_BC_type = enums::Ethermal_BC_type::adiabatic;
    config.boundary_conditions[2].ethermal_BC_type = enums::Ethermal_BC_type::constant_temperature;
    config.boundary_conditions[2].Temperature = 500.0;

    config.max_iteration = 50;

    return config;
}

Cconfiguration shocktube_test()
{
    Cconfiguration config;
    config.openmp_number_of_threads = 1;
    config.eproblem_type = enums::Eproblem_type::Fluid_flow;
    config.eflow_type_visc = enums::Eflow_type_visc::Inviscid;
    config.etimestepping_type = enums::Etimestepping_type::Transient;
    config.ematrix_solver_type = enums::Ematrix_solver_type::CG;
    config.econv_flux_type = enums::Econv_flux_type::Blended_central;
    config.emesh_source = enums::Emesh_source::generate_1D;
    config.geo_parameters_1D[0] = 0.3048;
    config.geo_parameters_1D[1] = 0.03048;
    config.geo_parameters_1D[2] = 0.03048;
    config.geo_parameters_1D[3] = 200;

    config.global_delta_t = 2.11725E-04/200.0;
    config.inital_condition->einitial_condition_type = enums::Einitial_condition_type::Shock_tube;
    config.inital_condition->Temperature = 288.89;
    config.inital_condition->Temperature2 = 231.1;
    config.inital_condition->Pressure = 68947.5;
    config.inital_condition->Pressure2 = 6894.75;
    config.inital_condition->Velocity.val[0] = 0.0;
    config.inital_condition->Velocity.val[1] = 0.0;
    config.inital_condition->Velocity.val[2] = 0.0;
    config.set_BC_number(3);

    config.boundary_conditions[0].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[0].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.boundary_conditions[1].eflow_BC_type = enums::Eflow_BC_type::empty;

    config.boundary_conditions[2].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[2].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.max_iteration = 200;
    config.max_inner_iteration = 10;

    return config;
}


Cconfiguration sub_bump_test()
{
    Cconfiguration config;

    config.openmp_number_of_threads = 4;
    config.eproblem_type = enums::Eproblem_type::Fluid_flow;
    config.eflow_type_visc = enums::Eflow_type_visc::Inviscid;
    config.etimestepping_type = enums::Etimestepping_type::Steady;
    config.ematrix_solver_type = enums::Ematrix_solver_type::ILU_BICSTAB;
    config.econv_flux_type = enums::Econv_flux_type::Blended_central;
    config.emesh_source = enums::Emesh_source::from_file;
    config.mesh_filename = "sub_bump_wedge.su2";

    config.inital_condition->Temperature = 300.0;
    config.inital_condition->Pressure = 101325.0;
    config.inital_condition->Velocity.val[0] = 173.5;
    config.inital_condition->Velocity.val[1] = 0.0;
    config.inital_condition->Velocity.val[2] = 0.0;
    config.set_BC_number(5);

    config.boundary_conditions[0].eflow_BC_type = enums::Eflow_BC_type::empty;

    config.boundary_conditions[1].eflow_BC_type = enums::Eflow_BC_type::velocity_inlet;
    config.boundary_conditions[1].Velocity.val[0] = 173.5;
    config.boundary_conditions[1].Velocity.val[1] = 0.0;
    config.boundary_conditions[1].Velocity.val[2] = 0.0;
    config.boundary_conditions[1].Temperature = 300.0;


    config.boundary_conditions[2].eflow_BC_type = enums::Eflow_BC_type::pressure_outlet;
    config.boundary_conditions[2].Pressure = 101325.0;

    config.boundary_conditions[3].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[3].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.boundary_conditions[4].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[4].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.max_iteration = 6000;
    config.implicit_underrelax = 0.5;
    config.implicit_underrelax_max = 0.85;
    config.underrelax_ramp_iteration = 2000;
    config.VTK_file_frequency = 0;

    return config;

}


Cconfiguration trans_bump_test()
{
    Cconfiguration config;

    config.openmp_number_of_threads = 4;
    config.eproblem_type = enums::Eproblem_type::Fluid_flow;
    config.eflow_type_visc = enums::Eflow_type_visc::Inviscid;
    config.etimestepping_type = enums::Etimestepping_type::Steady;
    config.ematrix_solver_type = enums::Ematrix_solver_type::ILU_BICSTAB;
    config.econv_flux_type = enums::Econv_flux_type::First_order_upwind;
    config.emesh_source = enums::Emesh_source::from_file;
    config.mesh_filename = "bump_mesh.su2";

    config.inital_condition->Temperature = 300.0;
    config.inital_condition->Pressure = 101325.0;
    config.inital_condition->Velocity.val[0] = 234.49;
    config.inital_condition->Velocity.val[1] = 0.0;
    config.inital_condition->Velocity.val[2] = 0.0;
    config.set_BC_number(5);

    config.boundary_conditions[0].eflow_BC_type = enums::Eflow_BC_type::empty;

    config.boundary_conditions[1].eflow_BC_type = enums::Eflow_BC_type::velocity_inlet;
    config.boundary_conditions[1].Velocity.val[0] = 234.49;
    config.boundary_conditions[1].Velocity.val[1] = 0.0;
    config.boundary_conditions[1].Velocity.val[2] = 0.0;
    config.boundary_conditions[1].Temperature = 300.0;


    config.boundary_conditions[2].eflow_BC_type = enums::Eflow_BC_type::pressure_outlet;
    config.boundary_conditions[2].Pressure = 101325.0;

    config.boundary_conditions[3].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[3].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.boundary_conditions[4].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[4].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.max_iteration =15000;
    config.implicit_underrelax = 0.5;
    config.implicit_underrelax_max = 0.8;
    config.underrelax_ramp_iteration = 2000;
    config.VTK_file_frequency = 0;

    return config;

}


Cconfiguration supersonic_bump_test()
{
    Cconfiguration config;

    config.openmp_number_of_threads = 4;
    config.eproblem_type = enums::Eproblem_type::Fluid_flow;
    config.eflow_type_visc = enums::Eflow_type_visc::Inviscid;
    config.etimestepping_type = enums::Etimestepping_type::Steady;
    config.ematrix_solver_type = enums::Ematrix_solver_type::ILU_BICSTAB;
    config.econv_flux_type = enums::Econv_flux_type::First_order_upwind;
    config.emesh_source = enums::Emesh_source::from_file;
    config.mesh_filename = "supersonic_cart.su2";

    config.inital_condition->Temperature = 300.0;
    config.inital_condition->Pressure = 101325.0;
    config.inital_condition->Velocity.val[0] = 573.197;
    config.inital_condition->Velocity.val[1] = 0.0;
    config.inital_condition->Velocity.val[2] = 0.0;
    config.set_BC_number(5);

    config.boundary_conditions[0].eflow_BC_type = enums::Eflow_BC_type::empty;

    config.boundary_conditions[1].eflow_BC_type = enums::Eflow_BC_type::supersonic_velocity_inlet;
    config.boundary_conditions[1].Velocity.val[0] = 573.197;
    config.boundary_conditions[1].Velocity.val[1] = 0.0;
    config.boundary_conditions[1].Velocity.val[2] = 0.0;
    config.boundary_conditions[1].Pressure = 101325.0;
    config.boundary_conditions[1].Temperature = 300.0;


    config.boundary_conditions[2].eflow_BC_type = enums::Eflow_BC_type::pressure_outlet;
    config.boundary_conditions[2].Pressure = 101325.0;

    config.boundary_conditions[3].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[3].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.boundary_conditions[4].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[3].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.max_iteration =15000;
    config.implicit_underrelax = 0.01;
    config.implicit_underrelax_max = 0.2;
    config.underrelax_ramp_iteration = 4000;
    config.VTK_file_frequency = 0;

    return config;

}

Cconfiguration laminar_flatplate_test()
{
    Cconfiguration config;

    config.eproblem_type = enums::Eproblem_type::Fluid_flow;
    config.eflow_type_visc = enums::Eflow_type_visc::Viscous;
    config.etimestepping_type = enums::Etimestepping_type::Steady;
    config.ematrix_solver_type = enums::Ematrix_solver_type::ILU_BICSTAB;
    config.econv_flux_type = enums::Econv_flux_type::Blended_central;
    config.emesh_source = enums::Emesh_source::from_file;
    config.mesh_filename = "lam_flat.su2";

    config.inital_condition->Temperature = 300.0;
    config.inital_condition->Pressure = 97250.0;
    config.inital_condition->Velocity.val[0] = 69.48;
    config.inital_condition->Velocity.val[1] = 0.0;
    config.inital_condition->Velocity.val[2] = 0.0;
    config.set_BC_number(5);

    config.boundary_conditions[0].eflow_BC_type = enums::Eflow_BC_type::empty;

    config.boundary_conditions[1].eflow_BC_type = enums::Eflow_BC_type::velocity_inlet;
    config.boundary_conditions[1].Velocity.val[0] = 69.48;
    config.boundary_conditions[1].Velocity.val[1] = 0.0;
    config.boundary_conditions[1].Velocity.val[2] = 0.0;
    config.boundary_conditions[1].Temperature = 300.0;

    config.boundary_conditions[2].eflow_BC_type = enums::Eflow_BC_type::pressure_outlet;
    config.boundary_conditions[2].Pressure = 97250.0;

    config.boundary_conditions[3].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[3].ewall_BC_type = enums::Ewall_BC_type::slip;

    config.boundary_conditions[4].eflow_BC_type = enums::Eflow_BC_type::wall;
    config.boundary_conditions[4].ethermal_BC_type = enums::Ethermal_BC_type::adiabatic;

    config.max_iteration =15000;
    config.max_inner_iteration = 1;
    config.implicit_underrelax = 0.5;
    config.implicit_underrelax_max = 0.9;
    config.underrelax_ramp_iteration = 1000;
    config.VTK_file_frequency = 0;

    return config;
}





}

#endif // TESTS_H
