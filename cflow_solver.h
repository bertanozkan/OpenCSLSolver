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

#ifndef CFLOW_SOLVER_H
#define CFLOW_SOLVER_H

#include "cfvm.h"
#include "enums.h"
#include "tensor.h"
#include "memory"

class CFlow_solver
{
public:
    CFlow_solver();

    class CFlow_cell
    {
    public:
        Svector Velocity;
        Svector Velocity_last_t;
        Svector Velocity_grad[3];

        double Density = 0.0;
        double Density_last_t = 0.0;
        Svector Density_grad;

        double Pressure = 0.0;
        double Pressure_last_t = 0.0;
        Svector Pressure_grad;

        double Pressure_corr = 0.0;
        Svector Pressure_corr_grad;

        double Enthalpy = 0.0;
        double Enthalpy_last_t = 0.0;
        Svector Enthalpy_grad;

        double ap = 1.0;
        double delta_t = 0.001;
        Svector ap_vec;
        double Temperature = 0.0;
        double total_mdot = 0.0;
        double Mach = 0.0;
    };


    class CFlow_face
    {
    public:
        Svector stress_explicit_part;
        double Mdot = 0.0;
        double Density = 0.0;
        double Qy = 1.0;
        CFlow_cell* left_cell_ptr = nullptr;
        CFlow_cell* right_cell_ptr = nullptr;
    };


    class CFlow_BC_face
    {
    public:
        double Temperature = 0.0;
        double Pressure = 0.0;
        double Out_Pressure = 0.0;
        double Density = 0.0;
        double Mdot = 0.0;
        double Enthalpy = 0.0;
        double Qy = 1.0;
        Svector Velocity;
        enums::Eflow_BC_type eflow_BC_type;
        enums::Ethermal_BC_type eThermal_BC_type;
        enums::Ewall_BC_type ewall_BC_type;
        CFlow_cell* left_cell_ptr = nullptr;
        bool supersonic_flag = false;
    };


    std::vector<CFlow_cell> cells;
    std::vector<CFlow_face> faces;
    std::vector<CFlow_BC_face> BC_faces;


    // Functions for Momentum equation

    void calc_flux_explicit_visc_part(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    Svector calc_wall_stress(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    template<enums::Eflow_type_visc eflow_type_visc , enums::Econv_flux_type econv_flux_type >
    void calc_flux_momentum(int dimention, const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    template<enums::Eflow_type_visc eflow_type_visc>
    void calc_flux_momentum(const int dimention, const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    template<enums::Etimestepping_type etimestepping_type>
    void cell_loop_momentum(const int dimention, const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    // Support Functions
    template<enums::Econv_flux_type econv_flux_type >
    void calc_face_mdot(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face);

    void calc_face_mdot(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    template<enums::Econv_flux_type econv_flux_type >
    void correct_face_mdot(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void correct_face_mdot(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void calc_density_grad(std::shared_ptr<CFVM> &FVM,  std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void calc_vel_grad(int dimention, std::shared_ptr<CFVM> &FVM,  std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void calc_press_grad(const std::shared_ptr<CFVM> &FVM,  const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void calc_press_corr_grad(const std::shared_ptr<CFVM> &FVM,  const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void calc_enthalpy_grad(std::shared_ptr<CFVM> &FVM,  std::shared_ptr<Cfluid_properties> &Fluid_properties);

    //void calc_timestep(std::shared_ptr<CFVM> &FVM, std::shared_ptr<Cfluid_properties> &Fluid_properties);


    // Functions for Pressure Equation
    template<enums::Econv_flux_type econv_flux_type >
    void calc_flux_pressure(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void calc_flux_pressure(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    template<enums::Etimestepping_type etimestepping_type>
    void cell_loop_pressure(const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


    // Functions for Enthalpy equations
    template<enums::Econv_flux_type econv_flux_type>
    void calc_flux_enthalpy(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void calc_flux_enthalpy(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    template<enums::Etimestepping_type etimestepping_type>
    void cell_loop_enthalpy(const CFVM::CFVM_cell &FVM_cell,CFlow_cell &Flow_cell,
                            std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


    void set_boundary_conditions(const CFVM::CFVM_BC_face &FVM_BC_face , CFlow_BC_face &BC_face, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


    // Residuals
    Svector Momentum_global_residual;
    double Pressure_global_residual = 0.0;
    double Enthalpy_global_residual = 0.0;

    double implicit_underrelaxation = 0.5;
    const double explicit_underrelaxation_P = 0.5;
    const double explicit_underrelaxation_u = 0.8;


    // Solver Functions
    void initialize(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    template<enums::Etimestepping_type etimestepping_type , enums::Eflow_type_visc eflow_type_visc , enums::Econv_flux_type econv_flux_type >
    void iterate(const std::shared_ptr<CFVM> &FVM, std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

    void update_timestep();
};

#endif // CFLOW_SOLVER_H
