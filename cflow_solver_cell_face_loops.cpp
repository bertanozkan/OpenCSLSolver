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

#include "cflow_solver.h"
#include "cfluid_properties.h"
#include "tensor.h"

void CFlow_solver::calc_flux_explicit_visc_part(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face ,const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    Stensor dum;
    dum = gradient(Flow_face.left_cell_ptr->Velocity,Flow_face.right_cell_ptr->Velocity,FVM_face.cell_center_unit_vector, FVM_face.cell_center_distance);
    dum = transpose(dum);
    dum = dum - ((2.0/3.0) * divergence_I(Flow_face.left_cell_ptr->Velocity,Flow_face.right_cell_ptr->Velocity, FVM_face.cell_center_unit_vector,FVM_face.cell_center_distance));

    double T_face = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Temperature,Flow_face.right_cell_ptr->Temperature);
    dum = Fluid_properties->get_viscosity(T_face) * dum;
    Flow_face.stress_explicit_part = (dum * FVM_face.area_vector);
}


template<enums::Eflow_type_visc eflow_type_visc , enums::Econv_flux_type econv_flux_type >
void CFlow_solver::calc_flux_momentum(const int dimention, const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    solvers_common::Convective_flux<econv_flux_type>
            ( Matrix_solver,FVM_face, Flow_face.Mdot,Flow_face.left_cell_ptr->Velocity[dimention],Flow_face.right_cell_ptr->Velocity[dimention]);

    if constexpr(eflow_type_visc == enums::Eflow_type_visc::Viscous)
    {
        double T_face = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Temperature,Flow_face.right_cell_ptr->Temperature);
        double D =(Fluid_properties->get_viscosity(T_face));
        D = D * FVM_face.diffusion_ksi_A;
        solvers_common::Diffusive_flux(Matrix_solver,FVM_face, D);

        // Adding explicit part to left side
        Matrix_solver->b_vector.add_insert_element(FVM_face.left_cell_index,  Flow_face.stress_explicit_part[dimention]);

        // Adding explicit part to right side
        Matrix_solver->b_vector.add_insert_element(FVM_face.right_cell_index, - Flow_face.stress_explicit_part[dimention]);
    }
}

template void CFlow_solver::calc_flux_momentum<enums::Eflow_type_visc::Inviscid , enums::Econv_flux_type::First_order_upwind >
(const int dimention, const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::calc_flux_momentum<enums::Eflow_type_visc::Inviscid , enums::Econv_flux_type::Blended_central >
(const int dimention, const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::calc_flux_momentum<enums::Eflow_type_visc::Viscous , enums::Econv_flux_type::First_order_upwind >
(const int dimention, const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::calc_flux_momentum<enums::Eflow_type_visc::Viscous , enums::Econv_flux_type::Blended_central >
(const int dimention, const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


template<enums::Eflow_type_visc eflow_type_visc>
void CFlow_solver::calc_flux_momentum(const int dimention, const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    // If empty BC, do nothing
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::empty){

    }

    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::symmetry){

    }

    // If boundary wall
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::wall){
        if constexpr( eflow_type_visc == enums::Eflow_type_visc::Viscous)
        {
            if( Flow_BC_face.ewall_BC_type == enums::Ewall_BC_type::no_slip)
            {
                // No convective flux
                Svector norm = (1.0/FVM_BC_face.area) * FVM_BC_face.area_vector;

                double D = (Fluid_properties->get_viscosity(Flow_BC_face.left_cell_ptr->Temperature));
                D *= FVM_BC_face.area / FVM_BC_face.perp_dist;
                double dum = D * (1.0 - (norm[dimention] * norm[dimention]));
                Matrix_solver->matrix.add_insert_element(FVM_BC_face.left_cell_index,FVM_BC_face.left_cell_index,dum);

                dum = 0.0;
                for(int i = dimention + 1; i < dimention + 3 ; i++ )
                {
                    int i_ = i % 3;
                    dum += Flow_BC_face.left_cell_ptr->Velocity[i_] * norm[i_] * norm[dimention];
                }

                dum = D * dum;
                Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,dum);
            }
        }
    }

    // If boundary velocity inlet
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::velocity_inlet || Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::supersonic_velocity_inlet){

        //Calc RHS for left cell
        double value = Flow_BC_face.Mdot *Flow_BC_face.Velocity.val[dimention];

        // Fill RHS for the left cell
        Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,-value);
    }

    // If boundary is pressure outlet
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::pressure_outlet){

        Matrix_solver->matrix.add_insert_element(FVM_BC_face.left_cell_index,FVM_BC_face.left_cell_index,Flow_BC_face.Mdot);
    }
}

template void CFlow_solver::calc_flux_momentum<enums::Eflow_type_visc::Inviscid>
    (const int dimention, const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::calc_flux_momentum<enums::Eflow_type_visc::Viscous>
    (const int dimention, const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


template<enums::Etimestepping_type etimestepping_type>
void CFlow_solver::cell_loop_momentum(const int dimention, const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    double dum = 0.0;
    double dum2 = 0.0;
    if constexpr(etimestepping_type == enums::Etimestepping_type::Transient)
    {
        dum = (FVM_cell.volume / Flow_cell.delta_t);
        dum2 = dum*Flow_cell.Density;
        Matrix_solver->matrix.add_insert_element(FVM_cell.index,FVM_cell.index,dum2);
        dum2 = dum*(Flow_cell.Velocity_last_t.val[dimention]*Flow_cell.Density_last_t);
    }
    dum2 = dum2 - Flow_cell.Pressure_grad.val[dimention]*FVM_cell.volume;
    Matrix_solver->b_vector.add_insert_element(FVM_cell.index,dum2);



    if(etimestepping_type == enums::Etimestepping_type::Steady)
    {
        dum = Matrix_solver->matrix.get_element(FVM_cell.index,FVM_cell.index);
        Matrix_solver->matrix.insert_element(FVM_cell.index,FVM_cell.index,dum/this->implicit_underrelaxation);
        dum = (1.0-this->implicit_underrelaxation)/this->implicit_underrelaxation*dum*Flow_cell.Velocity_last_t.val[dimention];
        Matrix_solver->b_vector.add_insert_element(FVM_cell.index,dum);
    }

    Flow_cell.ap += Matrix_solver->matrix.get_element(FVM_cell.index,FVM_cell.index)/3.0;
    Flow_cell.ap_vec.val[dimention] = Matrix_solver->matrix.get_element(FVM_cell.index,FVM_cell.index);

    Matrix_solver->x_vector.insert_element(FVM_cell.index , Flow_cell.Velocity[dimention]);

    Flow_cell.total_mdot = 0.0;
}

template void CFlow_solver::cell_loop_momentum<enums::Etimestepping_type::Steady>
(const int dimention, const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::cell_loop_momentum<enums::Etimestepping_type::Transient>
(const int dimention, const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                                      std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


template<enums::Econv_flux_type econv_flux_type >
void CFlow_solver::calc_flux_pressure(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    // insert uncorrected Mass flows to RHS
    Matrix_solver->b_vector.add_insert_element(FVM_face.left_cell_index,-Flow_face.Mdot);
    Flow_face.left_cell_ptr->total_mdot += Flow_face.Mdot;
    Matrix_solver->b_vector.add_insert_element(FVM_face.right_cell_index,Flow_face.Mdot);
    Flow_face.right_cell_ptr->total_mdot -= Flow_face.Mdot;



    // Calculate and apply convective part of Pressure equation
    double dum = 0.0;


    //double T_face = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Temperature,Flow_face.right_cell_ptr->Temperature);

    double T_face = solvers_common::Convective_face_interp<econv_flux_type>
            (FVM_face, Flow_face.Mdot, Flow_face.left_cell_ptr->Temperature,Flow_face.right_cell_ptr->Temperature);

    dum = Fluid_properties->get_Pressure_Rho_derivative(T_face);
    dum = dum / Flow_face.Density;

    dum *= Flow_face.Mdot;
    solvers_common::Convective_flux<econv_flux_type>
            (Matrix_solver,FVM_face, dum, Flow_face.left_cell_ptr->Pressure_corr, Flow_face.right_cell_ptr->Pressure_corr);

    double avgD = Flow_face.Qy*Flow_face.Density;

    solvers_common::Diffusive_flux( Matrix_solver, FVM_face,avgD);

}

template void CFlow_solver::calc_flux_pressure<enums::Econv_flux_type::First_order_upwind >
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::calc_flux_pressure<enums::Econv_flux_type::Blended_central >
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


void CFlow_solver::calc_flux_pressure(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{

    // If symmetry BC, do nothing
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::empty){

    }

    // If symmetry BC, do nothing
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::symmetry){

    }

    // If boundary wall, do nothing
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::wall){

    }

    // If boundary velocity inlet
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::velocity_inlet){
        // insert uncorrected Mass flows to RHS
        Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,-Flow_BC_face.Mdot);
        Flow_BC_face.left_cell_ptr->total_mdot += Flow_BC_face.Mdot;

        // Calculate and apply convective part of Pressure equation
        double dum = 0.0;

        dum = Fluid_properties->get_Pressure_Rho_derivative(Flow_BC_face.Temperature)
                    / Flow_BC_face.left_cell_ptr->Density;
        dum *= Flow_BC_face.Mdot;

        solvers_common::Convective_flux_neumann_zero_BC(Matrix_solver,FVM_BC_face, dum);

    }

    // If boundary supersonic velocity inlet
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::supersonic_velocity_inlet){
        // insert uncorrected Mass flows to RHS
        Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,-Flow_BC_face.Mdot);
        Flow_BC_face.left_cell_ptr->total_mdot += Flow_BC_face.Mdot;

        double avgD = Flow_BC_face.Qy*Flow_BC_face.Density;
        solvers_common::Diffusive_flux_dirichlet_BC( Matrix_solver, FVM_BC_face, avgD,0.0);
    }

    // If boundary is pressure outlet
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::pressure_outlet){
        // insert uncorrected Mass flows to RHS
        Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,-Flow_BC_face.Mdot);
        Flow_BC_face.left_cell_ptr->total_mdot += Flow_BC_face.Mdot;

        if(Flow_BC_face.supersonic_flag)
        {
            double dum = 0.0;

            dum = Fluid_properties->get_Pressure_Rho_derivative(Flow_BC_face.left_cell_ptr->Temperature)
                        / Flow_BC_face.left_cell_ptr->Density;
            dum *= Flow_BC_face.Mdot;

            solvers_common::Convective_flux_neumann_zero_BC(Matrix_solver,FVM_BC_face, dum);
        }
        else
        {
            double avgD = Flow_BC_face.Qy*Flow_BC_face.left_cell_ptr->Density;
            solvers_common::Diffusive_flux_dirichlet_BC( Matrix_solver, FVM_BC_face, avgD,0.0);
        }
    }
}

template<enums::Etimestepping_type etimestepping_type>
void CFlow_solver::cell_loop_pressure(const CFVM::CFVM_cell &FVM_cell,CFlow_cell &Flow_cell,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    if constexpr(etimestepping_type == enums::Etimestepping_type::Transient)
    {
        double dum = (FVM_cell.volume / Flow_cell.delta_t);

        double dum2 = dum * (Flow_cell.Density - Flow_cell.Density_last_t);
        Matrix_solver->b_vector.add_insert_element(FVM_cell.index,-dum2);

        dum2 = dum * Fluid_properties->get_Pressure_Rho_derivative(Flow_cell.Temperature);
        Matrix_solver->matrix.add_insert_element(FVM_cell.index,FVM_cell.index,dum2);
    }

    if constexpr(etimestepping_type == enums::Etimestepping_type::Steady)
    {

        Matrix_solver->matrix.insert_element(FVM_cell.index,FVM_cell.index,
                                             Matrix_solver->matrix.get_element(FVM_cell.index,FVM_cell.index)/this->implicit_underrelaxation);
    }

    Matrix_solver->x_vector.insert_element(FVM_cell.index , Flow_cell.Pressure_corr);
}

template void CFlow_solver::cell_loop_pressure<enums::Etimestepping_type::Steady>
(const CFVM::CFVM_cell &FVM_cell,CFlow_cell &Flow_cell,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::cell_loop_pressure<enums::Etimestepping_type::Transient>
(const CFVM::CFVM_cell &FVM_cell,CFlow_cell &Flow_cell,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);



template<enums::Econv_flux_type econv_flux_type >
void CFlow_solver::calc_flux_enthalpy(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    solvers_common::Convective_flux<econv_flux_type>
            ( Matrix_solver,FVM_face, Flow_face.Mdot,Flow_face.left_cell_ptr->Enthalpy, Flow_face.right_cell_ptr->Enthalpy);

    double T_face = solvers_common::face_average(FVM_face,Flow_face.left_cell_ptr->Temperature,Flow_face.right_cell_ptr->Temperature);
    double avgD = Fluid_properties->get_viscosity(T_face);
    avgD = avgD / Fluid_properties->Prandlt_number * FVM_face.diffusion_ksi_A;

    solvers_common::Diffusive_flux( Matrix_solver, FVM_face, avgD);
}

template void CFlow_solver::calc_flux_enthalpy<enums::Econv_flux_type::First_order_upwind >
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::calc_flux_enthalpy<enums::Econv_flux_type::Blended_central >
(const CFVM::CFVM_face &FVM_face, CFlow_face &Flow_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver, const std::shared_ptr<Cfluid_properties> &Fluid_properties);


void CFlow_solver::calc_flux_enthalpy(const CFVM::CFVM_BC_face &FVM_BC_face, CFlow_BC_face &Flow_BC_face,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver,const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    // If empty BC, do nothing
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::empty){

    }

    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::symmetry){

    }

    // If boundary wall
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::wall){
        // Adiabatic do nothing
        if(Flow_BC_face.eThermal_BC_type == enums::Ethermal_BC_type::adiabatic)
        {

        }
        // Constant temp will be written
        if(Flow_BC_face.eThermal_BC_type == enums::Ethermal_BC_type::constant_temperature)
        {

        }
    }

    // If boundary velocity inlet
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::velocity_inlet){

        double value = 0.0;
        //Calc RHS for left cell
        value = Flow_BC_face.Mdot *Flow_BC_face.Enthalpy;

        // Fill RHS for the left cell
        Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,-value);

        double avgD = (Fluid_properties->get_viscosity(Flow_BC_face.Temperature));
        avgD *= Fluid_properties->Prandlt_number * FVM_BC_face.diffusion_ksi_A;

        solvers_common::Diffusive_flux_dirichlet_BC( Matrix_solver, FVM_BC_face, avgD, Flow_BC_face.Enthalpy);

    }

    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::supersonic_velocity_inlet){

        double value = 0.0;
        //Calc RHS for left cell
        value = Flow_BC_face.Mdot *Flow_BC_face.Enthalpy;

        // Fill RHS for the left cell
        Matrix_solver->b_vector.add_insert_element(FVM_BC_face.left_cell_index,-value);

        double avgD = (Fluid_properties->get_viscosity(Flow_BC_face.Temperature));
        avgD *= Fluid_properties->Prandlt_number * FVM_BC_face.diffusion_ksi_A;

        solvers_common::Diffusive_flux_dirichlet_BC( Matrix_solver, FVM_BC_face, avgD, Flow_BC_face.Enthalpy);

    }

    // If boundary is pressure outlet
    if(Flow_BC_face.eflow_BC_type == enums::Eflow_BC_type::pressure_outlet){

        solvers_common::Convective_flux_neumann_zero_BC(Matrix_solver,FVM_BC_face, Flow_BC_face.Mdot);
    }
}

template<enums::Etimestepping_type etimestepping_type>
void CFlow_solver::cell_loop_enthalpy(const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver,const std::shared_ptr<Cfluid_properties> &Fluid_properties)
{
    if constexpr(etimestepping_type == enums::Etimestepping_type::Transient)
    {
        double dum = (FVM_cell.volume / Flow_cell.delta_t);
        double dum2 = dum*Flow_cell.Density;
        Matrix_solver->matrix.add_insert_element(FVM_cell.index,FVM_cell.index,dum2);
        dum2 = dum*(Flow_cell.Enthalpy_last_t*Flow_cell.Density_last_t);
        dum2 = dum2 + (Flow_cell.Pressure - Flow_cell.Pressure_last_t)/Flow_cell.delta_t*FVM_cell.volume;
        Matrix_solver->b_vector.add_insert_element(FVM_cell.index,dum2);
    }

    if constexpr(etimestepping_type == enums::Etimestepping_type::Steady)
    {
        double dum = Matrix_solver->matrix.get_element(FVM_cell.index,FVM_cell.index);
        Matrix_solver->matrix.insert_element(FVM_cell.index,FVM_cell.index,dum/this->implicit_underrelaxation);
        dum = (1.0-this->implicit_underrelaxation)/this->implicit_underrelaxation*dum*Flow_cell.Enthalpy_last_t;
        Matrix_solver->b_vector.add_insert_element(FVM_cell.index,dum);
    }

    Matrix_solver->x_vector.insert_element(FVM_cell.index , Flow_cell.Enthalpy);
}

template void CFlow_solver::cell_loop_enthalpy<enums::Etimestepping_type::Steady>
(const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver,const std::shared_ptr<Cfluid_properties> &Fluid_properties);

template void CFlow_solver::cell_loop_enthalpy<enums::Etimestepping_type::Transient>
(const CFVM::CFVM_cell &FVM_cell, CFlow_cell &Flow_cell,
                        std::shared_ptr<CMatrix_solver> &Matrix_solver,const std::shared_ptr<Cfluid_properties> &Fluid_properties);
