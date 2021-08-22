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

#ifndef ENUMS_H
#define ENUMS_H

namespace enums {


// Geometry enums
enum class Eface_type{tri,quad};
enum class Ecell_type{tetra,wedge,hexahedron};
enum class Emesh_source{from_file,generate_1D};

// Configuration enums
enum class Eproblem_type{Heat_conduction,Fluid_flow};
enum class Einitial_condition_type{Uniform,Shock_tube};
enum class Etimestepping_type{Steady,Transient};
enum class Econv_flux_type{First_order_upwind,Blended_central};
enum class Eflow_type_visc{Viscous,Inviscid};


// Boundary enums
enum class Ethermal_BC_type{constant_temperature,adiabatic};
enum class Eflow_BC_type{empty,symmetry,wall,velocity_inlet,supersonic_velocity_inlet,pressure_outlet};
enum class Ewall_BC_type{no_slip, slip};

// Matrix solver enums
enum class Ematrix_solver_type{CG, ILU_BICSTAB};


}

#endif // ENUMS_H
