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

#ifndef CMATRIX_SOLVER_H
#define CMATRIX_SOLVER_H

#include "vector"
#include "math.h"

class CMatrix_solver
{
public:
    CMatrix_solver();
    CMatrix_solver(size_t n_rows, size_t n_columns, size_t n_nnz);

    class CCSR
    {
    public:
        size_t n_rows; // Number of rows
        size_t n_columns; // Number of columns
        size_t n_nnz; // Number of non-zeros
        std::vector<double> value; // Value array
        std::vector<size_t> row_ptr; // Row pointer array
        std::vector<size_t> column_index; // Column index array

        CCSR(){}
        CCSR(size_t n_rows, size_t n_columns, size_t n_nnz);

        double get_element(const size_t i , const size_t j);
        bool insert_element(const size_t i , const size_t j , const double value_inp);
        bool add_insert_element(const size_t i , const size_t j, const double value_inp);
        void clear();
        void copy_same_structure(const CCSR &original);
    };

    class Cvector
    {
    public:
        size_t n_values;
        std::vector<double> value;

        Cvector(){}
        Cvector(const int n_values);

        double get_element(const size_t i);
        void insert_element(const size_t i, const double value_inp);
        void add_insert_element(size_t i, const double value_inp);
        void clear();
    };

    double norm(const Cvector &inp1, const Cvector &inp2);
    static double sum(const Cvector &inp);
    static double ABSsum(const Cvector &inp);
    void product(Cvector &inp1, Cvector &inp2, Cvector &out);
    void axpy(double a, const Cvector &inp1, Cvector &out);
    void axpy(const double a, const Cvector &inp1, const double b , const Cvector &inp2, Cvector &out);
    void gemv(const double a, const CCSR &inp1, const Cvector &inp2, Cvector &out);
    void gemv(const double a, const CCSR &inp1, const Cvector &inp2, const double b, const Cvector &inp3, Cvector &out);

    void extract_diag_indexes();
    void ILU_decompose(const CCSR &A, CCSR &M, const std::vector<int> &uptr);
    void ILU_forsub(const CCSR &M, const std::vector<int> &uptr, const Cvector &b, Cvector &y);
    void ILU_backsub(const CCSR &M, const std::vector<int> &uptr, const Cvector &y, Cvector &x);

    double CG_solve();
    double ILU_BICGSTAB_solve();
    void solve();

    // Linear eaquation system memory
    CCSR matrix;
    Cvector b_vector;
    Cvector x_vector;

    // ILU related memory
    std::vector<int> diag_index;
    CCSR M;

};

#endif // CMATRIX_SOLVER_H
