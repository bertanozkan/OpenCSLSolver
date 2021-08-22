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

#include "cmatrix_solver.h"

#include "croot.h"

#include <chrono>
#ifdef COMPILE_WITH_OPENMP
    #include <omp.h>
#endif

CMatrix_solver::CMatrix_solver(size_t n_rows, size_t n_columns, size_t n_nnz)
    : matrix(n_rows, n_columns, n_nnz), M(n_rows, n_columns, n_nnz) , b_vector(n_rows) , x_vector(n_rows)
{
    diag_index.resize(n_rows);
}

// Class constructor for pre allocation
CMatrix_solver::CCSR::CCSR(size_t n_rows, size_t n_columns, size_t n_nnz)
{
    this->n_rows = n_rows;
    this->n_columns = n_columns;
    this->n_nnz = n_nnz;
    value.reserve(n_nnz);
    row_ptr.resize(n_rows+1,0);
    column_index.reserve(n_nnz);
}

double CMatrix_solver::CCSR::get_element(const size_t i , const size_t j)
{
    for(size_t col = row_ptr[i] ; col < row_ptr[i+1] ; col++)
    {
        if(column_index[col] == j)
            return value[col];
    }
    return 0.0;
}

bool CMatrix_solver::CCSR::insert_element(const size_t i, const size_t j,const double value_inp)
{
    // Check if element is allready inserted
    for(size_t col = row_ptr[i] ; col < row_ptr[i+1] ; col++)
    {
        if(column_index[col] == j)
        {
            value[col] = value_inp;
            return true;
        }
    }
    // If element will be inserted first time to the row
    if(row_ptr[i] == row_ptr[i+1])
    {
        value.insert(value.begin()+row_ptr[i],value_inp);
        column_index.insert(column_index.begin()+row_ptr[i],j);
        for(size_t i_row = i+1; i_row < row_ptr.size() ; i_row++)
        {
            row_ptr[i_row] += 1;
        }
        return true;
    }else{
        for(size_t col = row_ptr[i] ; col < row_ptr[i+1] ; col++)
        {
            // If not insert the value
            if(j < column_index[col])
            {
                value.insert(value.begin()+col,value_inp);
                column_index.insert(column_index.begin()+col,j);
                for(size_t i_row = i+1; i_row < row_ptr.size() ; i_row++)
                {
                    row_ptr[i_row] += 1;
                }
                return true;
            }
        }
        size_t col = row_ptr[i+1] - 1;
        value.insert(value.begin()+col+1,value_inp);
        column_index.insert(column_index.begin()+col+1,j);
        for(size_t i_row = i+1; i_row < row_ptr.size() ; i_row++)
        {
            row_ptr[i_row] += 1;
        }
        return true;
    }
    return false;
}

bool CMatrix_solver::CCSR::add_insert_element(const size_t i,const size_t j, const double value_inp)
{
    // Check if element is allready inserted
    for(size_t col = row_ptr[i] ; col < row_ptr[i+1] ; col++)
    {
        if(column_index[col] == j)
        {
            #pragma omp atomic update
            value[col] += value_inp;
            return true;
        }
    }
    // If element will be inserted first time to the row
    if(row_ptr[i] == row_ptr[i+1])
    {
        value.insert(value.begin()+row_ptr[i],value_inp);
        column_index.insert(column_index.begin()+row_ptr[i],j);
        for(size_t i_row = i+1; i_row < row_ptr.size() ; i_row++)
        {
            row_ptr[i_row] += 1;
        }
        return true;
    }else{
        for(size_t col = row_ptr[i] ; col < row_ptr[i+1] ; col++)
        {
            // If not insert the value
            if(j < column_index[col])
            {
                value.insert(value.begin()+col,value_inp);
                column_index.insert(column_index.begin()+col,j);
                for(size_t i_row = i+1; i_row < row_ptr.size() ; i_row++)
                {
                    row_ptr[i_row] += 1;
                }
                return true;
            }
        }
        size_t col = row_ptr[i+1] - 1;
        value.insert(value.begin()+col+1,value_inp);
        column_index.insert(column_index.begin()+col+1,j);
        for(size_t i_row = i+1; i_row < row_ptr.size() ; i_row++)
        {
            row_ptr[i_row] += 1;
        }
        return true;
    }
    return false;
}

void CMatrix_solver::CCSR::clear()
{
    #ifdef COMPILE_WITH_OPENMP
    #pragma omp parallel
    {
        int thread_count = omp_get_num_threads();
        int thread_num   = omp_get_thread_num();
        size_t chunk_size= value.size() / thread_count;
        auto begin = value.begin();
        std::advance(begin, thread_num * chunk_size);
        auto end = begin;
        if(thread_num == (thread_count - 1)) // last thread iterates the remaining sequence
            end = value.end();
        else
            std::advance(end, chunk_size);
        #pragma omp barrier
        std::fill(begin, end, 0.0);
    }
    #else
    std::fill(value.begin(), value.end(), 0.0);
    #endif
}


void CMatrix_solver::CCSR::copy_same_structure(const CCSR &original)
{
    size_t i;
    #pragma omp parallel for
    for(i = 0 ; i < value.size() ; i++)
        value[i] = original.value[i];
}

// Class constructor for pre allocation
CMatrix_solver::Cvector::Cvector(const int n_values)
{
    this->n_values = n_values;
    value.resize(n_values,0.0);
}

double CMatrix_solver::Cvector::get_element(const size_t i)
{
    return value[i];
}

void CMatrix_solver::Cvector::insert_element(const size_t i,const double value_inp)
{
    value[i] = value_inp;
}

void CMatrix_solver::Cvector::add_insert_element(const size_t i, const double value_inp)
{
    #pragma omp atomic update
    value[i] += value_inp;
}

void CMatrix_solver::Cvector::clear()
{
    #ifdef COMPILE_WITH_OPENMP
    #pragma omp parallel
    {
        int thread_count = omp_get_num_threads();
        int thread_num   = omp_get_thread_num();
        size_t chunk_size= value.size() / thread_count;
        auto begin = value.begin();
        std::advance(begin, thread_num * chunk_size);
        auto end = begin;
        if(thread_num == (thread_count - 1)) // last thread iterates the remaining sequence
            end = value.end();
        else
            std::advance(end, chunk_size);
        #pragma omp barrier
        std::fill(begin, end, 0.0);
    }
    #else
    std::fill(value.begin(), value.end(), 0.0);
    #endif
}

// Operation for return = inp1^T * inp2
double CMatrix_solver::norm(const Cvector &inp1, const Cvector &inp2)
{
    double res = 0.0;
    size_t i;
    #pragma omp parallel for reduction (+:res)
    for(i = 0 ; i < inp1.n_values ; i++)
        res += inp1.value[i] * inp2.value[i];
    return res;
}

// Summation of all elements of vector inp
double CMatrix_solver::sum( const Cvector &inp)
{
    double res = 0.0;
    size_t i;
    #pragma omp parallel for reduction (+:res)
    for(i = 0 ; i < inp.n_values ; i++)
        res += inp.value[i];
    return res;
}

// Summation of ABS value of all elements of vector inp
double CMatrix_solver::ABSsum(const Cvector &inp)
{
    double res = 0.0;
    size_t i;
    #pragma omp parallel for reduction (+:res)
    for(i = 0 ; i < inp.n_values ; i++)
        res += abs(inp.value[i]);
    return res;
}

// Operation for out = inp1 * inp2
void CMatrix_solver::product(Cvector &inp1, Cvector &inp2, Cvector &out)
{
    size_t i;
    #pragma omp parallel for
    for(i = 0 ; i < inp1.n_values ; i++)
        out.value[i] = inp1.value[i] * inp2.value[i];
}

// Operation for out = a * inp1
void CMatrix_solver::axpy(const double a, const Cvector &inp1, Cvector &out)
{
    size_t i;
    #pragma omp parallel for
    for(i = 0 ; i < inp1.n_values ; i++)
        out.value[i] = a * inp1.value[i];
}

// Operation for out = a * inp1 + b * inp2
void CMatrix_solver::axpy(const double a,const Cvector &inp1,const double b , const Cvector &inp2, Cvector &out)
{
    if(b == 1.0)
    {
        if(a == 1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_values ; i++)
                out.value[i] = inp1.value[i] + inp2.value[i];
        }else if(a == -1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_values ; i++)
                out.value[i] = inp2.value[i] - inp1.value[i];
        }else
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_values ; i++)
                out.value[i] = a * inp1.value[i] + inp2.value[i];
        }
    }else if(b == -1.0)
    {
        if(a == 1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_values ; i++)
                out.value[i] = inp1.value[i] - inp2.value[i];
        }else if(a == -1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_values ; i++)
                out.value[i] = - inp2.value[i] - inp1.value[i];
        }else
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_values ; i++)
                out.value[i] = a * inp1.value[i] - inp2.value[i];
        }
    }else
    {
        size_t i;
        #pragma omp parallel for
        for(i = 0 ; i < inp1.n_values ; i++)
            out.value[i] = a * inp1.value[i] + b * inp2.value[i];
    }
}

// Operation for out = a * inp1 * inp2
void CMatrix_solver::gemv(const double a, const CCSR &inp1, const Cvector &inp2, Cvector &out)
{
    if(a == 1.0)
    {
        size_t i;
        #pragma omp parallel for
        for(i = 0 ; i < inp1.n_rows ; i++)
        {
            double dum = 0.0;
            for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
            {
                dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
            }
            out.value[i] = dum;
        }
    }else{
        size_t i;
        #pragma omp parallel for
        for(i = 0 ; i < inp1.n_rows ; i++)
        {
            double dum = 0.0;
            for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
            {
                dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
            }
            out.value[i] = a * dum;
        }
    }
}

// Operation for out = a * inp1 * inp2 + b * inp3
void CMatrix_solver::gemv(const double a, const CCSR &inp1, const Cvector &inp2,const double b, const Cvector &inp3, Cvector &out)
{
    if(a == 1.0)
    {
        if(b == 1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_rows ; i++)
            {
                double dum = 0.0;
                for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
                {
                    dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
                }
                out.value[i] = dum + inp3.value[i];
            }
        }else if(b == -1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_rows ; i++)
            {
                double dum = 0.0;
                for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
                {
                    dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
                }
                out.value[i] = dum - inp3.value[i];
            }
        }else {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_rows ; i++)
            {
                double dum = 0.0;
                for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
                {
                    dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
                }
                out.value[i] = dum + b * inp3.value[i];
            }
        }
    }else{
        if(b == 1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_rows ; i++)
            {
                double dum = 0.0;
                for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
                {
                    dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
                }
                out.value[i] = a * dum + inp3.value[i];
            }
        }else if(b == -1.0)
        {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_rows ; i++)
            {
                double dum = 0.0;
                for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
                {
                    dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
                }
                out.value[i] = a * dum - inp3.value[i];
            }
        }else {
            size_t i;
            #pragma omp parallel for
            for(i = 0 ; i < inp1.n_rows ; i++)
            {
                double dum = 0.0;
                for(size_t j = inp1.row_ptr[i] ; j < inp1.row_ptr[i+1] ; j++)
                {
                    dum += inp1.value[j] * inp2.value[inp1.column_index[j]];
                }
                out.value[i] = a * dum + b * inp3.value[i];
            }
        }
    }
}

void CMatrix_solver::extract_diag_indexes()
{
    size_t i;
    for( i = 0 ; i < M.n_rows ; i++)
    {
        for(size_t j = M.row_ptr[i] ; j < M.row_ptr[i+1] ; j++)
        {
            if(M.column_index[j] == i)
            {
                diag_index[i] = j;
                break;
            }
        }
    }
}


void CMatrix_solver::ILU_decompose(const CCSR &A, CCSR &M, const std::vector<int> &diag_index)
{
    /*int N = M.n_rows;
    std::vector<int> iw(N);

    int jw;
    int jrow;
    int j1;
    int j2;
    double tl;

    for (int k = 0; k < N; k++)
    {
        j1 = M.row_ptr[k];
        j2 = M.row_ptr[k+1];
        for (int j = j1; j < j2; j++)
        {
            iw[M.column_index[j]] = j;
        }
        for (int j = j1; j < j2; j++)
        {
            jrow = M.column_index[j];
            if(jrow >= k)
            {
                //diag_index[k] = j;
                if( (jrow != k) || (M.value[j] == 0.0) )
                {
                    std::cout << "Error! ILU decomposition is failed at step" << k << std::endl;
                    return;
                }
                std::fill(iw.begin(),iw.end(),0);
                break;
            }
            tl = M.value[j] /(M.value[diag_index[jrow]]);
            M.value[j] = tl;
            for (int jj = (diag_index[jrow] + 1); jj < int(M.row_ptr[jrow + 1]); jj++)
            {
                jw = iw[M.column_index[jj]];
                if (jw != 0)
                {
                    M.value[jw] -= tl*(M.value[jj] );
                }
            }
        }
        std::fill(iw.begin(),iw.end(),0);
    }
    */
    // Parallel Ilu
    CCSR Mnew = M;
    for(int i_swep = 0 ; i_swep < 3 ; i_swep++)
    {
        size_t i;
        #pragma omp parallel for
        for(i = 0 ; i < A.n_rows ; i++)
        {
            for(size_t j = A.row_ptr[i] ; j < A.row_ptr[i+1] ; j++)
            {
                if(i > A.column_index[j])
                {
                    double dum = 0.0;
                    for(size_t k = A.row_ptr[i] ; k < j  ; k++)
                    {
                        dum += M.value[k] * M.get_element(M.column_index[k],M.column_index[j]);
                    }
                    Mnew.value[j] = (A.value[j] - dum) / M.value[diag_index[A.column_index[j]]];
                }
                // Diagonal is on upper part
                else
                {
                    double dum = 0.0;
                    for(int k = A.row_ptr[i] ; k < diag_index[i] ; k++)
                    {
                        dum += M.value[k] * M.get_element(M.column_index[k],M.column_index[j]);
                    }
                    Mnew.value[j] = (A.value[j] - dum);
                }
            }
        }
        M = Mnew;
    }
}


void CMatrix_solver::ILU_forsub(const CCSR &M, const std::vector<int> &diag_index, const Cvector &b, Cvector &y)
{

    int N = M.n_rows;
    y.value[0] = b.value[0];

    double sum;
    for (int i = 1; i < N; i++)
    {
        sum = 0;
        for (int j_ = M.row_ptr[i]; j_ < diag_index[i]; j_++)
        {

            sum += M.value[j_] * y.value[M.column_index[j_]];
        }
        y.value[i] = b.value[i] - sum;
    }

/*
    // parallel version
    int N = M.n_rows;
    Cvector ynew = y;
    for(int i_swep = 0 ; i_swep < 3 ; i_swep++)
    {
        int i;
        #pragma omp parallel for
        for ( i = 0; i < N; i++)
        {
            double sum = 0.0;
            for (int j_ = M.row_ptr[i]; j_ < diag_index[i]; j_++)
            {
                sum +=  M.value[j_] * y.value[M.column_index[j_]];
            }
            ynew.value[i] = 0.9*(b.value[i] - sum) + 0.1 * y.value[i];
        }

        y = ynew;
    }*/
}

void CMatrix_solver::ILU_backsub(const CCSR &M, const std::vector<int> &uptr, const Cvector &y, Cvector &x)
{

    int N = M.n_rows;
    x.value[N-1] = y.value[N-1] / M.value[uptr[N-1]];

    double sum;
    for (int i = (int)N - 1; 0 <= i; i--)
    {
        sum = 0;
        for (size_t j_ = uptr[i] + 1; j_ < M.row_ptr[i + 1]; j_++)
        {
            sum += M.value[j_] * x.value[M.column_index[j_]];
        }
        x.value[i] = (y.value[i] - sum)/M.value[uptr[i]];
    }

/*
    // parallel version
    int N = M.n_rows;
    Cvector xnew = x;
    for(int i_swep = 0 ; i_swep < 3 ; i_swep++)
    {
        int i;
        #pragma omp parallel for
        for (i = 0; i < N; i++)
        {
            double sum = 0.0;
            double inv_diagonal = 1.0 / M.value[uptr[i]];
            for (size_t j_ = uptr[i] + 1; j_ < M.row_ptr[i + 1]; j_++)
            {
                sum +=  M.value[j_] * x.value[M.column_index[j_]];
            }
            xnew.value[i] = 0.9*inv_diagonal*(y.value[i] - sum) + 0.1 * x.value[i];
        }
        x = xnew;
    }*/
}


double CMatrix_solver::CG_solve()
{
    int n_iter = 4;
    double res_limit = 1.0e-6;

    double res = 0.0;
    size_t n_row = b_vector.n_values;
    double a = 0.0;
    double b = 0.0;

    Cvector r(n_row);
    Cvector r_next(n_row);
    Cvector p(n_row);
    Cvector p_next(n_row);
    Cvector dum(n_row);

    // r = b - Ax
    gemv(-1.0,matrix,x_vector,1.0,b_vector,r);

    // p = r
    p = r;

    for(int i_iter = 0 ; i_iter < n_iter ; i_iter++)
    {
        gemv(1.0,matrix,p,dum);
        a = norm(r,r) / (norm(p,dum) + 1.0e-40);

        dum = x_vector;
        axpy(1.0,dum,a,p,x_vector);

        gemv((-1.0*a),matrix,p,1.0,r,r_next);

        res = norm(r_next,r_next);

        if(res < res_limit)
            return res;

        b = norm(r_next,r_next) / (norm(r,r) + 1.0e-40);

        axpy(1.0,r_next,b,p,p_next);

        r = r_next;
        p = p_next;
    }
    return res;
}

double CMatrix_solver::ILU_BICGSTAB_solve()
{
    int n_iter = 3;
    double res_limit = 1.0e-20;

    double res = 0.0;
    size_t n_row = b_vector.n_values;
    double rho = 1.0 ,alpha = 1.0 , omega = 1.0, rhon = 1.0 , beta , dum1;
    Cvector r(n_row);
    Cvector p(n_row);
    Cvector v(n_row);
    Cvector r0(n_row);
    Cvector s(n_row);
    Cvector dum(n_row);
    Cvector h(n_row);
    Cvector t(n_row);
    Cvector y(n_row);
    Cvector z(n_row);
    Cvector dum2(n_row);

    M.copy_same_structure(matrix);

    ILU_decompose(matrix,M,diag_index);


    // r = b - Ax
    gemv(-1.0,matrix,x_vector,1.0,b_vector,r);

    r0 = r;

    for(int i_iter = 0 ; i_iter < n_iter ; i_iter++)
    {
        rhon = norm(r0,r);
        beta = rhon / (rho + std::numeric_limits<double>::min()) * alpha/(omega + std::numeric_limits<double>::min());

        // dum = p - v*omega
        axpy(1.0,p,-1.0*omega,v,dum);
        // p = r + dum*beta
        axpy(1.0,r,omega,dum,p);

        ILU_forsub(M,diag_index,p,dum);
        ILU_backsub(M,diag_index,dum,y);

        gemv(1.0,matrix,y,v);

        dum1 = norm(r0,v);
        alpha = rhon / (dum1 + std::numeric_limits<double>::min());

        // h = x_vector + alpha * y
        axpy(1.0,x_vector,alpha,y,h);

        res = norm(r,r) / norm(b_vector,b_vector) + std::numeric_limits<double>::min();
        if(res < res_limit)
            return res;

        //s = r - alpha*v
        axpy(-1.0*alpha,v,1.0,r,s);


        ILU_forsub(M,diag_index,s,dum);
        ILU_backsub(M,diag_index,dum,z);

        gemv(1.0,matrix,z,t);

        ILU_forsub(M,diag_index,t,dum);
        ILU_forsub(M,diag_index,s,dum2);


        omega = norm(dum,dum2) / (norm(dum,dum) + 1.0e-40);

        //x_vector = h + omega * z
        axpy(1.0,h,omega,z,x_vector);

        // r = s - omega * t
        axpy(-1.0*omega,t,1.0,s,r);

    }

    return res;

}

void CMatrix_solver::solve()
{
    std::chrono::steady_clock::time_point begin,end;
    begin = std::chrono::steady_clock::now();

    if (root.configuration->ematrix_solver_type ==  enums::Ematrix_solver_type::CG)
    {
        CG_solve();
    }
    else if (root.configuration->ematrix_solver_type ==  enums::Ematrix_solver_type::ILU_BICSTAB)
    {
        ILU_BICGSTAB_solve();
    }

    end = std::chrono::steady_clock::now();
    root.total_elapsed_matrix_solution_time += std::chrono::duration_cast<std::chrono::microseconds>(end-begin).count();
}
