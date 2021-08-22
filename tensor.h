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

#ifndef TENSOR_H
#define TENSOR_H

#include "vector"
#include "math.h"
#include <limits>
#include <iostream>

struct Svector
{
    double val[3];
    Svector()
    {
        val[0] = 0.0 ; val[1] = 0.0 ; val[2] = 0.0 ;
    }

    Svector(double x, double y, double z)
    {
        val[0] = x ; val[1] = y; val[2] = z ;
    }

    double operator [](int i)
    {
        return val[i];
    }

    void clear()
    {
        for(int i = 0 ; i < 3 ; i++)
            this->val[i] = 0.0;
    }



    Svector& operator +=(const Svector& rhs)
    {
        for(int i = 0 ; i < 3 ; i++)
            this->val[i] += rhs.val[i];
        return *this;
    }

    Svector& operator -=(const Svector& rhs)
    {
        for(int i = 0 ; i < 3 ; i++)
            this->val[i] -= rhs.val[i];
        return *this;
    }

    Svector& operator *=(const double rhs)
    {
        for(int i = 0 ; i < 3 ; i++)
            this->val[i] *= rhs;
        return *this;
    }

};

struct Stensor
{
    double val[3][3];
    Stensor()
    {
        for(int i = 0 ; i < 3 ; i++)
            for(int j = 0 ; j < 3 ; j++)
                val[i][j] = 0.0;
    }
    double operator ()(int i, int j)
    {
        return val[i][j];
    }
};

inline Svector operator *(double a,Svector B)
{
    Svector result;
    for(int i = 0 ; i < 3 ; i++)
        result.val[i] = a * B[i];
    return result;
}

inline Svector operator +(Svector A,Svector B)
{
    Svector result;
    for(int i = 0 ; i < 3 ; i++)
        result.val[i] = A[i] + B[i];
    return result;
}

inline Svector operator -(Svector A,Svector B)
{
    Svector result;
    for(int i = 0 ; i < 3 ; i++)
        result.val[i] = A[i] - B[i];
    return result;
}

inline Svector operator *(Svector A,Svector B)
{
    Svector result;
    result.val[0] = A[1] * B[2] - A[2] * B[1];
    result.val[1] = A[2] * B[0] - A[0] * B[2];
    result.val[2] = A[0] * B[1] - A[1] * B[0];
    return result;
}



inline double dot(Svector A,Svector B)
{
    double result = 0.0;
    for(int i = 0 ; i < 3 ; i++)
        result += A[i] * B[i];
    return result;
}

inline Svector dot_div(Svector A,Svector B)
{
    Svector result;
    for(int i = 0 ; i < 3 ; i++)
    {
        if(B[i] == 0.0)
            B.val[i] = std::numeric_limits<double>::min();
        result.val[i] = A.val[i] / B.val[i];
    }
    return result;
}


inline double mag(Svector A)
{
    double result = 0.0;
    for(int i = 0 ; i < 3 ; i++)
        result += A[i] * A[i];
    return sqrt(result);
}


inline Stensor operator *(double a,Stensor B)
{
    Stensor result;
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j < 3 ; j++)
            result.val[i][j] = a * B.val[i][j];
    return result;
}

inline Stensor operator +(Stensor A,Stensor B)
{
    Stensor result;
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j < 3 ; j++)
            result.val[i][j] = A.val[i][j] + B.val[i][j];
    return result;
}

inline Stensor operator -(Stensor A,Stensor B)
{
    Stensor result;
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j < 3 ; j++)
            result.val[i][j] = A.val[i][j] - B.val[i][j];
    return result;
}

inline Svector operator *(Stensor A,Svector B)
{
    Svector result;
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j < 3 ; j++)
            result.val[i] += A.val[i][j] * B.val[j];
    return result;
}

inline Stensor gradient( Svector left, Svector right, Svector distance_unit_vec, double distance)
{
    Stensor res;
    Svector dum = right - left;
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j < 3 ; j++)
            res.val[i][j] = dum.val[i] * distance_unit_vec.val[j] / distance;

    return res;
}

inline Stensor divergence_I( Svector left, Svector right, Svector distance_unit_vec, double distance)
{
    Stensor res;
    double dumd = 0.0;
    Svector dum = right - left;
    for(int i = 0 ; i < 3 ; i++)
        dumd += dum.val[i]* distance_unit_vec.val[i] / distance;

    for(int i = 0 ; i < 3 ; i++)
        res.val[i][i] = dumd;

    return res;
}

inline Stensor transpose (Stensor A)
{
    Stensor res;
    for(int i = 0 ; i < 3 ; i++)
        for(int j = 0 ; j < 3 ; j++)
            res.val[i][j] = A.val[j][i];

    return res;
}




#endif // TENSOR_H
