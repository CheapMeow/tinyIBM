#pragma once

#include <omp.h>

enum BoundaryType
{
    Dirichlet,
    Neumann,
    Periodic,
    NoBound
};

class field
{
public:
    unsigned int Nx, Ny;
    BoundaryType boundary_type_left, boundary_type_right, boundary_type_up, boundary_type_down;

    double* value;

    field(unsigned int _Nx,
          unsigned int _Ny,
          BoundaryType _boundary_type_left,
          BoundaryType _boundary_type_right,
          BoundaryType _boundary_type_up,
          BoundaryType _boundary_type_down)
    {
        Nx                  = _Nx;
        Ny                  = _Ny;
        boundary_type_left  = _boundary_type_left;
        boundary_type_right = _boundary_type_right;
        boundary_type_up    = _boundary_type_up;
        boundary_type_down  = _boundary_type_down;
        value               = new double[_Nx * _Ny];
        for (unsigned int i = 0; i < (_Nx * _Ny); i++)
            value[i] = 0.;
    }

    virtual ~field() { delete[] value; }

    double& to(unsigned int i, unsigned int j) { return value[i * Ny + j]; }

    double operator()(int i, int j)
    {
        if (i < 0)
            if (boundary_type_left == Dirichlet)
                return -1.0 * to(-1 * i, j);
            else if (boundary_type_left == Neumann)
                return to(-1 * i, j);
            else if (boundary_type_left == Periodic)
                return to(Nx + i, j);

        if (i >= Nx)
            if (boundary_type_right == Dirichlet)
                return -1.0 * to(2 * Nx - 2 - i, j);
            else if (boundary_type_right == Neumann)
                return to(2 * Nx - 2 - i, j);
            else if (boundary_type_right == Periodic)
                return to(i - Nx, j);

        if (j < 0)
            if (boundary_type_down == Dirichlet)
                return -1.0 * to(i, -1 * j);
            else if (boundary_type_down == Neumann)
                return to(i, -1 * j);
            else if (boundary_type_down == Periodic)
                return to(i, Ny + j);

        if (j >= Ny)
            if (boundary_type_up == Dirichlet)
                return -1.0 * to(i, 2 * Ny - 2 - j);
            else if (boundary_type_up == Neumann)
                return to(i, 2 * Ny - 2 - j);
            else if (boundary_type_up == Periodic)
                return to(i, i - Ny);

        // NoBound will cause error when passing negetive index i, j
        return to(i, j);
    }

    field& operator=(field& f)
    {
        if (this == &f)
        {
            return *this;
        }

        if (this->Nx != f.Nx || this->Ny != f.Ny)
        {
            return *this;
        }

#pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                this->to(i, j) = f.to(i, j);
            }
        }

        return *this;
    }
};