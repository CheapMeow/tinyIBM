#pragma once

#include "field.hpp"

inline double ddx(field& f, int i, int j, double dx) { return (f(i + 1, j) - f(i, j)) / dx; }

inline double ddy(field& f, int i, int j, double dy) { return (f(i, j + 1) - f(i, j)) / dy; }

inline double laplacian(field& f, int i, int j, double dx, double dy)
{
    return (f(i + 1, j) - 2.0 * f(i, j) + f(i - 1, j)) / dx / dx +
           (f(i, j + 1) - 2.0 * f(i, j) + f(i, j - 1)) / dy / dy;
}

inline double div(field& u, field& v, int i, int j, double dx, double dy) { return ddx(u, i, j, dx) + ddy(v, i, j, dy); }

inline double conv_diff_u(field& u, field& v, int i, int j, double nu, double dx, double dy)
{
    return -u(i, j) * ddx(u, i, j, dx) - u(i, j) * ddy(v, i, j, dy) + nu * laplacian(u, i, j, dx, dy);
}

inline double conv_diff_v(field& u, field& v, int i, int j, double nu, double dx, double dy)
{
    return -u(i, j) * ddx(v, i, j, dx) - v(i, j) * ddy(v, i, j, dy) + nu * laplacian(v, i, j, dx, dy);
}

double vorticity(field& u, field& v, int i, int j, double dx, double dy) { return ddx(v, i, j, dx) - ddy(u, i, j, dy); }