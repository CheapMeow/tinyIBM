#pragma once

#include "field.hpp"

double ddx(field f, int i, int j, double dx);

double ddy(field f, int i, int j, double dy);

double laplacian(field f, int i, int j, double dx, double dy);

double div(field u, field v, int i, int j, double dx, double dy);

double conv_diff_u(field u, field v, int i, int j, double nu, double dx, double dy);

double conv_diff_v(field u, field v, int i, int j, double nu, double dx, double dy);

double vorticity(field u, field v, int i, int j, double dx, double dy);