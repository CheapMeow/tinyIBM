#include "output_helper.h"
#include "scope_timer.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    const double cx = 1.85, cy = 4, r = 0.15;

    const double dx = 1.0 / 64.0;
    const double dy = dx;

    const int time_step  = 16000;
    const int output_gap = 50;

    const double xmin = 0, xmax = 8;
    const double ymin = 0, ymax = 8;

    const int Nx = (xmax - xmin) / dx + 1;
    const int Ny = (ymax - ymin) / dx + 1;

    // Create spatial grids
    std::vector<double> x(Nx), y(Ny);

    for (int i = 0; i < Nx; i++)
        x[i] = xmin + i * dx;

    for (int i = 0; i < Ny; i++)
        y[i] = ymin + i * dy;

    const double F                   = 1.00;
    const double dt                  = 0.01;
    const double Re                  = 100.0;
    const double rho                 = 1.00;
    const double dynamic_viscosity   = 2.0 * r / Re;
    const double kinematic_viscosity = dynamic_viscosity / rho;
    const double nu                  = kinematic_viscosity;

    std::vector<std::vector<double>> b(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> p(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> u(Nx, std::vector<double>(Ny));
    std::vector<std::vector<double>> v(Nx, std::vector<double>(Ny));

    std::vector<std::vector<double>> zeta(Nx, std::vector<double>(Ny));

    std::stringstream ss;
    ss << "./ib_cylinder_Re" << Re << "_Nx" << Nx << "_Ny" << Ny << "_H" << dx << "_cx" << cx << "_cy" << cy << "_r"
       << r;
    std::string root_dir = ss.str();

    if (fs::exists(root_dir))
    {
        uintmax_t remove_count = fs::remove_all(root_dir);
        std::cout << "remove_all(" << std::quoted(root_dir) << "): remove " << remove_count << "files" << std::endl;
    }

    std::error_code err;
    bool            create_root_dir_okey = fs::create_directory(root_dir, err);
    std::cout << "create_directory(" << std::quoted(root_dir) << "), result = " << create_root_dir_okey << "\n";
    std::cout << "err.value() = " << err.value() << "\n"
              << " err.message() = " << err.message() << "\n";

    // Lid-driven cavity condition
    for (int i = 0; i < Ny; i++)
    {
        u[Nx - 1][i] = 1.0;
    }

    for (int n = 0; n < time_step; n++)
    {   
        std::stringstream ss_iter;
        ss_iter << "iteration " << n << " / " << time_step; 
        BeginScopeTimer(ss_iter.str());

// Source term b calculation
#pragma omp parallel for schedule(dynamic)
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                b[i][j] =
                    rho *
                    (1.0 / dt * ((u[i + 1][j] - u[i - 1][j]) / (2.0 * dx) + (v[i][j + 1] - v[i][j - 1]) / (2.0 * dy)) -
                     ((u[i + 1][j] - u[i - 1][j]) / (2.0 * dx)) * ((u[i + 1][j] - u[i - 1][j]) / (2.0 * dx)) -
                     2.0 * ((u[i][j + 1] - u[i][j - 1]) / (2.0 * dy) * (v[i + 1][j] - v[i - 1][j]) / (2.0 * dx)) -
                     ((v[i][j + 1] - v[i][j - 1]) / (2.0 * dy)) * ((v[i][j + 1] - v[i][j - 1]) / (2.0 * dy)));
            }
        }

        // Boundary Condition
        for (int j = 1; j < Ny - 1; j++)
            b[0][j] =
                rho * (1.0 / dt * ((u[1][j] - u[Nx - 1][j]) / (2.0 * dx) + (v[0][j + 1] - v[0][j - 1]) / (2.0 * dy)) -
                       ((u[1][j] - u[Nx - 1][j]) / (2.0 * dx)) * ((u[1][j] - u[Nx - 1][j]) / (2.0 * dx)) -
                       2.0 * ((u[0][j + 1] - u[0][j - 1]) / (2.0 * dy) * (v[1][j] - v[Nx - 1][j]) / (2.0 * dx)) -
                       ((v[0][j + 1] - v[0][j - 1]) / (2.0 * dy)) * ((v[0][j + 1] - v[0][j - 1]) / (2.0 * dy)));
        for (int j = 1; j < Ny - 1; j++)
            b[Nx - 1][j] =
                rho *
                (1.0 / dt *
                     ((u[0][j] - u[Nx - 2][j]) / (2.0 * dx) + (v[Nx - 1][j + 1] - v[Nx - 1][j - 1]) / (2.0 * dy)) -
                 ((u[0][j] - u[Nx - 2][j]) / (2.0 * dx)) * ((u[0][j] - u[Nx - 2][j]) / (2.0 * dx)) -
                 2.0 * ((u[Nx - 1][j + 1] - u[Nx - 1][j - 1]) / (2.0 * dy) * (v[0][j] - v[Nx - 2][j]) / (2.0 * dx)) -
                 ((v[Nx - 1][j + 1] - v[Nx - 1][j - 1]) / (2.0 * dy)) *
                     ((v[Nx - 1][j + 1] - v[Nx - 1][j - 1]) / (2.0 * dy)));

        // Pressure Poisson equation
        for (int it = 1; it <= 50; it++)
        {
            std::vector<std::vector<double>> pn = p;

#pragma omp parallel for schedule(dynamic)
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    p[i][j] = ((pn[i + 1][j] + pn[i - 1][j]) * dy * dy + (pn[i][j + 1] + pn[i][j - 1]) * dx * dx -
                               b[i][j] * dx * dx * dy * dy) /
                              (2. * (dx * dx + dy * dy));
                }
            }

            // Boundary Condition
            for (int j = 1; j < Ny - 1; j++)
                p[0][j] = ((pn[1][j] + pn[Nx - 1][j]) * dy * dy + (pn[0][j + 1] + pn[0][j - 1]) * dx * dx -
                           b[0][j] * dx * dx * dy * dy) /
                          (2. * (dx * dx + dy * dy));
            for (int j = 1; j < Ny - 1; j++)
                p[Nx - 1][j] = ((pn[0][j] + pn[Nx - 2][j]) * dy * dy +
                                (pn[Nx - 1][j + 1] + pn[Nx - 1][j - 1]) * dx * dx - b[Nx - 1][j] * dx * dx * dy * dy) /
                               (2. * (dx * dx + dy * dy));

            // Wall boundary conditions
            for (int i = 0; i < Ny; i++)
            {
                p[i][0]      = p[i][1];
                p[i][Nx - 1] = p[i][Nx - 2];
            }
        }

        std::vector<std::vector<double>> un = u;
        std::vector<std::vector<double>> vn = v;

// Velocity update
#pragma omp parallel for schedule(dynamic)
        for (int i = 1; i < Nx - 1; i++)
        {
            for (int j = 1; j < Ny - 1; j++)
            {
                u[i][j] = un[i][j] - un[i][j] * dt / dx * (un[i][j] - un[i - 1][j]) -
                          vn[i][j] * dt / dy * (un[i][j] - un[i][j - 1]) -
                          dt / (2 * rho * dx) * (p[i + 1][j] - p[i - 1][j]) +
                          nu * dt *
                              ((un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / (dx * dx) +
                               (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / (dy * dy)) +
                          F * dt;
                v[i][j] = vn[i][j] - un[i][j] * dt / dx * (vn[i][j] - vn[i - 1][j]) -
                          vn[i][j] * dt / dy * (vn[i][j] - vn[i][j - 1]) -
                          dt / (2 * rho * dy) * (p[i][j + 1] - p[i][j - 1]) +
                          nu * dt *
                              ((vn[i + 1][j] - 2 * vn[i][j] + vn[i - 1][j]) / (dx * dx) +
                               (vn[i][j + 1] - 2 * vn[i][j] + vn[i][j - 1]) / (dy * dy));
            }
        }

        // Boundary Condition
        for (int j = 1; j < Ny - 1; j++)
        {
            u[0][j] = un[0][j] - un[0][j] * dt / dx * (un[0][j] - un[Nx - 1][j]) -
                      vn[0][j] * dt / dy * (un[0][j] - un[0][j - 1]) - dt / (2 * rho * dx) * (p[1][j] - p[Nx - 1][j]) +
                      nu * dt *
                          ((un[1][j] - 2 * un[0][j] + un[Nx - 1][j]) / (dx * dx) +
                           (un[0][j + 1] - 2 * un[0][j] + un[0][j - 1]) / (dy * dy)) +
                      F * dt;
            v[0][j] = vn[0][j] - un[0][j] * dt / dx * (vn[0][j] - vn[Nx - 1][j]) -
                      vn[0][j] * dt / dy * (vn[0][j] - vn[0][j - 1]) -
                      dt / (2 * rho * dy) * (p[0][j + 1] - p[0][j - 1]) +
                      nu * dt *
                          ((vn[1][j] - 2 * vn[0][j] + vn[Nx - 1][j]) / (dx * dx) +
                           (vn[0][j + 1] - 2 * vn[0][j] + vn[0][j - 1]) / (dy * dy));
        }

        for (int j = 1; j < Ny - 1; j++)
        {
            u[Nx - 1][j] = un[Nx - 1][j] - un[Nx - 1][j] * dt / dx * (un[Nx - 1][j] - un[Nx - 2][j]) -
                           vn[Nx - 1][j] * dt / dy * (un[Nx - 1][j] - un[Nx - 1][j - 1]) -
                           dt / (2 * rho * dx) * (p[0][j] - p[Nx - 2][j]) +
                           nu * dt *
                               ((un[0][j] - 2 * un[Nx - 1][j] + un[Nx - 2][j]) / (dx * dx) +
                                (un[Nx - 1][j + 1] - 2 * un[Nx - 1][j] + un[Nx - 1][j - 1]) / (dy * dy)) +
                           F * dt;
            v[Nx - 1][j] = vn[Nx - 1][j] - un[Nx - 1][j] * dt / dx * (vn[Nx - 1][j] - vn[Nx - 2][j]) -
                           vn[Nx - 1][j] * dt / dy * (vn[Nx - 1][j] - vn[Nx - 1][j - 1]) -
                           dt / (2 * rho * dy) * (p[Nx - 1][j + 1] - p[Nx - 1][j - 1]) +
                           nu * dt *
                               ((vn[0][j] - 2 * vn[Nx - 1][j] + vn[Nx - 2][j]) / (dx * dx) +
                                (vn[Nx - 1][j + 1] - 2 * vn[Nx - 1][j] + vn[Nx - 1][j - 1]) / (dy * dy));
        }

        for (int i = 0; i < Nx; i++)
        {
            u[i][0]      = 0;
            u[i][Ny - 1] = 0;
            v[i][0]      = 0;
            v[i][Ny - 1] = 0;
        }

        if (n % output_gap == 0)
        {
#pragma omp parallel for schedule(dynamic)
            for (int i = 1; i < Nx - 1; i++)
            {
                for (int j = 1; j < Ny - 1; j++)
                {
                    zeta[i][j] = (v[i + 1][j] - v[i - 1][j] - u[i][j + 1] + u[i][j - 1]) / 2.0 / dx;
                }
            }
            field_output(u, root_dir, "u", n);
            field_output(v, root_dir, "v", n);
            field_output(zeta, root_dir, "zeta", n);
        }
    }

    return 0;
}