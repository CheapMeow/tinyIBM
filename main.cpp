#include "operator.h"
#include "field.hpp"
#include "output_helper.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    const double cx = 1.85, cy = 4, r = 0.15;

    const double dx = 1.0 / 64.0;
    const double dy = dx;

    const int time_step  = 100000;
    const int output_gap = 500;

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

    field b(Nx, Ny, Neumann, Neumann, Dirichlet, Dirichlet);
    field p(Nx, Ny, Neumann, Neumann, Neumann, Neumann);
    field pn(Nx, Ny, Neumann, Neumann, Neumann, Neumann);
    field u(Nx, Ny, Neumann, Neumann, Dirichlet, Dirichlet);
    field v(Nx, Ny, Neumann, Neumann, Dirichlet, Dirichlet);
    field un(Nx, Ny, Neumann, Neumann, Dirichlet, Dirichlet);
    field vn(Nx, Ny, Neumann, Neumann, Dirichlet, Dirichlet);

    // debug
    field zeta(Nx, Ny, Neumann, Neumann, Neumann, Neumann);

    std::stringstream ss;
    ss << "./ib_cylinder_Re" << Re << "_Nx" << Nx << "_Ny" << Ny << "_H" << dx << "_cx" << cx << "_cy" << cy
       << "_r" << r;
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

    // Left inlet
    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < Ny; j++)
    {
        u.to(0, j) = 1.0;
    }

    // Up and down wall
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < Nx; i++)
    {
        u.to(i, 0)      = 0;
        u.to(i, Ny - 1) = 0;
        v.to(i, 0)      = 0;
        v.to(i, Ny - 1) = 0;
    }

    for (int n = 0; n <= time_step; n++)
    {
        // Source term b calculation
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                b.to(i, j) =
                    rho *
                    (1.0 / dt * ((u(i + 1, j) - u(i - 1, j)) / (2.0 * dx) + (v(i, j + 1) - v(i, j - 1)) / (2.0 * dy)) -
                     ((u(i + 1, j) - u(i - 1, j)) / (2.0 * dx)) * ((u(i + 1, j) - u(i - 1, j)) / (2.0 * dx)) -
                     2.0 * ((u(i, j + 1) - u(i, j - 1)) / (2.0 * dy) * (v(i + 1, j) - v(i - 1, j)) / (2.0 * dx)) -
                     ((v(i, j + 1) - v(i, j - 1)) / (2.0 * dy)) * ((v(i, j + 1) - v(i, j - 1)) / (2.0 * dy)));
            }
        }

        // Pressure Poisson equation

        for (int it = 1; it <= 50; it++)
        {
            pn = p;
            
            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    p.to(i, j) = ((pn(i + 1, j) + pn(i - 1, j)) * dy * dy + (pn(i, j + 1) + pn(i, j - 1)) * dx * dx -
                                  b(i, j) * dx * dx * dy * dy) /
                                 (2. * (dx * dx + dy * dy));
                }
            }
        }

        un = u;
        vn = v;

        // Velocity update
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                if ((x[i] - cx) * (x[i] - cx) + (y[i] - cy) + (y[i] - cy) < r * r)
                {
                    u.to(i, j) = 0.0;
                    v.to(i, j) = 0.0;
                }
            }
        }

        // Left inlet
        #pragma omp parallel for schedule(dynamic)
        for (int j = 0; j < Ny; j++)
        {
            u.to(0, j) = 1.0;
        }

        // Velocity update
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                u.to(i, j) = un(i, j) + dt * conv_diff_u(un, vn, i, j, nu, dx, dy) +
                             dt / (2 * rho * dx) * (p(i + 1, j) - p(i - 1, j));
                v.to(i, j) = vn(i, j) + dt * conv_diff_v(un, vn, i, j, nu, dx, dy) +
                             dt / (2 * rho * dy) * (p(i, j + 1) - p(i, j - 1));
            }
        }

        // Up and down wall
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < Nx; i++)
        {
            u.to(i, 0)      = 0;
            u.to(i, Ny - 1) = 0;
            v.to(i, 0)      = 0;
            v.to(i, Ny - 1) = 0;
        }

        if (n % output_gap == 0)
        {
            #pragma omp parallel for schedule(dynamic)
            for (int i = 0; i < Nx; i++)
            {
                for (int j = 0; j < Ny; j++)
                {
                    zeta.to(i, j) = vorticity(u, v, i, j, dx, dy);
                }
            }
            field_output(u, root_dir, "u", n);
            field_output(v, root_dir, "v", n);
            field_output(zeta, root_dir, "zeta", n);
        }
    }

    return 0;
}