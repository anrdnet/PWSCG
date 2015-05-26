/* Copyright 2015 Anders Matheson
* This file is part of PWSCG.
* 
* PWSCG is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* Foobar is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "solver.hpp"
#include "util.hpp"
#include "exp_fem.hpp"
#include <cmath>
#include <cstdlib>

void slice_c(getfem::stored_mesh_slice &sl, getfem::mesh_fem &mf_re,
        getfem::mesh_fem &mf_im, const plain_cvector &C,
        plain_cvector &Uslice)
{
    plain_vector Utmp(sl.nb_points());
    //sl.interpolate does not accept part_vector input
    plain_vector Cr(gmm::vect_size(C));
    plain_vector Ci(gmm::vect_size(C));
    gmm::copy(gmm::real_part(C), Cr);
    gmm::copy(gmm::imag_part(C), Ci);

    // c_re * phi_re
    sl.interpolate(mf_re, Cr, Utmp);

    gmm::add(Utmp, gmm::real_part(Uslice));
    gmm::clear(Utmp);

    // -c_im * phi_im
    sl.interpolate(mf_im, Ci, Utmp);

    gmm::add(gmm::scaled(Utmp, -1), gmm::real_part(Uslice));
    gmm::clear(Utmp);

    // i c_im * phi_re
    sl.interpolate(mf_re, Ci, Utmp);

    gmm::add(Utmp, gmm::imag_part(Uslice));
    gmm::clear(Utmp);

    // i c_re * phi_im
    sl.interpolate(mf_im, Cr, Utmp);

    gmm::add(Utmp, gmm::imag_part(Uslice));
}

void write_complex(getfem::vtk_export &exp, std::string name,
        const plain_cvector &U)
{
    // write_sliced_point_data does not accept part_vector as input
    plain_vector Utmp(gmm::vect_size(U));

    gmm::copy(gmm::real_part(U), Utmp);
    exp.write_sliced_point_data(Utmp, name + "_real");

    gmm::copy(gmm::imag_part(U), Utmp);
    exp.write_sliced_point_data(Utmp, name + "_imag");

}

bool parse_complex(std::string arg, std::string value, std::string prefix,
        complex_expression &expr)
{
    if(arg == "-"+prefix+"r")
    {
        if(expr.r != "0")
            std::cerr<<"Warning: multiple values for "<<
                prefix<<"r given"<<std::endl;
        expr.r = value;
        return true;
    }
    else if(arg == "-"+prefix+"i")
    {
        if(expr.i != "0")
            std::cerr<<"Warning: multiple values for "<<
                prefix<<"i given"<<std::endl;
        expr.i = value;
        return true;
    }
    return false;
}

void print_usage()
{
    std::cerr<<"pwcg {-p} {-n <subdiv>} {-s {<n>}+}+ -{k{,1,2,3},f,u}{r,i} <expr>"<<std::endl;
}

int main(int argc, char *argv[])
{
    int subdiv = 0;
    bool use_exp = true;
    bool skip_norms = false;
    std::vector< std::vector<size_type> > sizes;

    complex_expression k_expr;
    complex_expression u_expr;
    complex_expression f_expr;

    complex_expression k1_expr;
    complex_expression k2_expr;
    complex_expression k3_expr;

    for(int i = 1; i < argc; )
    {
        std::string arg(argv[i]);
        if(arg != "-p" && arg != "-z" && i + 1 >= argc)
        {
            print_usage();
            return 1;
        }
        if(arg == "-p")
        {
            use_exp = false;
            i += 1;
        }
        else if(arg == "-z")
        {
            skip_norms = true;
            i += 1;
        }
        else if(arg == "-n")
        {
            if(subdiv != 0)
                std::cerr<<"Warning: multiple values for subdivide"<<std::endl;
            subdiv = std::strtol(argv[i+1], NULL, 10);
            i += 2;
        }
        else if(arg == "-s")
        {
            char *end = NULL;
            std::vector<size_type> size;

            do
            {
                i++;
                size_type st = std::strtol(argv[i], &end, 10);
                if(end == argv[i] || *end != 0)
                    break;
                size.push_back(st);
            }
            while(i + 1 < argc || (i++, false));

            if(size.size() < 1 || size.size() > 3)
            {
                print_usage();
                return 1;
            }
            sizes.push_back(size);
        }
        else if(parse_complex(arg, argv[i+1], "k", k_expr) ||
                parse_complex(arg, argv[i+1], "u", u_expr) ||
                parse_complex(arg, argv[i+1], "f", f_expr) ||
                parse_complex(arg, argv[i+1], "k1", k1_expr) ||
                parse_complex(arg, argv[i+1], "k2", k2_expr) ||
                parse_complex(arg, argv[i+1], "k3", k3_expr))
        {
            i += 2;
        }
        else
        {
            print_usage();
            return 1;
        }
    }

    double lastL2 = 0.0;

    size_type prev_dim = 0;
    if(!skip_norms)
    {
        std::cout<<"            L2     (h^p)"<<std::endl;
        std::cout<<std::fixed;
        std::cout<<std::setfill(' ');
    }
    for(size_t i = 0; i != sizes.size(); i++)
    {
        solver s(sizes[i], k_expr, u_expr, f_expr, k1_expr, k2_expr, k3_expr, use_exp);

        if(prev_dim != 0 && prev_dim != sizes[i].size())
        {
            std::cout<<"-----------------------"<<std::endl;
            lastL2 = 0.0;
        }
        prev_dim = sizes[i].size();

        plain_cvector U(s.get_ndofs());
        if(!s.solve(U))
        {
            std::cerr<<"Solution did not converge."<<std::endl;
            return 1;
        }

        if(!skip_norms)
        {
            scalar_type L2err = s.L2_dist(U, u_expr);

            std::cout<<std::setprecision(8)<<std::setw(14)<<L2err<<"  ";
            if(lastL2 > 1e-8 && L2err > 1e-8)
                std::cout<<'('<<std::setprecision(2)<<std::setw(8)<<log2(lastL2/L2err)<<')';
            else
                std::cout<<"          ";
            lastL2 = L2err;

            std::cout<<std::endl;
        }

        if(subdiv != 0)
        {
            plain_cvector Uexact_d(s.get_ndofs());
            plain_cvector Uexact(s.get_ndofs());

            s.complex_interpolate_d(u_expr, Uexact_d);
            s.complex_interpolate_u(u_expr, Uexact);

            getfem::stored_mesh_slice sl;
            getfem::slicer_build_stored_mesh_slice storer(sl);
            getfem::mesh_slicer slicer(s.get_mesh());
            slicer.push_back_action(storer);
            slicer.exec(subdiv);

            std::ostringstream fname;
            fname<<"out_"<<i<<".vtk";
            getfem::vtk_export exp(fname.str());
            exp.exporting(sl);

            plain_cvector Uslice(sl.nb_points());
            slice_c(sl, s.get_mf_u_re(), s.get_mf_u_im(), U, Uslice);

            write_complex(exp, "a_u", Uslice);

            gmm::clear(Uslice);
            slice_c(sl, s.get_mf_u_re(), s.get_mf_u_im(), Uexact, Uslice);

            write_complex(exp, "e_u", Uslice);

            exp.write_point_data(s.get_mf_d(),
                    gmm::real_part(Uexact_d), "e_d_real");
            exp.write_point_data(s.get_mf_d(),
                    gmm::imag_part(Uexact_d), "e_d_imag");
        }
    }
    return 0;
}
