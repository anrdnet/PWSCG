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

solver::solver(std::vector<size_type> &sizes, complex_expression k_expr,
        complex_expression u_expr, complex_expression f_expr,
        complex_expression k1_expr, complex_expression k2_expr,
        complex_expression k3_expr, bool use_exp)
    : mf_u_re(mesh), mf_u_im(mesh), mf_d(mesh), mim(mesh),
    k_expr(k_expr), u_expr(u_expr), f_expr(f_expr)
{
    size_type dim = sizes.size();

    getfem::regular_unit_mesh(mesh, sizes, bgeot::parallelepiped_geotrans(dim,1));
    getfem::outer_faces_of_mesh(mesh, mesh.region(boundary_region));

    std::ostringstream QK_fem;
    QK_fem<<"FEM_QK("<<dim<<",1)";

    mf_d.set_finite_element(getfem::fem_descriptor(QK_fem.str()));

    if(!use_exp)
    {
        mf_u_re.set_finite_element(getfem::fem_descriptor(QK_fem.str()));
        mf_u_im.set_finite_element(getfem::fem_descriptor(QK_fem.str()));
    }
    else
    {
        plain_cvector k1(mf_d.nb_dof());
        plain_cvector k2(mf_d.nb_dof());
        plain_cvector k3(mf_d.nb_dof());

        complex_interpolate_d(k1_expr, k1);
        complex_interpolate_d(k2_expr, k2);
        complex_interpolate_d(k3_expr, k3);

        dal::bit_vector cs = mesh.convex_index();
        std::vector<domain_cvector> k(cs.size());

        for(dal::bv_visitor c(cs); !c.finished(); ++c)
        {
            getfem::mesh_fem::ind_dof_ct dofs = mf_d.ind_basic_dof_of_element(c);
            size_type n = mf_d.fem_of_element(c)->nb_base(c);
            double dn = n;
            for(size_type i = 0; i < n; ++i)
            {
                k[c] += domain_cvector(k1[dofs[i]] / dn,
                        k2[dofs[i]] / dn, k3[dofs[i]] / dn);
            }
        }

        mf_u_re.set_finite_element(new_EXP_fem(dim, k, exp_real_part()));
        mf_u_im.set_finite_element(new_EXP_fem(dim, k, exp_imag_part()));
    }

    std::ostringstream gauss_im;
    gauss_im<<"IM_GAUSS_PARALLELEPIPED("<<dim<<",8)";
    mim.set_integration_method(getfem::int_method_descriptor(gauss_im.str()));

    GMM_ASSERT1(mf_u_re.nb_dof() == mf_u_im.nb_dof(), "Imaginary and real part "
            "must have the same number of dofs");

    for(size_type i = 0; i != get_ndofs(); i++)
    {
        for(size_type k = 0; k < dim; k++)
            GMM_ASSERT1(std::abs(mf_u_re.point_of_basic_dof(i)[k] -
                    mf_u_im.point_of_basic_dof(i)[k]) < 1e-8, "Imaginary and "
                    "real part must have same location of dofs");
    }
}

void solver::mass_matrix(sparse_cmatrix &M, const getfem::mesh_region &rg)
{
    gmm::part_return<sparse_cmatrix*, gmm::linalg_real_part>::return_type Mr(
            gmm::real_part(M));
    gmm::part_return<sparse_cmatrix*, gmm::linalg_imag_part>::return_type Mi(
            gmm::imag_part(M));

    getfem::generic_assembly ga;
    ga.push_mi(mim);
    ga.push_mf(mf_u_re);
    ga.push_mf(mf_u_im);
    ga.push_mat(Mr);
    ga.push_mat(Mi);

    ga.set("M$1(#1,#1)+=comp(Base(#1).Base(#1)) + comp(Base(#2).Base(#2));"
            "M$2(#1,#1)+=comp(Base(#1).Base(#2)) - comp(Base(#2).Base(#1))");
    ga.assembly(rg);
}

void solver::dirichlet_constraints(sparse_cmatrix &H, plain_cvector &R,
        const plain_cvector &R0)
{
    mass_matrix(H, boundary_region);

    gmm::part_return<plain_cvector*, gmm::linalg_real_part>::return_type Rr =
        gmm::real_part(R);
    gmm::part_return<plain_cvector*, gmm::linalg_imag_part>::return_type Ri =
        gmm::imag_part(R);

    gmm::part_return<const plain_cvector*, gmm::linalg_real_part>::return_type R0r =
        gmm::real_part(R0);
    gmm::part_return<const plain_cvector*, gmm::linalg_imag_part>::return_type R0i =
        gmm::imag_part(R0);

    getfem::generic_assembly ga;
    ga.push_mi(mim);
    ga.push_mf(mf_u_re);
    ga.push_mf(mf_u_im);
    ga.push_mf(mf_d);
    ga.push_data(R0r);
    ga.push_data(R0i);
    ga.push_vec(Rr);
    ga.push_vec(Ri);

    ga.set("R0r=data$1(#1); R0i=data$2(#1);"
            "V$1(#1)+="
            "  comp(Base(#1).Base(#1))(:,i).R0r(i)"
            "+ comp(Base(#2).Base(#2))(:,i).R0r(i)"
            "- comp(Base(#1).Base(#2))(:,i).R0i(i)"
            "+ comp(Base(#2).Base(#1))(:,i).R0i(i);"
            "V$2(#1)+="
            "  comp(Base(#1).Base(#1))(:,i).R0i(i)"
            "+ comp(Base(#2).Base(#2))(:,i).R0i(i)"
            "+ comp(Base(#1).Base(#2))(:,i).R0r(i)"
            "- comp(Base(#2).Base(#1))(:,i).R0r(i)");

    ga.assembly(boundary_region);
}

size_type solver::nullspace(const sparse_cmatrix &H, sparse_cmatrix &N,
        const plain_cvector &R, plain_cvector &C)
{
    std::vector<size_type> nulls;
    std::vector<size_type> range;

    scalar_type norminfH = gmm::mat_maxnorm(H);
    for(size_type i = 0; i < gmm::mat_ncols(H); i++)
    {
        (gmm::vect_norm2(gmm::mat_col(H, i)) < 1e-8 * norminfH ? nulls : range).push_back(i);
    }

    gmm::sub_index range_index(range);

    for(size_type i = 0; i < nulls.size(); i++)
    {
        N(nulls[i],i) = 1;
    }

    sparse_cmatrix Hr(range.size(), range.size());
    plain_cvector Rr(range.size());
    plain_cvector Cr(range.size());
    gmm::copy(gmm::sub_matrix(H, range_index, range_index), Hr);
    gmm::copy(gmm::sub_vector(R, range_index), Rr);

    la_solve(Hr, Cr, Rr);

    gmm::copy(Cr, gmm::sub_vector(C, range_index));

    return nulls.size();
}

sparse_cmatrix solver::apply_dirichlet(
        const sparse_cmatrix &A, const plain_cvector &B, const plain_cvector &C0,
        sparse_cmatrix &rA, plain_cvector &rB, plain_cvector &rC)
{
    size_type ndofs = mf_u_re.nb_dof();

    sparse_cmatrix H(ndofs, ndofs);
    sparse_cmatrix N(ndofs, ndofs);
    plain_cvector R(ndofs);

    dirichlet_constraints(H, R, C0);

    sparse_cmatrix Mvec(ndofs, 1);

    size_type nbcols = nullspace(H, N, R, rC);
    gmm::resize(N, ndofs, nbcols);

    gmm::sub_interval int_N(0, nbcols);

    std::cerr<<"ncols/ndofs: "<<nbcols<<"/"<<ndofs<<std::endl;

    sparse_cmatrix NtA(nbcols, ndofs);
    gmm::mult(gmm::transposed(N), A, NtA);
    gmm::mult(NtA, N, gmm::sub_matrix(rA, int_N, int_N));

    // Nt(B - A U)
    gmm::mult(gmm::transposed(N), B, gmm::sub_vector(rB, int_N));
    gmm::mult_add(NtA, gmm::scaled(rC,-1), gmm::sub_vector(rB, int_N));

    return N;
}

bool solver::la_solve(const sparse_cmatrix &A, plain_cvector &C,
        const plain_cvector &B, const plain_cvector C0)
{
    getfem::linear_solver_gmres_preconditioned_ilutp<sparse_cmatrix, plain_cvector> ls_solver;
    gmm::iteration iter(1.0e-6, 0, 100000);

    if(C0.empty())
    {
        ls_solver(A, C, B, iter);
    }
    else
    {
        sparse_cmatrix As(get_ndofs(), get_ndofs());
        plain_cvector Bs(get_ndofs());

        sparse_cmatrix N = apply_dirichlet(A, B, C0, As, Bs, C);

        size_type nbcols = gmm::mat_ncols(N);
        gmm::resize(As, nbcols, nbcols);
        gmm::resize(Bs, nbcols);

        if(nbcols == 0)
            return true;

        plain_cvector Cs(nbcols);

        ls_solver(As, Cs, Bs, iter);

        gmm::mult_add(N, Cs, C);
    }

    if(iter.converged())
        std::cerr<<"Solution converged in "<<
            iter.get_iteration()<<" iterations."<<std::endl;
    else
        return false;

    return true;
}

void solver::complex_interpolate_d(complex_expression expr, plain_cvector &C)
{
    plain_vector Cpart(get_ndofs());
    getfem::ga_interpolation_Lagrange_fem(model, expr.r, mf_d, Cpart);
    gmm::copy(Cpart, gmm::real_part(C));
    gmm::clear(Cpart);
    getfem::ga_interpolation_Lagrange_fem(model, expr.i, mf_d, Cpart);
    gmm::copy(Cpart, gmm::imag_part(C));
}

bool solver::complex_interpolate_u(complex_expression expr, plain_cvector &C)
{
    sparse_cmatrix M(get_ndofs(), get_ndofs());
    mass_matrix(M);

    plain_cvector B(get_ndofs());
    gmm::part_return<plain_cvector*, gmm::linalg_real_part>::return_type Br(
            gmm::real_part(B));
    gmm::part_return<plain_cvector*, gmm::linalg_imag_part>::return_type Bi(
            gmm::imag_part(B));

    std::ostringstream r_asm;
    r_asm<<"("<<expr.r<<") * Test_u_re + ("<<expr.i<<") * Test_u_im";
    std::ostringstream i_asm;
    i_asm<<"("<<expr.i<<") * Test_u_re + ("<<expr.r<<") * -Test_u_im";

    gmm::sub_interval si(0, get_ndofs());
    plain_vector Btmp(get_ndofs());
    getfem::ga_workspace ws;
    ws.add_fem_variable("u_re", mf_u_re, si, Btmp);
    ws.add_fem_variable("u_im", mf_u_im, si, Btmp);

    ws.set_assembled_vector(Btmp);
    ws.add_expression(r_asm.str(), mim);
    ws.assembly(1);
    gmm::copy(Btmp, Br);

    ws.clear_expressions();
    gmm::clear(Btmp);

    ws.set_assembled_vector(Btmp);
    ws.add_expression(i_asm.str(), mim);
    ws.assembly(1);
    gmm::copy(Btmp, Bi);

    if(!la_solve(M, C, B))
        return false;

    return true;
}

bool solver::solve(plain_cvector &C)
{
    plain_cvector F(get_ndofs());
    complex_interpolate_d(f_expr, F);

    plain_cvector K(get_ndofs());
    complex_interpolate_d(k_expr, K);

    plain_cvector C0(get_ndofs());
    complex_interpolate_u(u_expr, C0);

    sparse_cmatrix A(get_ndofs(), get_ndofs());
    plain_cvector B(get_ndofs());

    gmm::part_return<plain_cvector*, gmm::linalg_real_part>::return_type Fr = gmm::real_part(F);
    gmm::part_return<plain_cvector*, gmm::linalg_imag_part>::return_type Fi = gmm::imag_part(F);
    gmm::part_return<plain_cvector*, gmm::linalg_real_part>::return_type Kr = gmm::real_part(K);
    gmm::part_return<plain_cvector*, gmm::linalg_imag_part>::return_type Ki = gmm::imag_part(K);

    gmm::part_return<plain_cvector*, gmm::linalg_real_part>::return_type Br = gmm::real_part(B);
    gmm::part_return<plain_cvector*, gmm::linalg_imag_part>::return_type Bi = gmm::imag_part(B);
    gmm::part_return<sparse_cmatrix*, gmm::linalg_real_part>::return_type Ar = gmm::real_part(A);
    gmm::part_return<sparse_cmatrix*, gmm::linalg_imag_part>::return_type Ai = gmm::imag_part(A);

    getfem::generic_assembly ga;
    ga.push_mi(mim);
    ga.push_mf(mf_u_re);
    ga.push_mf(mf_u_im);
    ga.push_mf(mf_d);
    ga.push_data(Fr);
    ga.push_data(Fi);
    ga.push_data(Kr);
    ga.push_data(Ki);
    ga.push_mat(Ar);
    ga.push_mat(Ai);
    ga.push_vec(Br);
    ga.push_vec(Bi);

    ga.set("Fre=data$1(#3); Fim=data$2(#3); Kre=data$3(#3); Kim=data$4(#3);"
            "M$1(#1,#1)+="
            "- comp(Grad(#1).Grad(#1))(:, k, :, k)"
            "- comp(Grad(#2).Grad(#2))(:, k, :, k)"
            "+ 2*Kre(i).Kim(i).comp(Base(#3).Base(#2).Base(#1))(i, :, :)"
            "- 2*Kre(i).Kim(i).comp(Base(#3).Base(#1).Base(#2))(i, :, :)"
            "+ Kre(i).Kre(i).comp(Base(#3).Base(#1).Base(#1))(i, :, :)"
            "- Kim(i).Kim(i).comp(Base(#3).Base(#1).Base(#1))(i, :, :)"
            "+ Kre(i).Kre(i).comp(Base(#3).Base(#2).Base(#2))(i, :, :)"
            "- Kim(i).Kim(i).comp(Base(#3).Base(#2).Base(#2))(i, :, :);"
            "V$1(#1)+="
            "  comp(Base(#1).Base(#3))(:, i).Fre(i)"
            "+ comp(Base(#2).Base(#3))(:, i).Fim(i);"
            "M$2(#1,#1)+="
            "- comp(Grad(#1).Grad(#2))(:, k, :, k)"
            "+ comp(Grad(#2).Grad(#1))(:, k, :, k)"
            "+ 2*Kre(i).Kim(i).comp(Base(#3).Base(#1).Base(#1))(i, :, :)"
            "+ 2*Kre(i).Kim(i).comp(Base(#3).Base(#2).Base(#2))(i, :, :)"
            "- Kre(i).Kre(i).comp(Base(#3).Base(#2).Base(#1))(i, :, :)"
            "+ Kim(i).Kim(i).comp(Base(#3).Base(#2).Base(#1))(i, :, :)"
            "+ Kre(i).Kre(i).comp(Base(#3).Base(#1).Base(#2))(i, :, :)"
            "- Kim(i).Kim(i).comp(Base(#3).Base(#1).Base(#2))(i, :, :);"
            "V$2(#1)+="
            "  comp(Base(#1).Base(#3))(:,i).Fim(i)"
            "- comp(Base(#2).Base(#3))(:,i).Fre(i)");

    ga.assembly();

    return la_solve(A, C, B, C0);
}

scalar_type solver::L2_dist(const plain_cvector &C, complex_expression expr)
{
    plain_vector Cr(get_ndofs());
    plain_vector Ci(get_ndofs());

    gmm::copy(gmm::real_part(C), Cr);
    gmm::copy(gmm::imag_part(C), Ci);

    std::ostringstream l2asm;
    l2asm<<"(c_re_u_re - c_im_u_im)*(c_re_u_re - c_im_u_im)"
        " + (c_re_u_im + c_im_u_re)*(c_re_u_im + c_im_u_re)"
        " + 2 * -(c_re_u_re - c_im_u_im) * ("<<expr.r<<")"
        " + 2 * -(c_re_u_im + c_im_u_re) * ("<<expr.i<<")"
        " + ("<<expr.r<<") * ("<<expr.r<<") + ("<<expr.i<<") * ("<<expr.i<<")";

    getfem::ga_workspace ws;
    ws.add_fem_constant("c_re_u_re", mf_u_re, Cr);
    ws.add_fem_constant("c_re_u_im", mf_u_im, Cr);
    ws.add_fem_constant("c_im_u_re", mf_u_re, Ci);
    ws.add_fem_constant("c_im_u_im", mf_u_im, Ci);

    ws.add_expression(l2asm.str(), mim);
    ws.assembly(0);

    scalar_type l2s = ws.assembled_potential();

    return std::sqrt(l2s);
}
