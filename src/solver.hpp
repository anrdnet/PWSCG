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
#ifndef SOLVE_H_
#define SOLVE_H_

#include <vector>
#include <getfem/getfem_mesh.h>
#include <getfem/getfem_regular_meshes.h>
#include <getfem/getfem_mesh_fem.h>
#include <getfem/getfem_mesh_im.h>
#include <getfem/getfem_mesh_region.h>
#include <getfem/getfem_assembling.h>
#include <getfem/getfem_generic_assembly.h>
#include <getfem/getfem_interpolation.h>
#include <getfem/getfem_mesh_slicers.h>
#include <getfem/getfem_mesh_slice.h>
#include <getfem/getfem_export.h>
#include <getfem/getfem_integration.h>
#include <gmm/gmm_iter_solvers.h>
#include <gmm/gmm_inoutput.h>
#include <getfem/getfem_inter_element.h>
#include <getfem/getfem_import.h>
#include <getfem/getfem_model_solvers.h>

using getfem::scalar_type;
using getfem::size_type;
using getfem::short_type;
using getfem::short_type;
using getfem::base_node;

typedef getfem::model_real_sparse_vector sparse_vector;
typedef getfem::model_real_sparse_matrix sparse_matrix;
typedef getfem::model_real_plain_vector  plain_vector;
typedef getfem::model_complex_plain_vector  plain_cvector;
typedef getfem::model_complex_sparse_vector sparse_cvector;
typedef getfem::model_complex_sparse_matrix sparse_cmatrix;

struct complex_expression
{
    std::string r;
    std::string i;

    complex_expression() : r("0"), i("0") { }
    complex_expression(std::string r, std::string i) : r(r), i(i) { }
};

class solver
{
    static const size_type boundary_region = 1;

    getfem::model model;
    getfem::mesh mesh;
    getfem::mesh_fem mf_u_re, mf_u_im, mf_d;
    getfem::mesh_im mim;
    complex_expression k_expr,u_expr,f_expr;
    complex_expression k1_expr,k2_expr,k3_expr;

    /** Make mass matrix (mf_u, mf_u).
     */
    void mass_matrix(sparse_cmatrix &M, const getfem::mesh_region &rg
            = getfem::mesh_region::all_convexes());

    void mass_matrix_d(sparse_cmatrix &M, const getfem::mesh_region &rg
            = getfem::mesh_region::all_convexes());

    /** Assemble dirichlet constraints Hc = R.
     * \param R0 The boundary data in mf_u.
     */
    void dirichlet_constraints(sparse_cmatrix &H, plain_cvector &R,
            const plain_cvector &R0);

    /** Given condition HC=R returns a matrix N whos columns span the nullspace
     * and C which satisfies the condition.
     */
    size_type nullspace(const sparse_cmatrix &H, sparse_cmatrix &N,
        const plain_cvector &R, plain_cvector &C);

    /** Apply dirichlet condition.
     * \param A The original A.
     * \param B The original B.
     * \param C0 The values at the boundary in mf_u.
     * \param rA Is set to the reduced A.
     * \param rB Is set to the reduced B.
     * \param rU Is set to the particular U which should be added to the homogeneous solution.
     */
    sparse_cmatrix apply_dirichlet(
            const sparse_cmatrix &A, const plain_cvector &B, const plain_cvector &C0,
            sparse_cmatrix &rA, plain_cvector &rB, plain_cvector &rU);

    /** Sovle the complex system AC=B.
     * \param C0 If set, the values of the solution on the boundary in mf_u.
     */
    bool la_solve(const sparse_cmatrix &A, plain_cvector &C,
            const plain_cvector &B, const plain_cvector C0 = plain_cvector());

    public:
    solver(std::vector<size_type> &sizes, complex_expression k_expr,
        complex_expression u_expr, complex_expression f_expr,
        complex_expression k1_expr, complex_expression k2_expr,
        complex_expression k3_expr, bool use_exp);

    /** Interpolate data in expressions to mf_d.
     */ 
    void complex_interpolate_d(complex_expression expr, plain_cvector &C);
    /** Interpolate data in expressions to mf_u.
     */
    bool complex_interpolate_u(complex_expression expr, plain_cvector &C);

    /** Measure distance in L2-norm between a function in mf_u and expression.
     */
    scalar_type L2_dist(const plain_cvector &C, complex_expression expr);

    getfem::model &get_model() { return model; }
    getfem::mesh &get_mesh() { return mesh; }
    getfem::mesh_fem &get_mf_u_re() { return mf_u_re; }
    getfem::mesh_fem &get_mf_u_im() { return mf_u_im; }
    getfem::mesh_fem &get_mf_d() { return mf_d; }
    getfem::mesh_im &get_mim() { return mim; }
    size_type get_ndofs() { return mf_u_re.nb_dof(); }

    bool solve(plain_cvector &C);
};

#endif
