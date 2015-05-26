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

#include "exp_fem.hpp"
#include <sstream>
#include <gmm/gmm_condition_number.h>

void exponential::derivate(short_type wrt)
{
    GMM_ASSERT1(wrt <= 2, "Derivate wrt invalid variable");
    c[0] = 0;
    c[1] *= complex_type(0, 1) * (wrt == 0 ? k[0] : 0);
    if(c.size() >= 4)
    {
        c[2] *= complex_type(0, 1) * (wrt == 1 ? k[1] : 0);

        c[3] *= complex_type(0, 1) * ((wrt == 0 || wrt == 1) ? k[wrt] : 0);
    }
    if(c.size() == 8)
    {
        c[4] *= complex_type(0, 1) * (wrt == 2 ? k[2] : 0);

        c[5] *= complex_type(0, 1) * ((wrt == 0 || wrt == 2) ? k[wrt] : 0);
        c[6] *= complex_type(0, 1) * ((wrt == 1 || wrt == 2) ? k[wrt] : 0);

        c[7] *= complex_type(0, 1) * k[wrt];
    }
}

bool aequals(double a, double b)
{
    return std::abs(a-b) < 1e-8;
}


exponential EXP_fem_base::get_exponential(size_type i, const domain_cvector &k,
        const base_matrix &G) const
{
    size_type dim = gmm::mat_nrows(G);
    GMM_ASSERT1(nbase == (1u << dim), "Invalid dimension");
    GMM_ASSERT1(i < nbase, "Base out of range");


    gmm::dense_matrix<complex_type> A(nbase, nbase);
    plain_cvector b(nbase);
    b[i] = 1;
    plain_cvector row(nbase);

    for(size_type r = 0; r < nbase; r++)
    {
        exponential::term_values(gmm::mat_col(G, r), k, row);
        gmm::copy(row, gmm::mat_row(A, r));
    }

    exponential func(dim);
    func.k = k;

    gmm::lu_solve(A, func.c, b);

    return func;
}

std::string c_to_str(complex_type c)
{
    std::ostringstream s;
    s<<c;
    return s.str();
}

void EXP_fem_base::prepare_get(const getfem::fem_interpolation_context &c,
                             base_tensor &t, size_type dim) const
{
    GMM_ASSERT1(c.is_convex_num_valid(), "Need convex num");
    GMM_ASSERT1(c.have_G(), "Need G");
    GMM_ASSERT1(gmm::mat_ncols(c.G()) == nbase, "Must have 2 << dim points");

    bgeot::multi_index mi(2);
    mi[0] = nbase; mi[1] = dim;
    t.adjust_sizes(mi);
}

exponential EXP_fem_base::get_base(
        const getfem::fem_interpolation_context &c, size_type i) const
{
    size_type cv = c.convex_num();
    size_type ci = cv * nbase + i;
    if(cache[ci].empty()) cache[ci] = get_exponential(i, k[cv], c.G());
    return cache[ci];
}

void EXP_fem_base::real_base_value(const getfem::fem_interpolation_context &c,
                             base_tensor &t, exp_real_part) const
{
    prepare_get(c, t, 1);
    for(size_type i = 0; i < nbase; i++)
    {
        exponential ex = get_base(c, i);

        for(size_type j = 0; j < nbase; j++)
        {
            complex_type expected = j == i ?
                complex_type(1,0) : complex_type(0,0);

            complex_type value = ex.eval(gmm::mat_col(c.G(), j));

            GMM_ASSERT1(aequals(value.real(), expected.real()) &&
                    aequals(value.imag(), expected.imag()),
                    "Invalid base value (" + c_to_str(value) + " != " +
                    c_to_str(expected) + ")");
        }

        t[i] = std::real(ex.eval(c.xreal()));
    }
}

void EXP_fem_base::real_grad_base_value(const getfem::fem_interpolation_context &c,
                             base_tensor &t, exp_real_part) const
{
    size_type dim = gmm::mat_nrows(c.G());
    prepare_get(c, t, dim);

    for(size_type i = 0; i < nbase; i++)
    {
        for(size_type j = 0; j < dim; j++)
        {
            exponential ex = get_base(c, i);
            ex.derivate(j);
            t(i,j) = std::real(ex.eval(c.xreal()));
        }
    }
}

void EXP_fem_base::real_base_value(const getfem::fem_interpolation_context &c,
                             base_tensor &t, exp_imag_part) const
{
    prepare_get(c, t, 1);
    for(size_type i = 0; i < nbase; i++)
    {
        exponential ex = get_base(c, i);
        t[i] = std::imag(ex.eval(c.xreal()));
    }
}

void EXP_fem_base::real_grad_base_value(const getfem::fem_interpolation_context &c,
                             base_tensor &t, exp_imag_part) const
{
    size_type dim = gmm::mat_nrows(c.G());
    prepare_get(c, t, dim);
    for(size_type i = 0; i < nbase; i++)
    {
        for(size_type j = 0; j < dim; j++)
        {
            exponential ex = get_base(c, i);
            ex.derivate(j);
            t(i,j) = std::imag(ex.eval(c.xreal()));
        }
    }
}

EXP_fem_base::EXP_fem_base(size_type dim, const std::vector<domain_cvector> &k)
    : k(k), nbase(1u << dim), cache(k.size() * nbase, exponential(dim))
{
}

DAL_SIMPLE_KEY(EXP_fem_key, getfem::pfem);

getfem::pfem new_EXP_fem(size_type dim,
        const std::vector<domain_cvector> &k, exp_real_part)
{

    getfem::virtual_fem *p = new EXP_fem_<exp_real_part>(dim, k);
    dal::add_stored_object(new EXP_fem_key(p), p);
    return p;
}

getfem::pfem new_EXP_fem(size_type dim,
        const std::vector<domain_cvector> &k, exp_imag_part)
{

    getfem::virtual_fem *p = new EXP_fem_<exp_imag_part>(dim, k);
    dal::add_stored_object(new EXP_fem_key(p), p);
    return p;
}
