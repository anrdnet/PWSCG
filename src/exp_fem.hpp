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
#ifndef EXPONENT_H_
#define EXPONENT_H_

#include <getfem/getfem_config.h>
#include <getfem/getfem_fem.h>
#include <getfem/getfem_models.h>

using getfem::complex_type;
using getfem::scalar_type;
using getfem::base_vector;
using getfem::short_type;
using getfem::size_type;
using getfem::base_matrix;
using getfem::base_tensor;
using getfem::base_node;

typedef getfem::model_complex_plain_vector  plain_cvector;

struct exp_real_part {};
struct exp_imag_part {};

struct domain_cvector
{
    complex_type x[4];

    domain_cvector() { x[0] = 0; x[1] = 0; x[2] = 0; x[3] = 0; }
    domain_cvector(complex_type x_, complex_type y,
            complex_type z = 0, complex_type w = 0)
    { x[0] = x_; x[1] = y; x[2] = z; x[3] = w; }

    complex_type *begin() { return &x[0]; }
    const complex_type *begin() const { return &x[0]; }
    complex_type *end() { return &x[4]; }
    const complex_type *end() const { return &x[4]; }

    complex_type operator[] (size_type i) const { return x[i]; }

    void operator+= (const domain_cvector &o)
    { x[0] += o.x[0]; x[1] += o.x[1]; x[2] += o.x[2]; x[3] += o.x[3]; }
};

class exponential
{
    public:
    plain_cvector c;
    domain_cvector k;

    template<typename ITER>
    complex_type eval(ITER it) const;

    template<typename VEC>
    static void term_values(VEC p, const domain_cvector &k, plain_cvector &res);

    exponential(size_type dim) : c(1u << dim) { }

    void derivate(short_type wrt);
    bool empty() const
    { return gmm::vect_norm1(c) < 1e-8; }
};

template<typename VEC>
void exponential::term_values(VEC p, const domain_cvector &k, plain_cvector &res)
{
    res[0] = 1;
    res[1] = std::exp(complex_type(-k[0].imag()*p[0], k[0].real()*p[0]));

    if(res.size() >= 4)
    {
        res[2] = std::exp(complex_type(-k[1].imag()*p[1], k[1].real()*p[1]));
        res[3] = res[1] * res[2];
    }
    if(res.size() == 8)
    {
        res[4] = std::exp(complex_type(-k[2].imag()*p[2], k[2].real()*p[2]));
        res[5] = res[1] * res[4];
        res[6] = res[2] * res[4];
        res[7] = res[1] * res[2] * res[4];
    }
}

template<typename ITER>
complex_type exponential::eval(ITER it) const
{
    plain_cvector terms(c.size());
    term_values(it, k, terms);
    return gmm::vect_sp(c, terms);
}

class EXP_fem_base
{
    std::vector<domain_cvector> k;
    size_type nbase;
    mutable std::vector<exponential> cache;

    exponential get_exponential(size_type bnum, const domain_cvector &k,
            const base_matrix &G) const;

    void prepare_get(const getfem::fem_interpolation_context &c,
            base_tensor &t, size_type dim) const;
    exponential get_base(
            const getfem::fem_interpolation_context &c, size_type i) const;
    public:
    EXP_fem_base(size_type dim, const std::vector<domain_cvector> &k);

    void real_base_value(const getfem::fem_interpolation_context &c,
                                 base_tensor &t, exp_real_part) const;
    void real_grad_base_value(const getfem::fem_interpolation_context &c,
                                 base_tensor &t, exp_real_part) const;

    void real_base_value(const getfem::fem_interpolation_context &c,
                                 base_tensor &t, exp_imag_part) const;
    void real_grad_base_value(const getfem::fem_interpolation_context &c,
                                 base_tensor &t, exp_imag_part) const;
};

template<typename PART>
class EXP_fem_ : public getfem::virtual_fem
{
    EXP_fem_base base;
    public:
    EXP_fem_(size_type dim, const std::vector<domain_cvector> &k);

    virtual void base_value(const base_node &x, base_tensor &t) const
    { GMM_ASSERT1(false, "Use ref base value"); }
    virtual void grad_base_value(const base_node &x, base_tensor &t) const
    { GMM_ASSERT1(false, "Use ref base value"); }
    virtual void hess_base_value(const base_node &x, base_tensor &t) const
    { GMM_ASSERT1(false, "Use ref base value"); }

    virtual void real_base_value(const getfem::fem_interpolation_context &c,
                                 base_tensor &t, bool withM = true) const;
    virtual void real_grad_base_value(const getfem::fem_interpolation_context &c,
                                 base_tensor &t, bool withM = true) const;
    virtual void real_hess_base_value(const getfem::fem_interpolation_context &c,
                                 base_tensor &t, bool withM = true) const;
};

template<typename PART>
EXP_fem_<PART>::EXP_fem_(size_type dim, const std::vector<domain_cvector> &k)
    : base(dim, k)
{
    this->cvr = bgeot::parallelepiped_of_reference(dim);
    this->dim_ = this->cvr->structure()->dim();
    this->is_equiv = this->is_lag = true;
    this->is_pol = false;
    this->es_degree = 1;
    this->real_element_defined = true;

    this->init_cvs_node();
    size_type R = this->cvr->nb_points();
    for (size_type i = 0; i < R; ++i)
      this->add_node(getfem::lagrange_dof(dim), this->cvr->points()[i]);
}

template<typename PART>
void EXP_fem_<PART>::real_base_value(const getfem::fem_interpolation_context &c,
        base_tensor &t, bool withM) const
{
    base.real_base_value(c, t, PART());
}

template<typename PART>
void EXP_fem_<PART>::real_grad_base_value(const getfem::fem_interpolation_context &c,
        base_tensor &t, bool withM) const
{
    base.real_grad_base_value(c, t, PART());
}

template<typename PART>
void EXP_fem_<PART>::real_hess_base_value(const getfem::fem_interpolation_context &c,
        base_tensor &t, bool withM) const
{
    GMM_ASSERT1(false, "Hess of EXP_fem not implemented");
}

getfem::pfem new_EXP_fem(size_type dim,
        const std::vector<domain_cvector> &k, exp_real_part);

getfem::pfem new_EXP_fem(size_type dim,
        const std::vector<domain_cvector> &k, exp_imag_part);

#endif
