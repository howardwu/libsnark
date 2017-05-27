/** @file
*****************************************************************************

Implementation of interfaces for a ppzkSNARK for R1CS.

See r1cs_gg_ppzksnark.hpp .

*****************************************************************************
* @author     This file is part of libsnark, developed by SCIPR Lab
*             and contributors (see AUTHORS).
* @copyright  MIT license (see LICENSE file)
*****************************************************************************/

#ifndef R1CS_GG_PPZKSNARK_TCC_
#define R1CS_GG_PPZKSNARK_TCC_

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <sstream>

#include "common/profiling.hpp"
#include "common/utils.hpp"
#include "algebra/scalar_multiplication/multiexp.hpp"
#include "algebra/scalar_multiplication/kc_multiexp.hpp"
#include "reductions/r1cs_to_qap/r1cs_to_qap.hpp"

namespace libsnark {

template<typename ppT>
bool r1cs_gg_ppzksnark_proving_key<ppT>::operator==(const r1cs_gg_ppzksnark_proving_key<ppT> &other) const
{
    return (this->alpha_g1 == other.alpha_g1 &&
            this->beta_g1 == other.beta_g1 &&
            this->beta_g2 == other.beta_g2 &&
            this->delta_g1 == other.delta_g1 &&
            this->delta_g2 == other.delta_g2 &&
            this->A_query == other.A_query &&
            this->B_query == other.B_query &&
            this->H_query == other.H_query &&
            this->L_query == other.L_query &&
            this->constraint_system == other.constraint_system);
}

template<typename ppT>
std::ostream& operator<<(std::ostream &out, const r1cs_gg_ppzksnark_proving_key<ppT> &pk)
{
    out << pk.alpha_g1 << OUTPUT_NEWLINE;
    out << pk.beta_g1 << OUTPUT_NEWLINE;
    out << pk.beta_g2 << OUTPUT_NEWLINE;
    out << pk.delta_g1 << OUTPUT_NEWLINE;
    out << pk.delta_g2 << OUTPUT_NEWLINE;
    out << pk.A_query;
    out << pk.B_query;
    out << pk.H_query;
    out << pk.L_query;
    out << pk.constraint_system;

    return out;
}

template<typename ppT>
std::istream& operator>>(std::istream &in, r1cs_gg_ppzksnark_proving_key<ppT> &pk)
{
    in >> pk.alpha_g1;
    consume_OUTPUT_NEWLINE(in);
    in >> pk.beta_g1;
    consume_OUTPUT_NEWLINE(in);
    in >> pk.beta_g2;
    consume_OUTPUT_NEWLINE(in);
    in >> pk.delta_g1;
    consume_OUTPUT_NEWLINE(in);
    in >> pk.delta_g2;
    consume_OUTPUT_NEWLINE(in);
    in >> pk.A_query;
    in >> pk.B_query;
    in >> pk.H_query;
    in >> pk.L_query;
    in >> pk.constraint_system;

    return in;
}

template<typename ppT>
bool r1cs_gg_ppzksnark_verification_key<ppT>::operator==(const r1cs_gg_ppzksnark_verification_key<ppT> &other) const
{
    return (this->alpha_g1_beta_g2 == other.alpha_g1_beta_g2 &&
            this->gamma_g2 == other.gamma_g2 &&
            this->delta_g2 == other.delta_g2 &&
            this->gamma_ABC_g1 == other.gamma_ABC_g1);
}

template<typename ppT>
std::ostream& operator<<(std::ostream &out, const r1cs_gg_ppzksnark_verification_key<ppT> &vk)
{
    out << vk.alpha_g1_beta_g2 << OUTPUT_NEWLINE;
    out << vk.gamma_g2 << OUTPUT_NEWLINE;
    out << vk.delta_g2 << OUTPUT_NEWLINE;
    out << vk.gamma_ABC_g1 << OUTPUT_NEWLINE;

    return out;
}

template<typename ppT>
std::istream& operator>>(std::istream &in, r1cs_gg_ppzksnark_verification_key<ppT> &vk)
{
    in >> vk.alpha_g1_beta_g2;
    consume_OUTPUT_NEWLINE(in);
    in >> vk.gamma_g2;
    consume_OUTPUT_NEWLINE(in);
    in >> vk.delta_g2;
    consume_OUTPUT_NEWLINE(in);
    in >> vk.gamma_ABC_g1;
    consume_OUTPUT_NEWLINE(in);

    return in;
}

template<typename ppT>
bool r1cs_gg_ppzksnark_processed_verification_key<ppT>::operator==(const r1cs_gg_ppzksnark_processed_verification_key<ppT> &other) const
{
    return (this->vk_alpha_g1_beta_g2 == other.vk_alpha_g1_beta_g2 &&
            this->vk_gamma_g2_precomp == other.vk_gamma_g2_precomp &&
            this->vk_delta_g2_precomp == other.vk_delta_g2_precomp &&
            this->gamma_ABC_g1 == other.gamma_ABC_g1);
}

template<typename ppT>
std::ostream& operator<<(std::ostream &out, const r1cs_gg_ppzksnark_processed_verification_key<ppT> &pvk)
{
    out << pvk.vk_alpha_g1_beta_g2 << OUTPUT_NEWLINE;
    out << pvk.vk_gamma_g2_precomp << OUTPUT_NEWLINE;
    out << pvk.vk_delta_g2_precomp << OUTPUT_NEWLINE;
    out << pvk.gamma_ABC_g1 << OUTPUT_NEWLINE;

    return out;
}

template<typename ppT>
std::istream& operator>>(std::istream &in, r1cs_gg_ppzksnark_processed_verification_key<ppT> &pvk)
{
    in >> pvk.vk_alpha_g1_beta_g2;
    consume_OUTPUT_NEWLINE(in);
    in >> pvk.vk_gamma_g2_precomp;
    consume_OUTPUT_NEWLINE(in);
    in >> pvk.vk_delta_g2_precomp;
    consume_OUTPUT_NEWLINE(in);
    in >> pvk.gamma_ABC_g1;
    consume_OUTPUT_NEWLINE(in);

    return in;
}

template<typename ppT>
bool r1cs_gg_ppzksnark_proof<ppT>::operator==(const r1cs_gg_ppzksnark_proof<ppT> &other) const
{
    return (this->g_A == other.g_A &&
            this->g_B == other.g_B &&
            this->g_C == other.g_C);
}

template<typename ppT>
std::ostream& operator<<(std::ostream &out, const r1cs_gg_ppzksnark_proof<ppT> &proof)
{
    out << proof.g_A << OUTPUT_NEWLINE;
    out << proof.g_B << OUTPUT_NEWLINE;
    out << proof.g_C << OUTPUT_NEWLINE;

    return out;
}

template<typename ppT>
std::istream& operator>>(std::istream &in, r1cs_gg_ppzksnark_proof<ppT> &proof)
{
    in >> proof.g_A;
    consume_OUTPUT_NEWLINE(in);
    in >> proof.g_B;
    consume_OUTPUT_NEWLINE(in);
    in >> proof.g_C;
    consume_OUTPUT_NEWLINE(in);

    return in;
}

template<typename ppT>
r1cs_gg_ppzksnark_verification_key<ppT> r1cs_gg_ppzksnark_verification_key<ppT>::dummy_verification_key(const size_t input_size)
{
    r1cs_gg_ppzksnark_verification_key<ppT> result;
    result.alpha_g1_beta_g2 = Fr<ppT>::random_element() * GT<ppT>::random_element();
    result.gamma_g2 = G2<ppT>::random_element();
    result.delta_g2 = G2<ppT>::random_element();

    G1<ppT> base = G1<ppT>::random_element();
    G1_vector<ppT> v;
    for (size_t i = 0; i < input_size; ++i)
    {
        v.emplace_back(G1<ppT>::random_element());
    }

    result.gamma_ABC_g1 = accumulation_vector<G1<ppT> >(std::move(base), std::move(v));

    return result;
}

template <typename ppT>
r1cs_gg_ppzksnark_keypair<ppT> r1cs_gg_ppzksnark_generator(const r1cs_gg_ppzksnark_constraint_system<ppT> &r1cs)
{
    enter_block("Call to r1cs_gg_ppzksnark_generator");

    /* Make the B_query "lighter" if possible */
    r1cs_gg_ppzksnark_constraint_system<ppT> r1cs_copy(r1cs);
    r1cs_copy.swap_AB_if_beneficial();

    /* Generate secret randomness */
    const Fr<ppT> t = Fr<ppT>::random_element();
    const Fr<ppT> alpha = Fr<ppT>::random_element();
    const Fr<ppT> beta = Fr<ppT>::random_element();
    const Fr<ppT> gamma = Fr<ppT>::random_element();
    const Fr<ppT> delta = Fr<ppT>::random_element();
    const Fr<ppT> gamma_inverse = gamma.inverse();
    const Fr<ppT> delta_inverse = delta.inverse();

    /* A quadratic arithmetic program evaluated at t. */
    qap_instance_evaluation<Fr<ppT> > qap = r1cs_to_qap_instance_map_with_evaluation(r1cs_copy, t);

    print_indent(); printf("* QAP number of variables: %zu\n", qap.num_variables());
    print_indent(); printf("* QAP pre degree: %zu\n", r1cs_copy.constraints.size());
    print_indent(); printf("* QAP degree: %zu\n", qap.degree());
    print_indent(); printf("* QAP number of input variables: %zu\n", qap.num_inputs());

    enter_block("Compute query densities");
    size_t non_zero_At = 0;
    size_t non_zero_Bt = 0;
    for (size_t i = 0; i < qap.num_variables() + 1; ++i)
    {
        if (!qap.At[i].is_zero())
        {
            ++non_zero_At;
        }
        if (!qap.Bt[i].is_zero())
        {
            ++non_zero_Bt;
        }
    }
    leave_block("Compute query densities");

    /* qap.{At,Bt,Ct,Ht} are now in unspecified state, but we do not use them later */
    Fr_vector<ppT> At = std::move(qap.At);
    Fr_vector<ppT> Bt = std::move(qap.Bt);
    Fr_vector<ppT> Ct = std::move(qap.Ct);
    Fr_vector<ppT> Ht = std::move(qap.Ht);

    /* The gamma inverse product component: (beta*A_i(t) + alpha*B_i(t) + C_i(t)) * gamma^{-1}. */
    enter_block("Compute gamma_ABC for R1CS verification key");
    Fr_vector<ppT> gamma_ABC;
    gamma_ABC.reserve(qap.num_inputs());
    
    const Fr<ppT> gamma_ABC_0 = (beta * At[0] + alpha * Bt[0] + Ct[0]) * gamma_inverse;
    for (size_t i = 1; i < qap.num_inputs() + 1; ++i)
    {
        gamma_ABC.emplace_back((beta * At[i] + alpha * Bt[i] + Ct[i]) * gamma_inverse);
    }
    leave_block("Compute gamma_ABC for R1CS verification key");

    /* The delta inverse product component: (beta*A_i(t) + alpha*B_i(t) + C_i(t)) * delta^{-1}. */
    enter_block("Compute L query for R1CS proving key");
    Fr_vector<ppT> Lt;
    Lt.reserve(qap.num_variables() - qap.num_inputs());
    
    const size_t Lt_offset = qap.num_inputs() + 1;
    for (size_t i = 0; i < qap.num_variables() - qap.num_inputs(); ++i)
    {
        Lt.emplace_back((beta * At[Lt_offset + i] + alpha * Bt[Lt_offset + i] + Ct[Lt_offset + i]) * delta_inverse);
    }
    leave_block("Compute L query for R1CS proving key");

    /**
     * Note that H for Groth's proof system is degree d-2, but the QAP
     * reduction returns coefficients for degree d polynomial H (in
     * style of PGHR-type proof systems)
     */
    Ht.resize(Ht.size() - 2);

#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
#else
    const size_t chunks = 1;
#endif

    enter_block("Generating G1 MSM window table");
    const G1<ppT> g1_generator = G1<ppT>::random_element();
    const size_t g1_scalar_count = non_zero_At + non_zero_Bt + qap.num_variables();
    const size_t g1_scalar_size = Fr<ppT>::size_in_bits();
    const size_t g1_window_size = get_exp_window_size<G1<ppT> >(g1_scalar_count);

    print_indent(); printf("* G1 window: %zu\n", g1_window_size);
    window_table<G1<ppT> > g1_table = get_window_table(g1_scalar_size, g1_window_size, g1_generator);
    leave_block("Generating G1 MSM window table");

    enter_block("Generating G2 MSM window table");
    const G2<ppT> G2_gen = G2<ppT>::random_element();
    const size_t g2_scalar_count = non_zero_Bt;
    const size_t g2_scalar_size = Fr<ppT>::size_in_bits();
    size_t g2_window_size = get_exp_window_size<G2<ppT> >(g2_scalar_count);

    print_indent(); printf("* G2 window: %zu\n", g2_window_size);
    window_table<G2<ppT> > g2_table = get_window_table(g2_scalar_size, g2_window_size, G2_gen);
    leave_block("Generating G2 MSM window table");

    enter_block("Generate R1CS proving key");
    G1<ppT> alpha_g1 = alpha * g1_generator;
    G1<ppT> beta_g1 = beta * g1_generator;
    G2<ppT> beta_g2 = beta * G2_gen;
    G1<ppT> delta_g1 = delta * g1_generator;
    G2<ppT> delta_g2 = delta * G2_gen;

    enter_block("Generate queries");
    enter_block("Compute the A-query", false);
    G1_vector<ppT> A_query = batch_exp(g1_scalar_size, g1_window_size, g1_table, At);
#ifdef USE_MIXED_ADDITION
    batch_to_special<G1<ppT> >(A_query);
#endif
    leave_block("Compute the A-query", false);

    enter_block("Compute the B-query", false);
    knowledge_commitment_vector<G2<ppT>, G1<ppT> > B_query = kc_batch_exp(Fr<ppT>::size_in_bits(), g2_window_size, g1_window_size, g2_table, g1_table, Fr<ppT>::one(), Fr<ppT>::one(), Bt, chunks);
    leave_block("Compute the B-query", false);

    enter_block("Compute the H-query", false);
    G1_vector<ppT> H_query = batch_exp_with_coeff(g1_scalar_size, g1_window_size, g1_table, qap.Zt * delta_inverse, Ht);
    leave_block("Compute the H-query", false);

    enter_block("Compute the L-query", false);
    G1_vector<ppT> L_query = batch_exp(g1_scalar_size, g1_window_size, g1_table, Lt);
#ifdef USE_MIXED_ADDITION
    batch_to_special<G1<ppT> >(L_query);
#endif
    leave_block("Compute the L-query", false);
    leave_block("Generate queries");

    leave_block("Generate R1CS proving key");

    enter_block("Generate R1CS verification key");
    GT<ppT> alpha_g1_beta_g2 = ppT::reduced_pairing(alpha_g1, beta_g2);
    G2<ppT> gamma_g2 = gamma * G2_gen;

    enter_block("Encode gamma_ABC for R1CS verification key");
    G1<ppT> gamma_ABC_g1_0 = gamma_ABC_0 * g1_generator;
    G1_vector<ppT> gamma_ABC_g1_values = batch_exp(g1_scalar_size, g1_window_size, g1_table, gamma_ABC);
    leave_block("Encode gamma_ABC for R1CS verification key");
    leave_block("Generate R1CS verification key");

    leave_block("Call to r1cs_gg_ppzksnark_generator");

    accumulation_vector<G1<ppT> > gamma_ABC_g1(std::move(gamma_ABC_g1_0), std::move(gamma_ABC_g1_values));

    r1cs_gg_ppzksnark_verification_key<ppT> vk = r1cs_gg_ppzksnark_verification_key<ppT>(alpha_g1_beta_g2,
                                                                                         gamma_g2,
                                                                                         delta_g2,
                                                                                         gamma_ABC_g1);

    r1cs_gg_ppzksnark_proving_key<ppT> pk = r1cs_gg_ppzksnark_proving_key<ppT>(std::move(alpha_g1),
                                                                               std::move(beta_g1),
                                                                               std::move(beta_g2),
                                                                               std::move(delta_g1),
                                                                               std::move(delta_g2),
                                                                               std::move(A_query),
                                                                               std::move(B_query),
                                                                               std::move(H_query),
                                                                               std::move(L_query),
                                                                               std::move(r1cs_copy));

    pk.print_size();
    vk.print_size();

    return r1cs_gg_ppzksnark_keypair<ppT>(std::move(pk), std::move(vk));
}

template <typename ppT>
r1cs_gg_ppzksnark_proof<ppT> r1cs_gg_ppzksnark_prover(const r1cs_gg_ppzksnark_proving_key<ppT> &pk,
                                                      const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input,
                                                      const r1cs_gg_ppzksnark_auxiliary_input<ppT> &auxiliary_input)
{
    enter_block("Call to r1cs_gg_ppzksnark_prover");

#ifdef DEBUG
    assert(pk.constraint_system.is_satisfied(primary_input, auxiliary_input));
#endif

    enter_block("Compute the polynomial H");
    const qap_witness<Fr<ppT> > qap_wit = r1cs_to_qap_witness_map(pk.constraint_system, primary_input, auxiliary_input, Fr<ppT>::zero(), Fr<ppT>::zero(), Fr<ppT>::zero());

    /* We are dividing degree 2(d-1) polynomial by degree d polynomial
       and not adding a PGHR-style ZK-patch, so our H is degree d-2 */
    assert(!qap_wit.coefficients_for_H[qap_wit.degree()-2].is_zero());
    assert(qap_wit.coefficients_for_H[qap_wit.degree()-1].is_zero());
    assert(qap_wit.coefficients_for_H[qap_wit.degree()].is_zero());
    leave_block("Compute the polynomial H");

#ifdef DEBUG
    const Fr<ppT> t = Fr<ppT>::random_element();
    qap_instance_evaluation<Fr<ppT> > qap_inst = r1cs_to_qap_instance_map_with_evaluation(pk.constraint_system, t);
    assert(qap_inst.is_satisfied(qap_wit));
#endif

    /* Choose two random field elements for prover zero-knowledge. */
    const Fr<ppT> r = Fr<ppT>::random_element();
    const Fr<ppT> s = Fr<ppT>::random_element();

#ifdef DEBUG
    assert(qap_wit.coefficients_for_ABCs.size() == qap_wit.num_variables());
    assert(pk.A_query.size() == qap_wit.num_variables()+1);
    assert(pk.B_query.domain_size() == qap_wit.num_variables()+1);
    assert(pk.H_query.size() == qap_wit.degree() - 1);
    assert(pk.L_query.size() == qap_wit.num_variables() - qap_wit.num_inputs());
#endif

#ifdef MULTICORE
    const size_t chunks = omp_get_max_threads(); // to override, set OMP_NUM_THREADS env var or call omp_set_num_threads()
#else
    const size_t chunks = 1;
#endif

    enter_block("Compute the proof");

    enter_block("Compute evaluation to A-query", false);
    // TODO: sort out indexing
    Fr_vector<ppT> const_padded_assignment(1, Fr<ppT>::one());
    const_padded_assignment.insert(const_padded_assignment.end(), qap_wit.coefficients_for_ABCs.begin(), qap_wit.coefficients_for_ABCs.end());

    G1<ppT> evaluation_At = multi_exp_with_mixed_addition<G1<ppT>, Fr<ppT> >(
        pk.A_query.begin(),
        pk.A_query.begin() + qap_wit.num_variables() + 1,
        const_padded_assignment.begin(),
        const_padded_assignment.begin() + qap_wit.num_variables() + 1,
        chunks,
        true);
    leave_block("Compute evaluation to A-query", false);

    enter_block("Compute evaluation to B-query", false);
    knowledge_commitment<G2<ppT>, G1<ppT> > evaluation_Bt = kc_multi_exp_with_mixed_addition<G2<ppT>, G1<ppT>, Fr<ppT> >(
        pk.B_query,
        0,
        qap_wit.num_variables() + 1,
        const_padded_assignment.begin(),
        const_padded_assignment.begin() + qap_wit.num_variables() + 1,
        chunks,
        true);
    leave_block("Compute evaluation to B-query", false);

    enter_block("Compute evaluation to H-query", false);
    G1<ppT> evaluation_Ht = multi_exp<G1<ppT>, Fr<ppT> >(
        pk.H_query.begin(),
        pk.H_query.begin() + (qap_wit.degree() - 1),
        qap_wit.coefficients_for_H.begin(),
        qap_wit.coefficients_for_H.begin() + (qap_wit.degree() - 1),
        chunks,
        true);
    leave_block("Compute evaluation to H-query", false);

    enter_block("Compute evaluation to L-query", false);
    G1<ppT> evaluation_Lt = multi_exp_with_mixed_addition<G1<ppT>, Fr<ppT> >(
        pk.L_query.begin(),
        pk.L_query.end(),
        const_padded_assignment.begin() + qap_wit.num_inputs() + 1,
        const_padded_assignment.begin() + qap_wit.num_variables() + 1,
        chunks,
        true);
    leave_block("Compute evaluation to L-query", false);

    /* A = alpha + sum_i(a_i*A_i(t)) + r*delta */
    G1<ppT> g1_A = pk.alpha_g1 + evaluation_At + r * pk.delta_g1;

    /* B = beta + sum_i(a_i*B_i(t)) + s*delta */
    G1<ppT> g1_B = pk.beta_g1 + evaluation_Bt.h + s * pk.delta_g1;
    G2<ppT> g2_B = pk.beta_g2 + evaluation_Bt.g + s * pk.delta_g2;

    /* C = sum_i(a_i*((beta*A_i(t) + alpha*B_i(t) + C_i(t)) + H(t)*Z(t))/delta) + A*s + r*b - r*s*delta */
    G1<ppT> g1_C = evaluation_Ht + evaluation_Lt + s *  g1_A + r * g1_B - (r * s) * pk.delta_g1;

    leave_block("Compute the proof");

    leave_block("Call to r1cs_gg_ppzksnark_prover");

    r1cs_gg_ppzksnark_proof<ppT> proof = r1cs_gg_ppzksnark_proof<ppT>(std::move(g1_A), std::move(g2_B), std::move(g1_C));
    proof.print_size();

    return proof;
}

template <typename ppT>
r1cs_gg_ppzksnark_processed_verification_key<ppT> r1cs_gg_ppzksnark_verifier_process_vk(const r1cs_gg_ppzksnark_verification_key<ppT> &vk)
{
    enter_block("Call to r1cs_gg_ppzksnark_verifier_process_vk");

    r1cs_gg_ppzksnark_processed_verification_key<ppT> pvk;
    pvk.vk_alpha_g1_beta_g2 = vk.alpha_g1_beta_g2;
    pvk.vk_gamma_g2_precomp = ppT::precompute_G2(vk.gamma_g2);
    pvk.vk_delta_g2_precomp = ppT::precompute_G2(vk.delta_g2);
    pvk.gamma_ABC_g1 = vk.gamma_ABC_g1;

    leave_block("Call to r1cs_gg_ppzksnark_verifier_process_vk");

    return pvk;
}

template <typename ppT>
bool r1cs_gg_ppzksnark_online_verifier_weak_IC(const r1cs_gg_ppzksnark_processed_verification_key<ppT> &pvk,
                                               const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input,
                                               const r1cs_gg_ppzksnark_proof<ppT> &proof)
{
    enter_block("Call to r1cs_gg_ppzksnark_online_verifier_weak_IC");
    assert(pvk.gamma_ABC_g1.domain_size() >= primary_input.size());

    enter_block("Accumulate input");
    const accumulation_vector<G1<ppT> > accumulated_IC = pvk.gamma_ABC_g1.template accumulate_chunk<Fr<ppT> >(primary_input.begin(), primary_input.end(), 0);
    const G1<ppT> &acc = accumulated_IC.first;
    leave_block("Accumulate input");

    bool result = true;

    enter_block("Check if the proof is well-formed");
    if (!proof.is_well_formed())
    {
        if (!inhibit_profiling_info)
        {
            print_indent(); printf("At least one of the proof elements does not lie on the curve.\n");
        }
        result = false;
    }
    leave_block("Check if the proof is well-formed");

    enter_block("Online pairing computations");
    enter_block("Check QAP divisibility");
    const G1_precomp<ppT> proof_g_A_precomp = ppT::precompute_G1(proof.g_A);
    const G2_precomp<ppT> proof_g_B_precomp = ppT::precompute_G2(proof.g_B);
    const G1_precomp<ppT> proof_g_C_precomp = ppT::precompute_G1(proof.g_C);
    const G1_precomp<ppT> acc_precomp = ppT::precompute_G1(acc);

    const Fqk<ppT> QAP1 = ppT::miller_loop(proof_g_A_precomp,  proof_g_B_precomp);
    const Fqk<ppT> QAP2 = ppT::double_miller_loop(
        acc_precomp, pvk.vk_gamma_g2_precomp,
        proof_g_C_precomp, pvk.vk_delta_g2_precomp);
    const GT<ppT> QAP = ppT::final_exponentiation(QAP1 * QAP2.unitary_inverse());

    if (QAP != pvk.vk_alpha_g1_beta_g2)
    {
        if (!inhibit_profiling_info)
        {
            print_indent(); printf("QAP divisibility check failed.\n");
        }
        result = false;
    }
    leave_block("Check QAP divisibility");
    leave_block("Online pairing computations");

    leave_block("Call to r1cs_gg_ppzksnark_online_verifier_weak_IC");

    return result;
}

template<typename ppT>
bool r1cs_gg_ppzksnark_verifier_weak_IC(const r1cs_gg_ppzksnark_verification_key<ppT> &vk,
                                        const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input,
                                        const r1cs_gg_ppzksnark_proof<ppT> &proof)
{
    enter_block("Call to r1cs_gg_ppzksnark_verifier_weak_IC");
    r1cs_gg_ppzksnark_processed_verification_key<ppT> pvk = r1cs_gg_ppzksnark_verifier_process_vk<ppT>(vk);
    bool result = r1cs_gg_ppzksnark_online_verifier_weak_IC<ppT>(pvk, primary_input, proof);
    leave_block("Call to r1cs_gg_ppzksnark_verifier_weak_IC");
    return result;
}

template<typename ppT>
bool r1cs_gg_ppzksnark_online_verifier_strong_IC(const r1cs_gg_ppzksnark_processed_verification_key<ppT> &pvk,
                                                 const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input,
                                                 const r1cs_gg_ppzksnark_proof<ppT> &proof)
{
    bool result = true;
    enter_block("Call to r1cs_gg_ppzksnark_online_verifier_strong_IC");

    if (pvk.gamma_ABC_g1.domain_size() != primary_input.size())
    {
        print_indent(); printf("Input length differs from expected (got %zu, expected %zu).\n", primary_input.size(), pvk.gamma_ABC_g1.domain_size());
        result = false;
    }
    else
    {
        result = r1cs_gg_ppzksnark_online_verifier_weak_IC(pvk, primary_input, proof);
    }

    leave_block("Call to r1cs_gg_ppzksnark_online_verifier_strong_IC");
    return result;
}

template<typename ppT>
bool r1cs_gg_ppzksnark_verifier_strong_IC(const r1cs_gg_ppzksnark_verification_key<ppT> &vk,
                                          const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input,
                                          const r1cs_gg_ppzksnark_proof<ppT> &proof)
{
    enter_block("Call to r1cs_gg_ppzksnark_verifier_strong_IC");
    r1cs_gg_ppzksnark_processed_verification_key<ppT> pvk = r1cs_gg_ppzksnark_verifier_process_vk<ppT>(vk);
    bool result = r1cs_gg_ppzksnark_online_verifier_strong_IC<ppT>(pvk, primary_input, proof);
    leave_block("Call to r1cs_gg_ppzksnark_verifier_strong_IC");
    return result;
}

template<typename ppT>
bool r1cs_gg_ppzksnark_affine_verifier_weak_IC(const r1cs_gg_ppzksnark_verification_key<ppT> &vk,
                                               const r1cs_gg_ppzksnark_primary_input<ppT> &primary_input,
                                               const r1cs_gg_ppzksnark_proof<ppT> &proof)
{
    enter_block("Call to r1cs_gg_ppzksnark_affine_verifier_weak_IC");
    assert(vk.gamma_ABC_g1.domain_size() >= primary_input.size());

    affine_ate_G2_precomp<ppT> pvk_vk_gamma_g2_precomp = ppT::affine_ate_precompute_G2(vk.gamma_g2);
    affine_ate_G2_precomp<ppT> pvk_vk_delta_g2_precomp = ppT::affine_ate_precompute_G2(vk.delta_g2);

    enter_block("Accumulate input");
    const accumulation_vector<G1<ppT> > accumulated_IC = vk.gamma_ABC_g1.template accumulate_chunk<Fr<ppT> >(primary_input.begin(), primary_input.end(), 0);
    const G1<ppT> &acc = accumulated_IC.first;
    leave_block("Accumulate input");

    bool result = true;

    enter_block("Check if the proof is well-formed");
    if (!proof.is_well_formed())
    {
        if (!inhibit_profiling_info)
        {
            print_indent(); printf("At least one of the proof elements does not lie on the curve.\n");
        }
        result = false;
    }
    leave_block("Check if the proof is well-formed");

    enter_block("Check QAP divisibility");
    const affine_ate_G1_precomp<ppT> proof_g_A_precomp = ppT::affine_ate_precompute_G1(proof.g_A);
    const affine_ate_G2_precomp<ppT> proof_g_B_precomp = ppT::affine_ate_precompute_G2(proof.g_B);
    const affine_ate_G1_precomp<ppT> proof_g_C_precomp = ppT::affine_ate_precompute_G1(proof.g_C);
    const affine_ate_G1_precomp<ppT> acc_precomp = ppT::affine_ate_precompute_G1(acc);

    const Fqk<ppT> QAP_miller = ppT::affine_ate_e_times_e_over_e_miller_loop(
        acc_precomp, pvk_vk_gamma_g2_precomp,
        proof_g_C_precomp, pvk_vk_delta_g2_precomp,
        proof_g_A_precomp,  proof_g_B_precomp);
    const GT<ppT> QAP = ppT::final_exponentiation(QAP_miller.unitary_inverse());

    if (QAP != vk.alpha_g1_beta_g2)
    {
        if (!inhibit_profiling_info)
        {
            print_indent(); printf("QAP divisibility check failed.\n");
        }
        result = false;
    }
    leave_block("Check QAP divisibility");

    leave_block("Call to r1cs_gg_ppzksnark_affine_verifier_weak_IC");

    return result;
}

} // libsnark
#endif // R1CS_GG_PPZKSNARK_TCC_
