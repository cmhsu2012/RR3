print(f'Computing b2 data:');

#Compute xi in Kz
adjustment_a1_vector=vector(list(SUQz_ell0.log(SUQz_ell0(adjustment_a1_Qz))));

#Construction delta-Kolylvagin-acton as linear transformation from SUKz_abs_ell0 to itself
delta_koly_action_vectors=[list(SUKz_over_Qz_ell0.log(Kz_abs_to_Kz_over_Qz(delta(D_1_operator(sigma,Kz_over_Qz_to_Kz_abs(u)))))) for u in SUKz_over_Qz_ell0.gens_values()];
delta_koly_action_matrix= Matrix(GF(p_value),delta_koly_action_vectors).transpose();

koly_action_vectors=[list(SUKz_over_Qz_ell0.log(Kz_abs_to_Kz_over_Qz(D_1_operator(sigma,Kz_over_Qz_to_Kz_abs(u))))) for u in SUKz_over_Qz_ell0.gens_values()];
koly_action_matrix= Matrix(GF(p_value),koly_action_vectors).transpose();

#Impose condition that xi satisfies delta equality
delta_action_matrix_xi = delta_koly_action_matrix - koly_action_matrix;

#Compute xi in Kz
matrix_for_finding_xi=Matrix(GF(p_value),relative_norm_matrix.rows()+delta_action_matrix_xi.rows());

#correction_vector=vector(GF(p_value),list(SUKz_over_Qz_ell0.log(Kz_abs_to_Kz_over_Qz(delta(b2_candidate_free)^(-2)*b2_candidate_free^2))));
b2_base = ((D2_gamma_free)^(-2)*(D1_gamma_free*adjustment_a1_Kz)^(-1)*delta(D2_gamma_free^(2)*D1_gamma_free*adjustment_a1_Kz))^(GF(p_value)(-2)^(-1));
c1_vector_with_rel_norm = vector(GF(p_value),list(SUQz_ell0.log(SUQz_ell0(adjustment_a1_Qz)))+list(SUKz_over_Qz_ell0.log(Kz_abs_to_Kz_over_Qz(b2_base))));
xi_vector=matrix_for_finding_xi.solve_right(c1_vector_with_rel_norm.column()).column(0);
xi_Kz_over_Qz=SUKz_over_Qz_ell0.exp(xi_vector);
xi_abs = mod_p_powers(SUKz_abs_ell0,Kz_over_Qz_to_Kz_abs(mod_p_powers(SUKz_over_Qz_ell0,xi_Kz_over_Qz)));

print("\txi_abs computed");

#Test that relative norm of xi is correct
assert SUQz_ell0.log(mod_p_powers(SUQz_ell0,mod_p_powers(SUKz_over_Qz_ell0,xi_Kz_over_Qz).relative_norm()))==SUQz_ell0.log(SUQz_ell0(adjustment_a1_Qz));

#Compute Kolyvagin derivative of xi
Dxi=mod_p_powers(SUKz_abs_ell0,D_1_operator(sigma,xi_abs));

#Define candidate for b2
b2_candidate = mod_p_powers(SUKz_abs_ell0,(D2_gamma_free*Dxi)^(-2)*(D1_gamma_free*adjustment_a1_Kz)^(-1));

#Check if b2_candidate is p-power free
if p_power_check(b2_candidate):
    print("\tb2_candidate is computed and p-power free");
    b2_candidate_free = b2_candidate;
else:
    b2_candidate_free = p_power_adjustment(b2_candidate);
    if p_power_check(b2_candidate_free):
        print("\tb2_candidate is computed and p-power free");
    else:
        print("\tFurther adjustment is needed");

#Check if b2_candidate is omega(0)-isotypic
assert omega_action_test(b2_candidate_free,omega(0)), "b2_candidate free does not have the correct omega action"