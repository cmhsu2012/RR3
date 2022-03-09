print(f'\nComputing Kolyvagin data:');

#Construct relative norm as linear transformation from SUKz_over_Qz_ell0 to SUQz_ell0 (note that these matrices live in GF(5))
relative_norm_vectors=[list((SUQz_ell0.log(u.relative_norm()))) for u in SUKz_over_Qz_ell0.gens_values()];
relative_norm_matrix= Matrix(GF(p_value),relative_norm_vectors).transpose();

#Construction delta-acton as linear transformation from SUKz_abs_ell0 to itself
delta_action_vectors=[list(SUKz_over_Qz_ell0.log(Kz_abs_to_Kz_over_Qz(delta(Kz_over_Qz_to_Kz_abs(u))))) for u in SUKz_over_Qz_ell0.gens_values()];
delta_action_matrix= Matrix(GF(p_value),delta_action_vectors).transpose();

#Impose condition that beta is omega^2-isotypic
delta_action_matrix_omega = delta_action_matrix - omega(2)*matrix.identity(ZZ(p_value*(p_value-1)/2+p_value*(p_value-1)));

#Compute gamma in Kz
matrix_for_finding_gamma=Matrix(GF(p_value),relative_norm_matrix.rows()+delta_action_matrix_omega.rows());
c1_vector=vector(GF(p_value),list(SUQz_ell0.log(SUQz_ell0(c1_element))));
c1_vector_with_omega_action = vector(GF(p_value),list(SUQz_ell0.log(SUQz_ell0(c1_element)))+list(0*delta_action_matrix.rows()[0]));
gamma_vector=matrix_for_finding_gamma.solve_right(c1_vector_with_omega_action.column()).column(0);
gamma_Kz_over_Qz=SUKz_over_Qz_ell0.exp(gamma_vector);
gamma_abs = mod_p_powers(SUKz_abs_ell0,Kz_over_Qz_to_Kz_abs(mod_p_powers(SUKz_over_Qz_ell0,gamma_Kz_over_Qz)));

#Check that gamma has omega^2-action
assert omega_action_test(gamma_abs,omega(2)), "gamma omega action";

#Define D0_gamma, D1_gamma, D2_gamma
D0_gamma = mod_p_powers(SUKz_abs_ell0,D_0_operator(sigma,gamma_abs));
D1_gamma = mod_p_powers(SUKz_abs_ell0,D_1_operator(sigma,gamma_abs));
D2_gamma = mod_p_powers(SUKz_abs_ell0,D_2_operator(sigma,gamma_abs));

#Check that relationship between D_0 and D_1
assert D0_gamma==mod_p_powers(SUKz_abs_ell0,Kz_over_Qz_to_Kz_abs(gamma_Kz_over_Qz.relative_norm()));
assert mod_p_powers(SUKz_abs_ell0,SUKz_abs_ell0.exp(vector(SUKz_abs_ell0.log(sigma(D1_gamma)/D1_gamma))+vector(SUKz_abs_ell0.log(D0_gamma))))==1;

#Check that D1_gamma has omega-action
assert omega_action_test(D1_gamma,omega(1)), "D1_gamma omega action";

#Check if D1_gamma is p-power free
if p_power_check(D1_gamma):
    print("\tD1_gamma is computed and p-power free");
    D1_gamma_free = D1_gamma;
else:
    D1_gamma_free = p_power_adjustment(D1_gamma);
    if p_power_check(D1_gamma_free):
        print("\tD1_gamma is computed and p-power free");
    else:
        print("\tFurther adjustment is needed");

#Check that D1_gamma_free has omega-action
assert omega_action_test(D1_gamma_free,omega(1)), "D1_gamma_free omega action";

#Check if D2_gamma is p-power free
if p_power_check(D2_gamma):
    print("\tD2_gamma is computed and p-power free");
    D2_gamma_free = D2_gamma;
else:
    D2_gamma_free = p_power_adjustment(D2_gamma);
    if p_power_check(D2_gamma_free):
       print("\tD2_gamma is computed and p-power free");
    else:
        print("\tFurther adjustment is needed");

#Check that the extension Kz(sqrt[p]{D1_gamma})/Q has correct splitting patterns for primes of characteristic ell0 (there should be 16 ramified primes and 4 primes that are either inert or split)
ramified_primes_D1_gamma, split_primes_D1_gamma,inert_primes_D1_gamma=Ell0_splitting_patterns(D1_gamma_free);
print(f'\tSplitting pattern without adjustments is {ramified_primes_D1_gamma,split_primes_D1_gamma,inert_primes_D1_gamma}');
assert len(ramified_primes_D1_gamma)== (p_value-1)^2;
assert len(split_primes_D1_gamma)==p_value-1 or len(inert_primes_D1_gamma)==p_value-1;