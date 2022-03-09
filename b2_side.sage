if p_ramification:
	#Setup in K/Q for the finite flat adjustment

	#Prime ideal of characteristic p in K/Q
	P_K=K_field.primes_above(p_value)[0];

	#Move b2 to Kz/K
	b2_candidate_base = b2_candidate_free;
	b2_K = Kz_abs_to_Kz_over_K(b2_candidate_base);

	#Impose condition that b2 is omega^0-isotypic (exactly, not up to pth powers)

	#Construction delta-acton as linear transformation from SUKz_abs_ell0 to itself
	delta_action_vectors_K=[list(SUKz_over_K_ell0.log(Kz_abs_to_Kz_over_K(delta(Kz_over_K_to_Kz_abs(u))))) for u in SUKz_over_K_ell0.gens_values()];
	delta_action_matrix_K= Matrix(ZZ,delta_action_vectors_K).transpose();
	delta_action_matrix_trivial = delta_action_matrix_K - omega(0)*matrix.identity(ZZ(p_value*(p_value-1)/2+p_value*(p_value-1)));
	test_matrix_trivial = omega(0)*matrix.identity(GF(p_value),ZZ(p_value*(p_value-1)/2+p_value*(p_value-1)));

	#Correction if needed to make vector integral (rather than fractional) 
	new_b2_K = b2_K;
	correction_term = Kz_abs_to_Kz_over_K(delta(Kz_over_K_to_Kz_abs(new_b2_K)))/(new_b2_K);
	correction_inv = correction_term^(-1);
	matrix_for_finding_correction=Matrix(ZZ,delta_action_matrix_trivial.rows());
	correction_vector=vector(ZZ,list(SUKz_over_K_ell0.log(SUKz_over_K_ell0(correction_inv))));
	square_vector=matrix_for_finding_correction.solve_right(correction_vector.column()).column(0);
	denom_in_vec = lcm(a.denominator() for a in square_vector);
	assert denom_in_vec%p_value!=0, "denominator is not coprime to p";
	power_correction = 1;
	if denom_in_vec >1:
		denom_correction = -1*(p_value^(-1))%denom_in_vec;
		assert (denom_correction*p_value +1)%denom_in_vec ==0, "error computing denominator correction mod denom";
		assert (denom_correction*p_value +1)%p_value ==1, "error computing denominator correction mod p";
		power_correction=denom_correction*p_value +1;
	b2_delta_correction = SUKz_over_K_ell0.exp(power_correction*square_vector);
	assert mod_p_powers(SUKz_over_K_ell0,b2_delta_correction)==1, "delta correction is not a p power";
	print("\tb2_candidate has correct Galois action");

	#Compute finite flat adjustment
	finite_flat_adjustment = [];
	p_free = Kz_abs_to_Kz_over_K(from_K_to_Kz_abs(p_power_adjustment(K_field(p_value))));
	P_K_new = Kz_abs_to_Kz_over_K(from_K_to_Kz_abs(P_K));
	new_b2_delta = b2_delta_correction*new_b2_K^power_correction;
	assert new_b2_delta in K_field, "Galois adjustment did not work";

	for p_index in range(p_value):
		if (P_K_new^2).small_residue(new_b2_delta*p_free^p_index) == 1 or (P_K_new^2).small_residue(new_b2_delta*p_free^p_index) == -1:
			finite_flat_adjustment= finite_flat_adjustment+ [p_index];
	assert len(finite_flat_adjustment)==1, "finite flat adjustment not working";
	finite_flat_exponent = finite_flat_adjustment[0];
	print(f'\tThe finite-flat adjustment exponent is {finite_flat_adjustment[0]}');

else:
	if not (p_value).divides(b2_candidate_free.norm()):
		finite_flat_exponent = 0;
		print(f'\tThe finite-flat adjustment exponent is {finite_flat_exponent}');
	else:
		print("\tAdjustment needed for condition at p. Check this");

b2_finite_flat = b2_candidate_free*p_value^finite_flat_exponent;

#Mazur-Tate derivative (in F_p)
zeta_MT = GF(p_value)(2^(-1))*(sum([GF(p_value)(bernoulli_polynomial(ZZ(i),2))*log_ell0(GF(ell0_value),i) for i in range(ell0_value) if i!=0]));
assert zeta_MT == GF(p_value)(2^(-1))*sum([(GF(p_value)(i)^2)*log_ell0(GF(ell0_value),i) for i in range(ell0_value) if i!=0]);

#Find primes of characteristic ell0 and that give the correction extension on the a1-side
exponents_dictionary = {};
for prime_ideal in Kz_abs.primes_above(ell0_value):
    distinguished_prime_ell0=prime_ideal;
    distinguished_prime_ell0_index = primes_over_ell0[prime_ideal];

    #Local method adjustment on a1-side for condition at ell0
    a1_c1_correction = local_correction_Ell0(distinguished_prime_ell0,D1_gamma_free*zeta_Kz_correction,a0_element);
    if a0_correction_exponent[0] == a1_c1_correction:
        assert distinguished_prime_ell0_index in distinguished_primes_ell0_a1_side, "problem with a0 adjustment using local method";

        #Local method adjustment on b2-side for condition at ell0
        adjustment_exponent_ell0_b2= local_correction_Ell0(distinguished_prime_ell0,b2_finite_flat,Kz_abs(ell0_value));
        exponents_dictionary[distinguished_prime_ell0_index] = [a1_c1_correction,adjustment_exponent_ell0_b2];

#Check that the correction exponents agree for entire distinguished orbit of primes
for exp_index in range(len(exponents_dictionary)):
    assert list(exponents_dictionary.values())[0] == list(exponents_dictionary.values())[exp_index], "error with correction exponents not matching";

#Check that distinguished orbit of primes matches on a(1) and b(2)-sides
assert set(distinguished_primes_ell0_a1_side) == set(exponents_dictionary.keys()), "distinguished orbit of primes doesn't match on a(1) and b(2)-sides";
distinguished_prime_ell0_b2_side = distinguished_primes_ell0_a1_side[0];

#print(f'\tThe distinguished prime at ell0={ell0_value} is the ideal numbered {distinguished_prime_ell0_b2_side}');
print(f'\tThe ell0 adjustment exponent is {adjustment_exponent_ell0_b2}, computed using the local method with the ideal numbered {distinguished_prime_ell0_b2_side}');

#Test the splitting behavior of each prime of characteristic ell0 on b2 side
beta_zero_check = None;
#print(f'Computing splitting behavior on b2-side using local method with ideal numbered {distinguished_prime_ell0_b2_side}:');
b2_adjusted=b2_finite_flat*(ell0_value)^exponents_dictionary[distinguished_prime_ell0_b2_side][1];
if p_power_check(b2_adjusted):
    print("\tb2_adjusted is computed and p-power free\n");
    b2_adjusted_free = b2_adjusted;
else:
    b2_adjusted_free = p_power_adjustment(b2_adjusted);
    if p_power_check(b2_adjusted_free):
        print("\tb2_adjusted is computed and p-power free\n");
    else:
        print("\tFurther adjustment is needed\n");

assert p_power_check(b2_adjusted_free), "b2_adjusted is not p-power free --  need to adjust";

b2_ramified_primes=[];
b2_split_primes=[];
b2_inert_primes=[];

for prime_ideal in Kz_abs.primes_above(ell0_value):
    if primes_over_ell0[prime_ideal] in exponents_dictionary:
        FF_ell0_check = prime_ideal.residue_field();
        if prime_ideal.divides(ideal(b2_adjusted_free)):
            b2_ramified_primes = b2_ramified_primes + [primes_over_ell0[prime_ideal]];
        elif Ell0_splits_check(prime_ideal,b2_adjusted_free):
            b2_split_primes = b2_split_primes + [primes_over_ell0[prime_ideal]];
            if primes_over_ell0[prime_ideal]== distinguished_prime_ell0_b2_side:
                beta_zero_check = True;
        else:
            b2_inert_primes = b2_inert_primes + [primes_over_ell0[prime_ideal]];

if beta_zero_check:
    assert galois_orbit_check(b2_split_primes), "galois split"
    print(f"Condition (ii): b2 splits at the distinguished prime at ell0={ell0_value}, so beta=0");
else:
    assert galois_orbit_check(b2_ramified_primes), "galois ramified"
    print(f"Condition (ii): b2 ramifies at the distinguished prime at ell0={ell0_value}, so beta does not equal 0");
