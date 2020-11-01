#Compute xi in Kz
adjustment_a1_vector=vector(list(SUQz_ell0.log(SUQz_ell0(adjustment_a1_Qz))));
xi_vector=relative_norm_matrix.solve_right(adjustment_a1_vector.column()).column(0);
xi_Kz_over_Qz=SUKz_over_Qz_ell0.exp(xi_vector);
xi_abs = mod_p_powers(SUKz_abs_ell0,Kz_over_Qz_to_Kz_abs(mod_p_powers(SUKz_over_Qz_ell0,xi_Kz_over_Qz)))

#Test that relative norm of xi is correct
assert SUQz_ell0.log(mod_p_powers(SUQz_ell0,mod_p_powers(SUKz_over_Qz_ell0,xi_Kz_over_Qz).relative_norm()))==SUQz_ell0.log(SUQz_ell0(adjustment_a1_Qz))

#Compute Kolyvagin derivative of xi
Dxi=mod_p_powers(SUKz_abs_ell0,D_1_operator(sigma,xi_abs))

#Define candidate for b2
b2_candidate = mod_p_powers(SUKz_abs_ell0,(D2_gamma_free*Dxi)^(-2)*(D1_gamma_free*adjustment_a1_Kz)^(-1))

#Check if b2_candidate is p-power free
if p_power_check(b2_candidate):
    print("b2_candidate is p-power free.")
    b2_candidate_free = b2_candidate
elif p_power_10_check(b2_candidate):
    print("D2_gamma is not p-power free -- adjusting now.")
    b2_candidate_free = p_power_adjustment_10(b2_candidate)
    assert p_power_check(b2_candidate_free)
    print("We've adjusted so that b2_candidate is p-power free.")
else:
    print("b2_candidate is not p-power free and further adjustment is needed.")

#Adjustment for condition at p
if not 5.divides(b2_candidate_free.norm()):
    print("No adjustment needed for condition at p.")
else:
    print("Adjustment needed for condition at p. Check this.")

#Mazur-Tate derivative (in F_p)
zeta_MT = GF(5)(3)*(sum([GF(5)(bernoulli_polynomial(ZZ(i),2))*log_ell0(GF(11),i) for i in range(11) if i!=0]))
print('zeta_MT is equal to',GF(5)(zeta_MT))
assert zeta_MT == GF(5)(3)*sum([(GF(5)(i)^2)*log_ell0(GF(11),i) for i in range(11) if i!=0])

#Find primes of characteristic ell0 and that give the correction extension on the a1-side
exponents_dictionary = {}
for prime_ideal in Kz_abs.primes_above(11):
    distinguished_prime_ell0=prime_ideal
    distinguished_prime_ell0_index = primes_over_ell0[prime_ideal]

    #Local method adjustment on a1-side for condition at ell0
    a1_c1_correction = local_correction_Ell0(distinguished_prime_ell0,D1_gamma_free*zeta_Kz_correction,a0_element)
    if a0_correction_exponent[0] == a1_c1_correction:
        assert distinguished_prime_ell0_index in distinguished_primes_ell0_a1_side

        #Local method adjustment on b2-side for condition at ell0
        adjustment_exponent_ell0_b2= local_correction_Ell0(distinguished_prime_ell0,b2_candidate_free,Kz_abs(11))
        exponents_dictionary[distinguished_prime_ell0_index] = [a1_c1_correction,adjustment_exponent_ell0_b2]
        print("For ideal numbered", distinguished_prime_ell0_index, "the correction exponents are", exponents_dictionary[distinguished_prime_ell0_index])

#Test the splitting behavior of each prime of characteristic ell0 on b2 side
beta_zero_check = None
for ideal_index in exponents_dictionary:
    b2_adjusted=mod_p_powers(SUKz_abs_ell0,b2_candidate_free*11^exponents_dictionary[ideal_index][1])
    assert p_power_check(b2_adjusted), "b2_adjusted is not p-power free --  need to adjust"

    print("Computing splitting behavior on b2-side using local method with ideal numbered", ideal_index)
    for prime_ideal in Kz_abs.primes_above(11):
        if primes_over_ell0[prime_ideal] in exponents_dictionary:
            FF_ell0_check = prime_ideal.residue_field()
            if prime_ideal.divides(ideal(b2_adjusted)):
                print("\t The ideal numbered", primes_over_ell0[prime_ideal],"ramifies.")
            elif Ell0_splits_check(prime_ideal,b2_adjusted):
                print("\t The ideal numbered", primes_over_ell0[prime_ideal], "splits.")
                if primes_over_ell0[prime_ideal]== ideal_index:
                    beta_zero_check = True
            else:
                print("\t The ideal numbered", primes_over_ell0[prime_ideal], "is inert.")

if beta_zero_check:
    print("We can arrange for beta=0.")
else:
    print("We cannot arrange for beta=0.")