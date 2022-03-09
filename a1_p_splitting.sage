#Check behavior of a(1) at p
a1_Kz =mod_p_powers(SUKz_abs_ell0,D1_gamma_free*zeta_Kz_correction*a0_element^(a0_correction_exponent[0]))
if not p_power_check(a1_Kz):
    a1_Kz_free = p_power_adjustment(a1_Kz)
else:
    a1_Kz_free = a1_Kz

#Find p-powers mod p^(p*p+1) in Kz (in the wild ramification case)
if p_ramification:
    P_Kz=Kz_abs.primes_above(p_value)[0];
    P_Kz_psquared_plusone = P_Kz^(p_value*p_value+1);
    for _ in range(100):
        possible_gen = P_Kz.random_element();
        if possible_gen.norm().ord(p_value)==1:
            prime_elt_Kz_local=possible_gen;
            #print("\tlocal generator found in wild ramification case");
            break;
    assert prime_elt_Kz_local.ord(P_Kz)==1;

    coeff_space_plusone = VectorSpace(GF(p_value),p_value+1);

    elts_all_res_plusone = [sum([Kz_abs(a)*prime_elt_Kz_local^i for i,a in enumerate(v)]) for v in coeff_space_plusone];
    elts_all_res_pth_powers_plusone = [a^p_value for a in elts_all_res_plusone];
    elts_all_res_pth_powers_mod_plusone = set([a.mod(P_Kz_psquared_plusone) for a in elts_all_res_pth_powers_plusone]);

    #print("\tp-powers computed in wild ramification case");

    if p_split_in_a1_check(a1_Kz_free):
        print(f"Difficulty check for finite-flat adjustment of b2: a1 splits at p={p_value}\n")
    else:
        print(f"Difficulty check for finite-flat adjustment of b2: a1 does not split at p={p_value}\n")
else:
    P_Kz_all = Kz_abs.primes_above(p_value);
    if p_split_in_a1_check_tame(a1_Kz_free):
        print(f"Difficulty check for finite-flat adjustment of b2: a1 splits at p={p_value}\n")
    else:
        print(f"Difficulty check for finite-flat adjustment of b2: a1 does not split at p={p_value}\n")
    