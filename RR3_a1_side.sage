#Adjustment for condition at p
roots_of_unity_Kz=[zeta_Kz^i for i in range(5)];
roots_of_unity_Kz_correction = None;

if p_ramification:
    roots_of_unity_Kz_correction=[z for z in roots_of_unity_Kz if wild_p_unramified_check(z,D1_gamma_free)]
else:
    roots_of_unity_Kz_correction=[z for z in roots_of_unity_Kz if tame_p_unramified_check(z,D1_gamma_free)]

assert len(roots_of_unity_Kz_correction)==1

zeta_Kz_correction =roots_of_unity_Kz_correction[0];

for i in range(5):
    if zeta_Kz_correction == zeta_Kz^i:
        zeta_correction = i
        zeta_Qz_correction =zeta_Qz^zeta_correction
        break
print("The zeta_p_correction_exponent on the a1 side is",zeta_correction)

#Compute a0 (as an element in Qz)
expvectors=[tuple(l) for l in ProjectiveSpace(GF(5),5)]
all_ell0_unit_lines=[SUQz_ell0.exp(vector_exp) for vector_exp in expvectors]

delta_Qz = None
for aut in Qz.galois_group():
    if aut(zeta_Qz) == zeta_Qz^delta_generator:
        delta_Qz = aut;
        break;
assert delta_Qz.multiplicative_order()==4

possible_a0_Qz =[a0_Qz for a0_Qz in all_ell0_unit_lines if mod_p_powers(SUQz_ell0,delta_Qz(a0_Qz))==mod_p_powers(SUQz_ell0,a0_Qz^omega(1))];
#possible_a0_Qz is a list of representatives in SUQz for the set of omega-isotypic lines in SUQz \otimes F_p
print("the number of possible a0_Qz lines is",len(possible_a0_Qz));

#Possible a0 corrections from cyclotomic fields method
Q11cyc=CyclotomicField(11);
Qcyc.<b>=Qz.extension(Q11cyc.polynomial());
zeta11=Qcyc.zeta(11);
Qcycplus.<c>=Qz.extension((zeta11+zeta11^(-1)).minpoly());
Qcycplus_abs=Qcycplus.absolute_field('d')
s=Qz['s'].0;
possible_a0_cyc=[a0_test for a0_test in possible_a0_Qz if Qz.extension(s^5-a0_test,'a').absolute_field('a').is_isomorphic(Qcycplus_abs)];

assert len(possible_a0_cyc)==1

a0_cyc=possible_a0_cyc[0];
a0_element = mod_p_powers(SUKz_abs_ell0,from_Qz_to_Kz_abs(a0_cyc))
assert SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,a0_element^omega(1)))==SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,delta(a0_element)))

###Adjust so that a0_correction_candidates only have powers divisible by 10 -- this section is done manually (for now)
if ell1_value ==67:
    a0_exp_0_candidate = D1_gamma_free
    a0_exp_1_candidate = D1_gamma_free*a0_element^1*((Kz_abs.primes_above(11)[1]*Kz_abs.primes_above(11)[5]*Kz_abs.primes_above(11)[12]*Kz_abs.primes_above(11)[18])).gens_reduced(proof = False)[0]^(-5)
    a0_exp_2_candidate = D1_gamma_free*a0_element^2*(SUKz_abs_ell0.gens_values()[21])^(-5)
    a0_exp_3_candidate = D1_gamma_free*a0_element^3*(SUKz_abs_ell0.gens_values()[20])^(-5)*(SUKz_abs_ell0.gens_values()[27])^(-5)
    a0_exp_4_candidate = D1_gamma_free*a0_element^4*(SUKz_abs_ell0.gens_values()[17])^(-5)*(SUKz_abs_ell0.gens_values()[23])^(-5)
else:
    a0_exp_0_candidate = D1_gamma_free
    a0_exp_1_candidate = mod_p_powers(SUKz_abs_ell0,D1_gamma_free*zeta_Kz_correction*a0_element^1)
    a0_exp_2_candidate = mod_p_powers(SUKz_abs_ell0,D1_gamma_free*zeta_Kz_correction*a0_element^2)
    a0_exp_3_candidate = mod_p_powers(SUKz_abs_ell0,D1_gamma_free*zeta_Kz_correction*a0_element^3)
    a0_exp_4_candidate = mod_p_powers(SUKz_abs_ell0,D1_gamma_free*zeta_Kz_correction*a0_element^4)

a1_candidates_10=[a0_exp_0_candidate,a0_exp_1_candidate,a0_exp_2_candidate,a0_exp_3_candidate,a0_exp_4_candidate]

#For each possible correction at ell0, test for splitting on a1 side
a0_correction_exponent = []
distinguished_primes_ell0_a1_side = []

for i in range(5):
    print("Computing splitting patterns for a_0 exponent", i)
    if p_power_check(a1_candidates_10[i]):
        print("\t a1_candidate is p-power free.")
        a1_candidate_free = a1_candidates_10[i]
    elif p_power_10_check(a1_candidates_10[i]):
        print("\t a1_candidate is not p-power free -- adjusting now.")
        a1_candidate_free = p_power_adjustment_10(a1_candidates_10[i])
        assert p_power_check(a1_candidate_free)
        print("\t We've adjusted so that a1_candidate is p-power free.")
    else:
        print("\t a1_candidate is not p-power free and further adjustment is needed.")
    ramified_primes_ell0_a1_side, split_primes_ell0_a1_side, inert_primes_ell0_a1_side= Ell0_splitting_patterns(a1_candidate_free)
    print("\t ", split_primes_ell0_a1_side,inert_primes_ell0_a1_side)
    assert len(ramified_primes_ell0_a1_side)==16, "ram"
    assert len(split_primes_ell0_a1_side)==4 or len(inert_primes_ell0_a1_side)==4, "inert or split"
    assert galois_orbit_check(split_primes_ell0_a1_side), "galois split"
    assert galois_orbit_check(inert_primes_ell0_a1_side), "galois inert"
    if len(split_primes_ell0_a1_side)!=0:
        a0_correction_exponent = a0_correction_exponent +[i]
        distinguished_primes_ell0_a1_side = distinguished_primes_ell0_a1_side + split_primes_ell0_a1_side

assert len(a0_correction_exponent)==1, "len"

print("The a0_correction_exponent on the a1 side is", a0_correction_exponent[0])
print("The distinguished Galois orbit of primes on the a1 side is", distinguished_primes_ell0_a1_side)
print("We have arranged for alpha=0.")

#Adjustment for condition at ell0
adjustment_a1_Kz =mod_p_powers(SUKz_abs_ell0,zeta_Kz_correction*a0_element^(a0_correction_exponent[0]))
adjustment_a1_Qz = mod_p_powers(SUQz_ell0,zeta_Qz_correction*a0_cyc^(a0_correction_exponent[0]))
assert mod_p_powers(SUKz_abs_ell0,from_Qz_to_Kz_abs(adjustment_a1_Qz))==mod_p_powers(SUKz_abs_ell0,adjustment_a1_Kz)