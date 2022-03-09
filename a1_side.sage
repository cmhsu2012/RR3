print(f'\nComputing a1 data:');

#Adjustment for condition at p
roots_of_unity_Kz=[zeta_Kz^i for i in range(p_value)];
roots_of_unity_Kz_correction = None;

if p_ramification:
    roots_of_unity_Kz_correction=[z for z in roots_of_unity_Kz if wild_p_unramified_check(z,D1_gamma_free)];
else:
    roots_of_unity_Kz_correction=[z for z in roots_of_unity_Kz if tame_p_unramified_check(z,D1_gamma_free)];

assert len(roots_of_unity_Kz_correction)==1, "error with adjustment for condition at p";

zeta_Kz_correction =roots_of_unity_Kz_correction[0];

for i in range(p_value):
    if zeta_Kz_correction == zeta_Kz^i:
        zeta_correction = i;
        zeta_Qz_correction =zeta_Qz^zeta_correction;
        break;
print(f'\tThe zeta_p adjustment exponent is {zeta_correction}');

#For each possible correction at ell0, test for splitting on a1 side
a0_correction_exponent = [];
distinguished_primes_ell0_a1_side = [];

for i in range(p_value):
    a1_candidate = D1_gamma_free*zeta_Kz_correction*a0_element^i;
    if p_power_check(a1_candidate):
        a1_candidate_free = a1_candidate;
    else:
        a1_candidate_free = p_power_adjustment(a1_candidate);
        assert p_power_check(a1_candidate_free), "Further adjustment is needed to make a1_candidate p-power free for a_0 exponent {i}"
    ramified_primes_ell0_a1_side, split_primes_ell0_a1_side, inert_primes_ell0_a1_side= Ell0_splitting_patterns(a1_candidate_free);
    assert len(ramified_primes_ell0_a1_side)==(p_value-1)^2, "ram";
    assert len(split_primes_ell0_a1_side)==(p_value-1) or len(inert_primes_ell0_a1_side)==(p_value-1), "inert or split";
    assert galois_orbit_check(split_primes_ell0_a1_side), "galois split";
    assert galois_orbit_check(inert_primes_ell0_a1_side), "galois inert";
    if len(split_primes_ell0_a1_side)!=0:
        a0_correction_exponent = a0_correction_exponent +[i];
        distinguished_primes_ell0_a1_side = distinguished_primes_ell0_a1_side + split_primes_ell0_a1_side;

assert len(a0_correction_exponent)==1, "len";

print(f'\tThe a0 adjustment exponent is {a0_correction_exponent[0]}');
print(f'\tThe distinguished Galois orbit of primes at ell0={ell0_value} is {distinguished_primes_ell0_a1_side}');
print(f'\tWe take the distinguished prime at ell0={ell0_value} to be the ideal numbered {distinguished_primes_ell0_a1_side[0]}');
print("\tWe have arranged for alpha=0\n");

#Adjustment for condition at ell0
adjustment_a1_Kz =mod_p_powers(SUKz_abs_ell0,zeta_Kz_correction*a0_element^(a0_correction_exponent[0]));
adjustment_a1_Qz = mod_p_powers(SUQz_ell0,zeta_Qz_correction*a0_cyc^(a0_correction_exponent[0]));
assert mod_p_powers(SUKz_abs_ell0,from_Qz_to_Kz_abs(adjustment_a1_Qz))==mod_p_powers(SUKz_abs_ell0,adjustment_a1_Kz);

#Check if a1 splits at ell1
if ell1_split_check(D1_gamma_free*adjustment_a1_Kz):
    print(f'Condition (i): a1 splits at ell1={ell1_value}\n')
else:
    print(f'Condition (i): a1 does not split at ell1={ell1_value}\n')

