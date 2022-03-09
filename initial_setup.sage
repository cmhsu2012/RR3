#Define cyclotomic field
Qz=CyclotomicField(p_value);
zeta_Qz=Qz.0;
assert zeta_Qz^p_value==1 and zeta_Qz !=1, "zeta_Qz is not a primitive root of unity";

#Construct Kz_abs as splitting field of x^5-ell1
R.<x>=ZZ[];
K_def_poly = x^p_value-ell1_value;
K_field.<sqrt_ell1>=NumberField(K_def_poly);
Kz_split=K_def_poly.splitting_field('kz_splitting');
Kz_poly = Kz_split.optimized_representation()[0].absolute_polynomial();
Kz_abs.<kz_gen> = NumberField(Kz_poly);

print("\tKz_abs constructed");

#Check class groups
Kz_abs_clgp_exp = Kz_abs.class_group(proof=False).exponent();
K_field_clgp_exp = K_field.class_group(proof=False).exponent(); 
assert Kz_abs_clgp_exp == K_field_clgp_exp, "class numbers are not equal";
assert gcd(Kz_abs_clgp_exp,p_value)==1, "class number is not coprime to p";
print(f'\tThe class group of Kz_abs has exponent {Kz_abs_clgp_exp}');

#Construct Kz as a relative extensions of Qz and K
from_Qz_to_Kz_abs = Qz.Hom(Kz_abs).an_element();
assert from_Qz_to_Kz_abs(1)!=0; #This checks that our chosen homomorphism is not trivial

from_K_to_Kz_abs = K_field.Hom(Kz_abs).an_element();
assert from_K_to_Kz_abs(1)!=0; #This checks that our chosen homomorphism is not trivial

Kz_over_Qz.<kz_qz>=Kz_abs.relativize(from_Qz_to_Kz_abs);
Kz_over_K.<kz_k>=Kz_abs.relativize(from_K_to_Kz_abs);

Kz_over_Qz_to_Kz_abs, Kz_abs_to_Kz_over_Qz = Kz_over_Qz.structure();
Kz_over_K_to_Kz_abs, Kz_abs_to_Kz_over_K = Kz_over_K.structure();

#Choose zeta_p in Kz
zeta_Kz = from_Qz_to_Kz_abs(zeta_Qz);

#Choose a generator of the Galois group of Kz/Qz
sigma = None;
for aut in Kz_abs.galois_group():
    if aut(zeta_Kz) == zeta_Kz and not aut.multiplicative_order()==1:
        sigma = aut;
        break;
assert sigma.multiplicative_order()==p_value, "error with order of sigma";

#Choose a generator of the Galois group of Kz/K
delta = None;
for aut in Kz_abs.galois_group():
    if aut(from_K_to_Kz_abs(K_field.0))==from_K_to_Kz_abs(K_field.0) and aut(zeta_Kz)!=zeta_Kz and aut.multiplicative_order()==(p_value-1):
        delta = aut;
        break;
assert delta.multiplicative_order()== (p_value-1), "error fixing delta generator";

#Label choice of generator for delta action in terms of exponent
for i in GF(p_value):
    if delta(zeta_Kz) == zeta_Kz^i:
        delta_generator = i;
        break;

#Label corresponding generator of the Galois group of Qz/Q
delta_Qz = None;
for aut in Qz.galois_group():
    if aut(zeta_Qz) == zeta_Qz^delta_generator:
        delta_Qz = aut;
        break;
assert delta_Qz.multiplicative_order()==(p_value-1);

print("\tGalois elements chosen");

#Construct unit groups
UKz_abs=UnitGroup(Kz_abs,proof=False);
assert UKz_abs.rank()==(p_value*(p_value-1)/2)-1, "error with unit group rank abs";

UKz_over_Qz = UnitGroup(Kz_over_Qz, proof=False);
assert UKz_over_Qz.rank()==(p_value*(p_value-1)/2)-1, "error with unit group rank rel Qz";

UQz = UnitGroup(Qz);

#Fix a numbering system for primes over ell0
primes_over_ell0={};
for index in range(p_value*(p_value-1)):
    primes_over_ell0[Kz_abs.primes_above(ell0_value)[index]]=index;
assert len(primes_over_ell0)==p_value*(p_value-1), "error with ell0 primes dictionary";

print(f'\tThe size of {ell0_value}-primes dictionary is {len(primes_over_ell0)}');

#Construct S-unit groups for ell0
S_Qz_ell0 =Qz.primes_above(ell0_value);
SUQz_ell0 =UnitGroup(Qz,S=tuple(S_Qz_ell0));

S_Kz_over_Qz_ell0=Kz_over_Qz.primes_above(ell0_value);
SUKz_over_Qz_ell0=UnitGroup(Kz_over_Qz,proof=False,S=tuple(S_Kz_over_Qz_ell0));

S_Kz_abs_ell0=Kz_abs.primes_above(ell0_value);
SUKz_abs_ell0=UnitGroup(Kz_abs,proof=False,S=tuple(S_Kz_abs_ell0));

#S-units in K/Q
S_Kz_over_K_ell0=Kz_over_K.primes_above(ell0_value);
SUKz_over_K_ell0=UnitGroup(Kz_over_K,proof=False,S=tuple(S_Kz_over_K_ell0));

print("\tS-unit groups constructed");

#Check ramification of p in K/Q
p_ramification=None; #this Boolean is "true" if p wildly ramifies in K/Q and false otherwise
if p_value.divides(K_field.primes_above(p_value)[0].ramification_index()):
    print(f'\tp={p_value} wildly ramifies in K/Q');
    p_ramification=True;
else:
    print(f'\tp={p_value} tamely ramifies in K/Q');
    p_ramification=False;

#Find p-powers mod p^(p*p) in Kz (in the wild ramification case)
if p_ramification:
    P_Kz=Kz_abs.primes_above(p_value)[0];
    P_Kz_psquared = P_Kz^(p_value*p_value);
    for _ in range(100):
        possible_gen = P_Kz.random_element();
        if possible_gen.norm().ord(p_value)==1:
            prime_elt_Kz_local=possible_gen;
            break;
    assert prime_elt_Kz_local.ord(P_Kz)==1;

    coeff_space = VectorSpace(GF(p_value),p_value);

    elts_all_res = [sum([Kz_abs(a)*prime_elt_Kz_local^i for i,a in enumerate(v)]) for v in coeff_space];
    elts_all_res_pth_powers = [a^p_value for a in elts_all_res];
    elts_all_res_pth_powers_mod = set([a.mod(P_Kz_psquared) for a in elts_all_res_pth_powers]);

else:
    P_Kz_all = Kz_abs.primes_above(p_value);

#Find p-powers in Qz
P_Qz=Qz.primes_above(p_value)[0];
P_Qz_pplusone = P_Qz^(p_value+1);
prime_elt_Qz_local=(1-zeta_Qz);
assert prime_elt_Qz_local.ord(P_Qz)==1;

coeff_space_Qz = VectorSpace(GF(p_value),2);

elts_all_res_Qz = [sum([Qz(a)*prime_elt_Qz_local^i for i,a in enumerate(v)]) for v in coeff_space_Qz];
elts_all_res_pth_powers_Qz = [a^p_value for a in elts_all_res_Qz];
elts_all_res_pth_powers_mod_Qz = set([a.mod(P_Qz_pplusone) for a in elts_all_res_pth_powers_Qz]);