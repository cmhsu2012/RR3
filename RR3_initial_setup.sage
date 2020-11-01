#Define cyclotomic field
Qz=CyclotomicField(5);
zeta_Qz=Qz.0;
assert zeta_Qz^5==1 and zeta_Qz !=1, "zeta_Qz is not a primitive root of unity";

#Construct Kz_abs as splitting field of x^5-ell1
R.<x>=ZZ[];
K_def_poly = x^5-ell1_value
K_field.<sqrt_ell1>=NumberField(K_def_poly)
Kz_split=K_def_poly.splitting_field('kz_splitting')
Kz_poly = Kz_split.optimized_representation()[0].absolute_polynomial()
Kz_abs.<kz_gen> = NumberField(Kz_poly);

#Construct Kz as a relative extensions of Qz and K
from_Qz_to_Kz_abs = Qz.Hom(Kz_abs).an_element()
assert from_Qz_to_Kz_abs(1)!=0 #This checks that our chosen homomorphism is not trivial

from_K_to_Kz_abs = K_field.Hom(Kz_abs).an_element()
assert from_K_to_Kz_abs(1)!=0 #This checks that our chosen homomorphism is not trivial

Kz_over_Qz.<kz_qz>=Kz_abs.relativize(from_Qz_to_Kz_abs)
Kz_over_K.<kz_k>=Kz_abs.relativize(from_K_to_Kz_abs)

Kz_over_Qz_to_Kz_abs, Kz_abs_to_Kz_over_Qz = Kz_over_Qz.structure()
Kz_over_K_to_Kz_abs, Kz_abs_to_Kz_over_K = Kz_over_K.structure()

#Choose zeta_p in Kz
zeta_Kz = from_Qz_to_Kz_abs(zeta_Qz)

#Choose a generator of the Galois group of Kz/Qz
sigma = None
for aut in Kz_abs.galois_group():
    if aut(zeta_Kz) == zeta_Kz and not aut.multiplicative_order()==1:
        sigma = aut
        break
assert sigma.multiplicative_order()==5;

#Choose a generator of the Galois group of Kz/K
delta = None
for aut in Kz_abs.galois_group():
    if aut(from_K_to_Kz_abs(K_field.0))==from_K_to_Kz_abs(K_field.0) and aut(zeta_Kz)!=zeta_Kz and aut(zeta_Kz) == zeta_Kz^2:
        delta = aut
        break
assert delta.multiplicative_order()==4, "no";

#Label choice of generator for delta action in terms of exponent
for i in GF(5):
    if delta(zeta_Kz) == zeta_Kz^i:
        delta_generator = i
        break
assert delta_generator==2

#Construct unit groups
UKz_abs=UnitGroup(Kz_abs,proof=False)
assert UKz_abs.rank()==9

UKz_over_Qz = UnitGroup(Kz_over_Qz, proof=False)
assert UKz_over_Qz.rank()==9

UQz = UnitGroup(Qz)

#Fix a numbering system for primes over ell0
primes_over_ell0={}
for index in range(20):
    primes_over_ell0[Kz_abs.primes_above(11)[index]]=index
assert len(primes_over_ell0)==20

#Construct S-unit groups for ell0
S_Qz_ell0 =Qz.primes_above(11)
SUQz_ell0 =UnitGroup(Qz,S=tuple(S_Qz_ell0))

S_Kz_over_Qz_ell0=Kz_over_Qz.primes_above(11)
SUKz_over_Qz_ell0=UnitGroup(Kz_over_Qz,proof=False,S=tuple(S_Kz_over_Qz_ell0))

S_Kz_abs_ell0=Kz_abs.primes_above(11)
SUKz_abs_ell0=UnitGroup(Kz_abs,proof=False,S=tuple(S_Kz_abs_ell0))

#Check ramification of p in K/Q
p_ramification=None; #this Boolean is "true" if 5 wildly ramifies in K/Q and false otherwise
if 25.divides(ell1_value^4-1):
    print("5 tamely ramifies in K/Q")
    p_ramification=false;
else:
    print("5 wildly ramifies in K/Q")
    p_ramification=true;

#Find p-powers (in the wild ramification case)
if p_ramification:
    P_Kz=Kz_abs.primes_above(5)[0]
    P_Kz_25 = P_Kz^25
    for _ in range(100):
        possible_gen = P_Kz.random_element()
        if possible_gen.norm().ord(5)==1:
            prime_elt_Kz_local=possible_gen
            print("local generator found in wild ramification case")
            break
    assert prime_elt_Kz_local.ord(P_Kz)==1

    coeff_space = VectorSpace(GF(5),5)

    elts_all_res = [sum([Kz_abs(a)*prime_elt_Kz_local^i for i,a in enumerate(v)]) for v in coeff_space]
    elts_all_res_pth_powers = [a^5 for a in elts_all_res]
    elts_all_res_pth_powers_mod = set([a.mod(P_Kz_25) for a in elts_all_res_pth_powers])

    print("p-powers computed in wild ramification case")