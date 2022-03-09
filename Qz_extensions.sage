#Compute lines in S-unit group of Qz
expvectors=[tuple(l) for l in ProjectiveSpace(GF(p_value),SUQz_ell0.rank())];
all_ell0_unit_lines=[SUQz_ell0.exp(vector_exp) for vector_exp in expvectors];

#Compute a0 (as an element in Qz)
possible_a0_Qz =[a0_Qz for a0_Qz in all_ell0_unit_lines if mod_p_powers(SUQz_ell0,delta_Qz(a0_Qz))==mod_p_powers(SUQz_ell0,a0_Qz^omega(1))];

#Possible a0 corrections from cyclotomic fields method
Qell0cyc=CyclotomicField(ell0_value);
Qcyc.<b>=Qz.extension(Qell0cyc.polynomial());
zetaell0=Qcyc.zeta(ell0_value);
Qcycplus_abs=Qell0cyc.subfields(p_value)[0][0].composite_fields(Qz)[0];
s=Qz['s'].0;
possible_a0_cyc=[a0_test for a0_test in possible_a0_Qz if Qz.extension(s^p_value-a0_test,'a').absolute_field('a').is_isomorphic(Qcycplus_abs)];

assert len(possible_a0_cyc)==1, "error in number of possible values for a0";
print(f'\tS-unit for a0 computed');


a0_cyc=possible_a0_cyc[0];
a0_element = mod_p_powers(SUKz_abs_ell0,from_Qz_to_Kz_abs(a0_cyc));
assert SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,a0_element^omega(1)))==SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,delta(a0_element)));

#Compute c1_element (as an element in Qz)
possible_c1_Qz =[c1_Qz for c1_Qz in all_ell0_unit_lines if mod_p_powers(SUQz_ell0,delta_Qz(c1_Qz))==mod_p_powers(SUQz_ell0,c1_Qz^omega(2))];
possible_c1_cyc=[c1_test for c1_test in possible_c1_Qz if p_split_in_Lz_check(c1_test)];
assert len(possible_c1_cyc)==1, "error in number of possible values for c1";
print(f'\tS-unit for c1 computed');

c1_element = possible_c1_cyc[0];

#In the case of p=5, ell0=11, we can check our computation of c1_element using elliptic curves
if p_value ==5 and ell0_value==11:
    E=EllipticCurve("11a3");
    Ez=E.base_extend(Qz);
    Lz_abs.<lz_ell>=Ez.division_field(5);

    Lz_Kummer.<lz_kummer> = NumberField(Qz.extension(x^p_value-c1_element,'c1_element').absolute_polynomial());

    assert Kummer_equal_check(Lz_Kummer,Lz_abs);