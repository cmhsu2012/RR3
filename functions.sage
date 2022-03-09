#This function tests if the given Kummer extensions are equal
def Kummer_equal_check(Field_1,Field_2):
    return Field_1.optimized_representation()[0].absolute_polynomial()==Field_2.optimized_representation()[0].absolute_polynomial();

#Define omega function in terms of chosen delta_generator
def omega(i):
    return ZZ(delta_generator^i);

#Mod p-powers in each entry of a SU-vector
def mod_p_powers(SU,alpha):
    #takes alpha in SU and returns the element x of SU such that x=alpha mod SU^p, and such that x is integral and p-power-free
    M=VectorSpace(GF(p_value),SU.ngens());
    return SU.exp(M(SU.log(alpha)));

#Check that an element is p-power free
def p_power_check(alpha):
    for i in range(len(list(ideal(alpha).factor()))):
        if list(ideal(alpha).factor())[i][1]%p_value==0:
            return False;
    return True; days

#Lists the prime factorization of an ideal:
def ideal_prime_factorization(ideal_alpha):
    p_power_ideals = [];
    p_power_ideals_exp = [];
    for i in range(len(list(ideal_alpha.factor()))):
        p_power_ideals = p_power_ideals + [list(ideal_alpha.factor())[i][0]];
        p_power_ideals_exp = p_power_ideals_exp + [list(ideal_alpha.factor())[i][1]];
    return(p_power_ideals,p_power_ideals_exp);

#p_power_adjustment helper function
def p_power_adjustment_helper(current_exp_index,list_of_ideals, list_of_exponents):
    current_exp = list_of_exponents[current_exp_index];
    exp_residue_mod_n = current_exp%Kz_abs_clgp_exp;
    modulus_p_residue = p_value*exp_residue_mod_n;
    modulus_p_n = p_value*Kz_abs_clgp_exp;
    crt_output = crt([0,current_exp],[modulus_p_residue,modulus_p_n]);
    assert crt_output>0, "chinese remainder theorem output is negative";
    correction_exp = crt_output/modulus_p_residue;
    assert correction_exp in ZZ, "correction exponent is not an integer";
    
    boost_dict = {};
    
    #Adjustment for the exponent we are currently removing
    boost_dict[list_of_ideals[current_exp_index]]=(crt_output/modulus_p_n).floor();
    adjustment_ideal_part = list_of_ideals[current_exp_index]^exp_residue_mod_n;
    
    #Adjustment for the rest of the exponents 
    for j in range(len(list_of_ideals)):
        if j != current_exp_index:
            if list_of_exponents[j]%Kz_abs_clgp_exp ==0:
                1==1
            elif (list_of_exponents[j])>=(p_value*correction_exp*(list_of_exponents[j]%Kz_abs_clgp_exp)):
                adjustment_ideal_part = adjustment_ideal_part*(list_of_ideals[j]^(list_of_exponents[j]%Kz_abs_clgp_exp));
            else:
                diff_to_boost = (((p_value*correction_exp*(list_of_exponents[j]%Kz_abs_clgp_exp)-list_of_exponents[j]))/modulus_p_n).ceil();
                boost_dict[list_of_ideals[j]]=diff_to_boost;
                adjustment_ideal_part = adjustment_ideal_part*(list_of_ideals[j]^(list_of_exponents[j]%Kz_abs_clgp_exp));
    
    boost_ideal_adjustment = 1;
    for boost_ideal in list(boost_dict.keys()):
        boost_ideal_adjustment = boost_ideal_adjustment*(boost_ideal^Kz_abs_clgp_exp).gens_reduced(proof=False)[0]^(p_value*boost_dict[boost_ideal]);
    
    if boost_ideal_adjustment != 1:
        assert adjustment_ideal_part.is_principal(proof=False), 'adjustment_ideal_part is not principal';
        alpha_adjustment = boost_ideal_adjustment*(adjustment_ideal_part.gens_reduced(proof = False)[0]^(-p_value*correction_exp));
    else:
        assert adjustment_ideal_part.is_principal(proof=False), 'adjustment_ideal_part is not principal';
        alpha_adjustment = (adjustment_ideal_part.gens_reduced(proof = False)[0]^(-p_value*correction_exp));
    return alpha_adjustment;

def p_power_adjustment(alpha):
    alpha_ideal = ideal(alpha);
    p_power_ideals = ideal_prime_factorization(alpha_ideal)[0];
    p_power_ideals_exp=ideal_prime_factorization(alpha_ideal)[1];
        
    first_adjustment= 1
    for j in range(len(p_power_ideals)):
        if p_power_ideals_exp[j]>=p_value*Kz_abs_clgp_exp:
            first_adjustment = first_adjustment*((p_power_ideals[j]^Kz_abs_clgp_exp).gens_reduced(proof=false)[0])^(-((p_power_ideals_exp[j]/(p_value*Kz_abs_clgp_exp)).floor())*p_value)
    alpha_first_adj = alpha*first_adjustment      

    alpha_ideal_first_adj = ideal(alpha_first_adj);
    assert alpha_ideal_first_adj.is_principal(proof=False), "alpha ideal not principal"
    if p_power_check(alpha_ideal_first_adj):
        alpha_free = alpha_first_adj
        if alpha_free in SUKz_abs_ell0 and alpha in SUKz_abs_ell0:
            assert mod_p_powers(SUKz_abs_ell0,alpha)==mod_p_powers(SUKz_abs_ell0,alpha_free), "output is not equivalent to alpha"
        return alpha_free;    
    else:
        p_power_ideals = ideal_prime_factorization(alpha_ideal_first_adj)[0];
        p_power_ideals_exp = ideal_prime_factorization(alpha_ideal_first_adj)[1]
        
        for j in range(len(p_power_ideals)):
            if p_power_ideals_exp[j]%p_value == 0:
                current_exp = j
                break
        
        new_alpha = alpha_first_adj*p_power_adjustment_helper(current_exp,p_power_ideals, p_power_ideals_exp)
        if p_power_check(new_alpha):
            alpha_free = new_alpha
            if alpha_free in SUKz_abs_ell0 and alpha in SUKz_abs_ell0:
                assert mod_p_powers(SUKz_abs_ell0,alpha)==mod_p_powers(SUKz_abs_ell0,alpha_free), "output is not equivalent to alpha"
            return alpha_free;
        else:
            next_alpha = p_power_adjustment(new_alpha)
            alpha_free = next_alpha
            assert p_power_check(alpha_free), "alpha_free not p-power free";
            if alpha_free in SUKz_abs_ell0 and alpha in SUKz_abs_ell0:
                assert mod_p_powers(SUKz_abs_ell0,alpha)==mod_p_powers(SUKz_abs_ell0,alpha_free), "output is not equivalent to alpha"
            return alpha_free;    

#Function to check omega action
def omega_action_test(alpha,omega_exp_test):
    return SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,delta(alpha)))==SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,(alpha^omega_exp_test)));

#Kolyvagin derivative operators
def D_0_operator(aut,alpha):
    return prod([(aut^(i))(alpha) for i in range(p_value)]);

def D_1_operator(aut,alpha):
    return prod([((aut^(i))(alpha))^(i) for i in range(p_value)]);

def D_2_operator(aut, alpha):
    return prod([((aut^(i))(alpha))^(binomial(i,2)) for i in range(p_value)]);

#This function checks if a prime Ell0 splits in a Kummer extension Kz(sqrt[p]{alpha})/Q
def Ell0_splits_check(Ell0,alpha):
    #returns "true" if the given prime Ell0 splits in Kz(sqrt[p]{alpha})/Kz, and "false" otherwise
    F_Ell0 = Ell0.residue_field();
    g=F_Ell0.primitive_element();
    if not p_value.divides(F_Ell0(alpha).log(g)): #This is checking if Ell0 is inert in Kz(sqrt[p]{alpha})/Kz
        return False;
    return True;

#This function sorts the primes of characteristic ell_0 into ramified, split, and inert primes in a Kummer extension Kz(sqrt[p]{alpha})/Q
def Ell0_splitting_patterns(alpha):
    ramified_ideals_check = [];
    split_ideals_check = [];
    inert_ideals_check =[];
    for prime_ideal in Kz_abs.primes_above(ell0_value):
        if prime_ideal.divides(ideal(alpha)):
            ramified_ideals_check = ramified_ideals_check+[primes_over_ell0[prime_ideal]];
        elif Ell0_splits_check(prime_ideal,alpha):
            split_ideals_check = split_ideals_check + [primes_over_ell0[prime_ideal]];
        else:
            inert_ideals_check = inert_ideals_check + [primes_over_ell0[prime_ideal]];

    return ramified_ideals_check, split_ideals_check, inert_ideals_check;

#Test that the prime of characteristic p are unramified in each extension when p is wildly ramified in K/Q
def wild_p_unramified_check(alpha,beta):
    #returns "true" if prime over p in Kz is unramified in Kz(sqrt[p]{beta*alpha})/Kz, and "false" otherwise
    assert p_power_check(beta*alpha);
    if (beta*alpha).mod(P_Kz_psquared) in elts_all_res_pth_powers_mod:
        return True;
    return False;

#Test that all primes of characteristic p are unramified in each extension when p is tamely ramified in K/Q
def tame_p_unramified_check(alpha,beta):
    #returns "true" if all primes over 5 in Kz are unramified in Kz(sqrt[p]{beta*alpha})/Kz, and "false" otherwise
    assert p_power_check(beta*alpha);
    for P in P_Kz_all:
        if not ((beta*alpha)^(p_value-1)-1) in P^p_value: #This is checking if P is ramified in Kz(sqrt[p]{beta*alpha})/Kz
            return False;
    return True;

def ell1_split_check(alpha):
    #returns "true" if all primes over ell1 in Kz are split in Kz(sqrt[p]{alpha})/Kz, and "false" otherwise
    unramified_primes_ell1=[prime for prime in Kz_abs.primes_above(ell1_value) if not prime.divides(ideal(alpha))];
    F_ell1s=[prime.residue_field() for prime in unramified_primes_ell1]
    if len(unramified_primes_ell1)==len(Kz_abs.primes_above(ell1_value)):
        if {p_value.divides(F_ell1(alpha).log(F_ell1.primitive_element())) for F_ell1 in F_ell1s}=={True}:
            return true;
        else:
            return false;
    else:
        return false;

#Test that all primes of characteristic p are unramified in each extension when p is tamely ramified in K/Q
def p_split_in_a1_check_tame(alpha):
    #returns "true" if all primes over 5 in Kz are split in Kz(sqrt[p]{alpha})/Kz, and "false" otherwise
    assert p_power_check(alpha);
    for P in P_Kz_all:
        if not (alpha^(p_value-1)-1) in P^(p_value+1): #This is checking if P is ramified in Kz(sqrt[p]{beta*alpha})/Kz
            return False;
    return True;

#Test that the prime of characteristic p is completely split in a Kummer extension Kz(sqrt[p]{alpha})/Q
def p_split_in_a1_check(alpha):
    #returns "true" if prime over p in Qz is completely split in Qz(sqrt[p]{alpha})/Qz
    assert p_power_check(alpha);
    if (alpha).mod(P_Kz_psquared_plusone) in elts_all_res_pth_powers_mod_plusone:
        return True;
    return False;

#Test that the prime of characteristic p is completely split in Lz/Qz
def p_split_in_Lz_check(alpha):
    #returns "true" if prime over p in Qz is completely split in Qz(sqrt[p]{alpha})/Qz
    assert p_power_check(alpha);
    if (alpha).mod(P_Qz_pplusone) in elts_all_res_pth_powers_mod_Qz:
        return True;
    return False;

#Function that tests if a collection of prime ideals are in a Galois orbit under the action of delta
def galois_orbit_check(prime_ideal_indices):
    for ideal_index in prime_ideal_indices:
        if primes_over_ell0[delta(Kz_abs.primes_above(ell0_value)[ideal_index])] not in prime_ideal_indices:
                return False;
    return True;

#Function that takes in an element, coerces the element into an appropriate finite field, computes the log (with 2 fixed as a generator), and then returns the result modulo 5
def log_ell0(FF,n):
    g=FF.primitive_element();
    return GF(p_value)(FF(n).log(g));

#Function that computes local correction (with correct coefficient)
def local_correction_Ell0(Ell0,alpha,beta):
    prime_ell0_elt_Kz_local=Ell0.0;
    assert prime_ell0_elt_Kz_local.ord(Ell0)==1;
    assert prime_ell0_elt_Kz_local ==ell0_value;
    FF_ell0 = Ell0.residue_field();
    return -(GF(p_value)(6)*zeta_MT*(alpha).ord(Ell0)+ log_ell0(FF_ell0,alpha/(prime_ell0_elt_Kz_local^((alpha).ord(Ell0)))))*(GF(p_value)(6)*zeta_MT*(beta).ord(Ell0)+ log_ell0(FF_ell0,beta/(prime_ell0_elt_Kz_local^((beta).ord(Ell0)))))^(-1);