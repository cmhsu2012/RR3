#This function tests if the given Kummer extensions are equal
def Kummer_equal_check(Field_1,Field_2):
    return Field_1.optimized_representation()[0].absolute_polynomial()==Field_2.optimized_representation()[0].absolute_polynomial()

#Define omega function in terms of chosen delta_generator
def omega(i):
    return ZZ(delta_generator^i);

###Passing to mod p-powers

#Mod p-powers in each entry of a SU-vector
def mod_p_powers(SU,alpha):
    #takes alpha in SU and returns the element x of SU such that x=alpha mod SU^p, and such that x is integral and p-power-free
    M=VectorSpace(GF(5),SU.ngens())
    return SU.exp(M(SU.log(alpha)))

#Check that an element is p-power free
def p_power_check(alpha):
    for i in range(len(list(ideal(alpha).factor()))):
        if list(ideal(alpha).factor())[i][1]%5==0:
            return False
    return True

###Checks that the factorization of the principal ideal generated by the given element does not have any exponents congruent to 5 mod 10
def p_power_10_check(alpha):
    for i in range(len(list(ideal(alpha).factor()))):
        if list(ideal(alpha).factor())[i][1]%10==5:
            return False
    return True

def p_power_adjustment_10(alpha):
    p_power_ideals = [];
    p_power_ideals_exp=[]
    for i in range(len(list(ideal(alpha).factor()))):
        if list(ideal(alpha).factor())[i][1]%10==0:
            p_power_ideals = p_power_ideals + [list(ideal(alpha).factor())[i][0]]
            p_power_ideals_exp = p_power_ideals_exp + [list(ideal(alpha).factor())[i][1]]
    p_power_generator = 1
    for j in range(len(p_power_ideals)):
        p_power_generator = p_power_generator*(p_power_ideals[j]^2).gens_reduced(proof = False)[0]^(-p_power_ideals_exp[j]/2)
    return alpha*p_power_generator

#Function to check omega action
def omega_action_test(alpha,omega_exp_test):
    return SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,delta(alpha)))==SUKz_abs_ell0.log(mod_p_powers(SUKz_abs_ell0,(alpha^omega_exp_test)))

#Kolyvagin derivative operators
def D_0_operator(aut,alpha):
    return prod([(aut^(i))(alpha) for i in range(5)])

def D_1_operator(aut,alpha):
    return prod([((aut^(i))(alpha))^(i) for i in range(5)])

def D_2_operator(aut, alpha):
    return prod([((aut^(i))(alpha))^(binomial(i,2)) for i in range(5)])

###Define functions to check splitting patterns in Kummer extensions

#This function checks if a prime Ell0 splits in a Kummer extension Kz(sqrt[p]{alpha})/Q
def Ell0_splits_check(Ell0,alpha):
    #returns "true" if the given prime Ell0 splits in Kz(sqrt[p]{alpha})/Kz, and "false" otherwise
    F_Ell0 = Ell0.residue_field();
    g=F_Ell0.primitive_element();
    if not 5.divides(F_Ell0(alpha).log(g)): #This is checking if Ell0 is inert in Kz(sqrt[p]{alpha})/Kz
        return False;
    return True;

#This function sorts the primes of characteristic ell_0 into ramified, split, and inert primes in a Kummer extension Kz(sqrt[p]{alpha})/Q
def Ell0_splitting_patterns(alpha):
    ramified_ideals_check = []
    split_ideals_check = []
    inert_ideals_check =[]
    for prime_ideal in Kz_abs.primes_above(11):
        if prime_ideal.divides(ideal(alpha)):
            ramified_ideals_check = ramified_ideals_check+[primes_over_ell0[prime_ideal]]
        elif Ell0_splits_check(prime_ideal,alpha):
            split_ideals_check = split_ideals_check + [primes_over_ell0[prime_ideal]]
        else:
            inert_ideals_check = inert_ideals_check + [primes_over_ell0[prime_ideal]]

    return ramified_ideals_check, split_ideals_check, inert_ideals_check

#Test that the prime of characteristic p are unramified in each extension when p is wildly ramified in K/Q
def wild_p_unramified_check(alpha,beta):
    #returns "true" if prime over 5 in Kz is unramified in Kz(sqrt[p]{beta*alpha})/Kz, and "false" otherwise
    assert p_power_check(beta*alpha)
    if (beta*alpha).mod(P_Kz_25) in elts_all_res_pth_powers_mod:
        return True
    return False

#Test that all primes of characteristic p are unramified in each extension when p is tamely ramified in K/Q
def tame_p_unramified_check(alpha,beta):
    #returns "true" if all primes over 5 in Kz are unramified in Kz(sqrt[p]{beta*alpha})/Kz, and "false" otherwise
    assert p_power_check(beta*alpha)
    for P in Kz_abs.primes_above(5):
        if not ((beta*alpha)^4-1) in P^5: #This is checking if P is ramified in Kz(sqrt[p]{beta*alpha})/Kz
            return False;
    return True;

#Function that tests if a collection of prime ideals are in a Galois orbit under the action of delta
def galois_orbit_check(prime_ideal_indices):
    for ideal_index in prime_ideal_indices:
        if primes_over_ell0[delta(Kz_abs.primes_above(11)[ideal_index])] not in prime_ideal_indices:
                return False
    return True

#Function that takes in an element, coerces the element into an appropriate finite field, computes the log (with 2 fixed as a generator), and then returns the result modulo 5
def log_ell0(FF,n):
    #note that sage automatically views elements of the finite field FF_ell0 as elements in GF(11), so our fixed generator is consistent between the two finite fields
    return GF(5)(FF(n).log(2))

#Function that computes local correction
def local_correction_Ell0(Ell0,alpha,beta):
    prime_ell0_elt_Kz_local=Ell0.0
    assert prime_ell0_elt_Kz_local.ord(Ell0)==1
    assert prime_ell0_elt_Kz_local ==11
    FF_ell0 = Ell0.residue_field()
    return -(zeta_MT*(alpha).ord(Ell0)+ log_ell0(FF_ell0,alpha/(prime_ell0_elt_Kz_local^((alpha).ord(Ell0)))))*(zeta_MT*(beta).ord(Ell0)+ log_ell0(FF_ell0,beta/(prime_ell0_elt_Kz_local^((beta).ord(Ell0)))))^(-1)