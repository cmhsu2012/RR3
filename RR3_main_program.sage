###Choose value of ell1
ell1_value = 67
assert not 5.divides(ell1_value-1);
assert not 5.divides(ell1_value+1);
assert GF(11)(ell1_value^2)==1;

load("RR3_initial_setup.sage")
load("RR3_functions.sage")
load("RR3_Kolyvagin_setup.sage")
load("RR3_a1_side.sage")
load("RR3_b2_side.sage")