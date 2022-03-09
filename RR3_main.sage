#To use this program, load this page and call the function RR3_comptuations().
#Input: p_value, ell0_value, ell1_value, do_a1_p_check=True/False, mode = "Full"/"Ell1"
#Output: Conclusion of Theorem 1.2.1 in RR3 Paper II, including all required computations

def RR3_computations(p_value,ell0_value,ell1_value,do_a1_p_check=True, mode="Full"):

	#Define input value p_value, ell0_value, ell1_value as global objects
	globals()["p_value"]=p_value;
	globals()["ell0_value"]=ell0_value;
	globals()["ell1_value"]=ell1_value;

	#Check that valid mode is given in input
	assert mode in ["Full", "Ell1"], "Invalid mode provided"

	#Check that Assumption 2.1.1 in RR3 Paper II is met
	assert p_value.is_prime(), "p_value is not prime"
	assert p_value >= 5, "p_value is not greater than or equal to 5"
	assert ell0_value.is_prime(), "ell0_value is not prime"
	assert ell1_value.is_prime(), "ell1_value is not prime"
	assert GF(p_value)(ell0_value) == 1, "Assumption (1) is not met"
	assert GF(p_value)(ell1_value)%p_value!=0 and (GF(p_value)(ell1_value)%p_value)^2!=1, "Assumption (2) is not met"
	assert GF(ell0_value)(ell1_value).log(GF(ell0_value).primitive_element())%p_value == 0, "Assumption (3) is not met"   
	assert GF(ell0_value)(prod([i^i for i in [1..(ell0_value-1)/2]])).log(GF(ell0_value).primitive_element())%5!=0, "Assumption (4) is not met" 

	print(f'We take p={p_value}, ell0={ell0_value}, and ell1={ell1_value}:\n');

	if mode == "Full":
		print(f'Computing initial set up using {mode} mode:');
		load("initial_setup.sage");
		load("functions.sage");
		load("Qz_extensions.sage");
		
	elif mode == "Ell1":
		print(f'Computing initial set up using {mode} mode:');
		load("initial_setup_ell1.sage");
		load("functions.sage");

	load("Kolyvagin_setup.sage");
	load("a1_side.sage");
	if do_a1_p_check:
		load("a1_p_splitting.sage");
	load("b2_prelim.sage");
	load("b2_side.sage");