Version 1.0

2022-03-08

Authors: Catherine Hsu, Preston Wake, Carl Wang-Erickson

This code is written for Sage, Version 9.2.

*****************************************************************************
       Copyright (C) 2022 Hsu--Wake--Wang-Erickson 

  Distributed under the terms of the GNU General Public License (GPL)

    This code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

  The full text of the GPL is available at:

                  http://www.gnu.org/licenses/
*****************************************************************************

To use this program, load RR3_main.sage. The main program function RR3_main takes the following inputs:

- values for p, ell0, and ell1 satisfying Assumption 2.1.1
- do_a1_p_check = True or False, True if you want the program to check if a1 is split at p, False otherwise
- 

Full mode --

1. initial_setup.sage: This implements the algorithms in Section 3.1.
2. functions.sage: This implements the functions used throughout Sections 3.2 to 3.5.
3. Qz_extensions.sage:
4. Kolyvagin_setup.sage: This implements the algorithms in Section 3.2.
5. a1_side.sage: This implements the algorithms in Section 3.3.
6. a1_p_splitting.sage
7. b2_prelim.sage: This implements the algorithms in Section 3.4 and 3.5.
8. b2_side.sage

Ell1 mode 
- initial_setup_ell1.sage:
- Qz_extensions.sage

