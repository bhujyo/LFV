# This repository is for parsing files in .lhe format
.lhe stands for "Les-Houches Events" a formating type used to output
event data from a simulator such as MadGraph. This repo hasthe following
python codes :

LHEP.py    : in this code are functions defined to analyze any .lhe file \\
Analyze.py : this code calls LHEP.py as a module and shows a scenario for
its implementation

## How to use the code

Leave LHEP.py unchanged for your specific scenario of .lhe file(s)
Implement your changes in Analyze.py suitable for your own needs

### Predefined functions 

In the current version, the following functions have been defined in
LHEP.py (Analyze.py uses all these functions):

1. event_info(event) : reformats the entire data within an event into a giant list of lists
2. particle_info(event, particle) : produces a list in the format [[px, py, pz, E, m], ...]
The output lists within the big list is arranged in a descending order by transverse momentum (p_t)
3. p_t(list) : generates the the transverse momentum pt from the list [px, py, pz, E, m].
4. rapidity, invariant_mass, and azimuthal_angle are similar functions as above
5. momentum_sum(p1, p2) : generates the four momentum p12 = p1 + p2 and returns it in the form
[p12x, p12y, p12z, E12]
6. miss_pt(event) : generates the missing transverse momentum in an event (any particle that is a
neutrino or has rapidity larger than 3.5 contributes to this missing pt.)

To be continued ...
