# This repository is for parsing files in .lhe format
.lhe stands for "Les-Houches Events" a formating type used to output
event data from a simulator such as MadGraph. This repo hasthe following
python codes :

LHEP.py    : in this code are functions defined to analyze any .lhe file
Analyze.py :this code calls LHEP.py as a module and shows a scenario for 
its implementation

## How to use the code

Leave LHEP.py unchanged for your specific scenario of .lhe file(s)
Implement your changes in Analyze.py suitable for your own needs

### Predefined functions 

The following functions have been defined in LHEP.py:

1. particle_info(event, particle) : produces a list [px, py, pz, E, m]
2. pt(event, particle) : generates the pt. If particle is not specified 
then it generates the missing pt for the event
3. invariant_mass(event, particle1, particle2) : produces the invariant 
mass for the 2 particles

To be continued ...
