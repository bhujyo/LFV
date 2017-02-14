############################################################################################
#        This code parses through .lhe files obtained from madgraph
#        The code produces a plot of number of events vs invariant mass of two final state
#         particles
#        Author :  B. Bhattacharya  (Wayne State University)
#        Date   :  February 2, 2017
############################################################################################
from __future__ import print_function, division
from numpy import *
############################################################################################
# pdg particle dictionary
############################################################################################
pdg_id_dict = {'d': 1.0, 'd~':-1.0,'u': 2.0, 'u~':-2.0, 's': 3.0, 's~':-3.0, 
'c': 4.0, 'c~':-4.0, 'b': 5.0, 'b~':-5.0,'t': 6.0, 't~':-6.0, 'e-': 11.0, 'e+':-11.0, 
've': 12.0, 've~':-12.0, 'mu-': 13.0,'mu+':-13.0, 'vm': 14.0, 'vm~':-14.0,'ta-': 15.0, 
'ta+':-15.0, 'vt': 16.0, 'vt~':-16.0, 'g':21.0, 'a':22.0, 'z':23.0, 'w+':24.0, 'w-':-24.0,
'h0':25.0}
############################################################################################
# function to appropriately format a string
############################################################################################
LHEf_cols = ['pdg_id', 'status', 'parent', 'x', 'c1', 'c2', 'px', 'py', 'pz', 'E', 'm', 'x',
     'helicity']
############################################################################################
# function to appropriately format a string
############################################################################################
def str_format(string):
    string = string.strip()
    try:
        return float(string)
    except:
        return string
############################################################################################
# functions to extract a particle's kinematic information, pt, rapidity, azimuthal angle
############################################################################################
def particle_info(event, particle):              # generates a list : [px, py, pz, E, m]
    entries = [[float(y) for y in x.strip().split()]
        for x in event.text.strip().split('\n')[1:-6]]
    info = []
    for line in entries:
        if line[0] == pdg_id_dict[particle]:
            info.append(line[6:11])
    return info

def p_t(particle_info):      # internal function only; input list = [px, py, pz, E, m]
    return sqrt(particle_info[0]**2 + particle_info[1]**2)

def rapidity(particle_info)  # internal function only; input list = [px, py, pz, E, m]
    p = sqrt(sum([x**2 for x in particle_info]))
    pz = particle_info[2]
    return 0.5 * log((p + pz)/(p - pz))

def invariant_mass(p1_info, p2_info): # internal function only; input list = [px, ...]
    p12_info = [p1_info[i] + p2_info[i] for i in range(4)]
    return sqrt(p12_info[3]**2 - p12_info[0]**2 - p12_info[1]**2 - p12_info[2]**2)
############################################################################################
# function to extract the azimuthal angle from an event
############################################################################################
#def azimuthal_angle(event, particle):
#    angle = particle_info(event, particle)
#    if len(angle) == 0:
#        return 0
#    return arctan(angle[1] / angle[0])
############################################################################################
# function to extract the rapidity from an event
############################################################################################
#def rapidity(event, particle):
#    eta = particle_info(event, particle)
#    if len(eta) == 0:
#        return 0
#    return ln( 2 eta[2] / eta[0])
############################################################################################
# function(s) to extract transverse momentum for a particle and the missing transverse 
# momentum from an event
############################################################################################
def miss_pt(event):
    entries = [[float(y) for y in x.strip().split()]
        for x in event.text.strip().split('\n')[1:-6]]
    pxy = [0, 0]
    for item in entries:
        if item[1] == 1.0:
            pxy[0] = pxy[0] + item[6]
            pxy[1] = pxy[1] + item[7]
    return -sqrt(pxy[0]**2 + pxy[1]**2)
def pt(event, particle='miss'):
    if particle == 'miss':
        return miss_pt(event)  
    try:
        p = particle_info(event, particle)[0]
    except:
        p = [0, 0]
    return sqrt(p[0]**2 + p[1]**2)
############################################################################################
