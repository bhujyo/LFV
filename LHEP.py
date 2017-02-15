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
# input numbers
############################################################################################
par_in = {'mt':[173.21,0.51,0.21]}
############################################################################################
# list analysis function
def centerr(in_list, sigma=1):
    cent = in_list[0]
    err = sigma * sqrt(sum([x**2 for x in in_list[1:]]))
    return [cent - err, cent + err]
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
# functions to extract info from an event
############################################################################################
def event_info(event):                    # list with information about particles in event
    return [[float(y) for y in x.strip().split()]
        for x in event.text.strip().split('\n')[1:-6]]
def particle_info(event, particle):        # generates a list : [[px, py, pz, E, m], ...]
    entries = event_info(event)            # ordered by largest p_t
    if abs(pdg_id_dict[particle]) in [12.0, 14.0, 16.0]:
        return [0, 0, 0, 0, 0]
    info = []
    for line in entries:
        if line[0] == pdg_id_dict[particle]:
            info.append(line[6:11])
        elif abs(line[0]) <= 6.0 and abs(line[0]) == pdg_id_dict[particle]:
            info.append(line[6:11])
    return sorted(info, key = lambda x : p_t(x), reverse = True)
############################################################################################
# functions to extract kinematic information, pt, rapidity, azimuthal angle, etc.
############################################################################################
def p_t(p_list):      # internal function only; input list = [px, py, pz, E, m]
    return sqrt(p_list[0]**2 + p_list[1]**2)

def rapidity(p_list):  # internal function only; input list = [px, py, pz, E, m]
    p = sqrt(sum([x**2 for x in p_list]))
    pz = p_list[2]
    try:
        eta = 0.5 * log((p + pz)/(p - pz))
    except:
        eta = 4.0
    return eta

def momentum_sum(p1, p2):      # internal function only; input 2 lists = [px, ...], ...
    return [p1[i] + p2[i] for i in range(4)]

def invariant_mass(p):         # internal function only; input list = [px, ...]
    return sqrt(p[3]**2 - p[0]**2 - p[1]**2 - p[2]**2)

def azimuthal_angle(particle_info):   # internal function only; input list = [px, ...]
    try:
        phi = arctan(particle_info[1]/particle_info[0])
    except:
        phi = pi/2
    return phi
############################################################################################
# function to extract the missing transverse momentum from an event
############################################################################################
def miss_pt(event):
    entries = event_info(event)
    pxy = [0, 0]
    for item in entries:
        if item[1] == 1.0 and rapidity(item[6:11]) <= 3.5:
            pxy[0] = pxy[0] + item[6]
            pxy[1] = pxy[1] + item[7]
    return -sqrt(pxy[0]**2 + pxy[1]**2)
############################################################################################
