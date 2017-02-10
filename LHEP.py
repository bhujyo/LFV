############################################################################################
#        This code parses through .lhe files obtained from madgraph
#        The code produces a plot of number of events vs invariant mass of two final state
#         particles
#        Author :  B. Bhattacharya  (Wayne State University)
#        Date   :  February 2, 2017
############################################################################################
from __future__ import print_function, division
import glob
from bs4 import BeautifulSoup
from numpy import *
import re
import matplotlib.pyplot as plt
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
# function to extract a particle's kinematic information
############################################################################################
def particle_info(event, particle):
    entries = [[float(y) for y in x.strip().split()]
        for x in event.text.strip().split('\n')[1:-6]]
    info = []
    for line in entries:
        if line[0] == pdg_id_dict[particle]:
            info.append(line[6:11]) 
    return info
############################################################################################
# function to extract the invariant mass of 2 particles
############################################################################################
def invariant_mass(event, particle1, particle2):
    p1_info = particle_info(event, particle1)
    p2_info = particle_info(event, particle2)
    if len(p1_info) == 0 or len(p2_info) == 0:
        return 0
    p12_info = [p1_info[0][i] + p2_info[0][i] for i in range(4)]
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
def pt(event, particle):
    if particle == 'miss':
        return miss_pt(event)  
    try:
        p = particle_info(event, particle)[0]
    except:
        p = [0, 0]
    return sqrt(p[0]**2 + p[1]**2)
############################################################################################
# Open the .lhe files in python
############################################################################################
file_handles = glob.glob("*.lhe")
if len(file_handles) == 0:
    print("I need at least one .lhe file to parse")
    exit()
############################################################################################
# Cuts on events
############################################################################################
#
#    
#    
#    
#
#
#
#
#
#
#
############################################################################################
# Open the .lhe file in python, construct an xml tree, and isolate the event generation info
############################################################################################
n_of_events = []
cs = []
invariant_mass_list = []
for handle in file_handles:
    data = BeautifulSoup(open(handle, 'r').read(), "lxml")
    gen_info = dict([[str_format(y) for y in x.split(':')]
        for x in data.mggenerationinfo.string.strip().split('\n')])
    n_of_events.append(gen_info['#  Number of Events'])
    cs.append(gen_info['#  Integrated weight (pb)'])
    all_events = data.find_all('event')
    for event in all_events:
        invariant_mass_list.append(invariant_mass(event, 'ta+', 'ta-'))
############################################################################################
# We will implement checks here : usually commented out
############################################################################################
# print(pt(all_events[-1], 'ta+'))
############################################################################################
# Output : generation information
############################################################################################
print('Total number of events analyzed:', sum(n_of_events))
print('Total cross section:', sum([n_of_events[i] * cs[i] for i in range(len(n_of_events))])
     /sum(n_of_events), 'pb')
print('Retrieving events ...')
############################################################################################
plt.hist(invariant_mass_list, bins=300)
plt.title('histogram for ta+ ta- events distribution')
plt.show()
