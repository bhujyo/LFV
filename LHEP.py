############################################################################################
#        This code parses through .lhe files obtained from madgraph
#        The code produces a plot of number of events vs invariant mass of two final state
#         particles
#        Author :  B. Bhattacharya, R. Morgan  (Wayne State University)
#        Date   :  February 17, 2017
############################################################################################
from __future__ import print_function, division
from numpy import *
import matplotlib.pyplot as plt
import glob
from bs4 import BeautifulSoup
import re
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

def invariant_mass(p_list):         # internal function only; input list = [px, ...]
    return sqrt(p_list[3]**2 - p_list[0]**2 - p_list[1]**2 - p_list[2]**2)

def azimuthal_angle(p_list):   # internal function only; input list = [px, ...]
    try:
        phi = arctan(p_list[1]/p_list[0])
    except:
        phi = pi/2
    return phi

def momentum_sum(p1, p2):      # internal function only; input 2 lists = [px, ...], ...
    return [p1[i] + p2[i] for i in range(4)]
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
# function to search events for flavor 1 and 2 leptons
############################################################################################
def lepton12_search(event):
    entry = 0
    muon_check = False
    muon = []
    positron = []
    antimuon = []
    electron = []
    for entry in range(len(event_info(event))):
        if event_info(event)[entry][0] == 13.0:
            muon = event_info(event)[entry]
            entry = len(event_info(event)) + 1
            muon_check = True
        else:
            entry += 1
    for entry in range(len(event_info(event))):
        if event_info(event)[entry][0] == -11.0:
            positron = event_info(event)[entry]
            entry = len(event_info(event)) + 1
        else:
            entry += 1
    for entry in range(len(event_info(event))):
        if event_info(event)[entry][0] == -13.0:
            antimuon = event_info(event)[entry]
            entry = len(event_info(event)) + 1
        else:
            entry += 1
    for entry in range(len(event_info(event))):
        if event_info(event)[entry][0] == 11.0:
            electron = event_info(event)[entry]
            entry = len(event_info(event)) + 1
        else:
            entry += 1

    if muon_check:
        return [muon, positron]
    else:
        return [antimuon, electron]
        
############################################################################################
# funciton to preform equation 3.5 cuts
############################################################################################
def cuts3_5_pass(event):
    failed_particles = 0
    entry = 0
    for entry in range(len(event_info(event))):
        #Transverse momentum cut on muons
        if event_info(event)[entry][0] == 13.0 and p_t(event_info(event)[entry][6:11]) < 20.0:
            failed_particles += 1
            entry += 1
        elif event_info(event)[entry][0] == -13.0 and p_t(event_info(event)[entry][6:11]) < 20.0:
            failed_particles += 1
            entry += 1
        #Transverse momentum cut on electrons
        elif event_info(event)[entry][0] == 11.0 and p_t(event_info(event)[entry][6:11]) < 20.0:
            failed_particles += 1
            entry += 1
        elif event_info(event)[entry][0] == -11.0 and p_t(event_info(event)[entry][6:11]) < 20.0:
            failed_particles += 1
            entry += 1
        #Rapidity cut on muons
        if event_info(event)[entry][0] == 13.0 and rapidity(event_info(event)[entry][6:11]) > 2.5:
            failed_particles += 1
            entry += 1
        elif event_info(event)[entry][0] == -13.0 and rapidity(event_info(event)[entry][6:11]) > 2.5:
            failed_particles += 1
            entry += 1
        #Rapidity cut on electrons
        elif event_info(event)[entry][0] == 11.0 and rapidity(event_info(event)[entry][6:11]) > 2.5:
            failed_particles += 1
            entry += 1
        elif event_info(event)[entry][0] == -11.0 and rapidity(event_info(event)[entry][6:11]) > 2.5:
            failed_particles += 1
            entry += 1
       
    if failed_particles == 0:
        return True
    else:
        return False

############################################################################################
# function to preform equation 3.7 cuts
############################################################################################
def cuts3_7_pass(event):
    failed_particles = 0
    if abs(azimuthal_angle(lepton12_search(event)[0][6:11]) - azimuthal_angle(lepton12_search(event)[1][6:11])) < 2.75:
        failed_particles = 1
    elif p_t(lepton12_search(event)[0][6:11]) < p_t(lepton12_search(event)[1][6:11]):
        failed_particles = 1
    #elif abs(azimuthal_angle(lepton12_search(event)[0][6:11]) - azimuthal_angle(lepton12_search(event)[1][6:11])) > 0.6:
    
    if failed_particles == 0:
        return True
    else:
        return False

############################################################################################
# function to preform equation 3.8 cut
############################################################################################
def cuts3_8_pass(event):
    failed_particles = 0
    if invariant_mass(event, 'mu+', 'ta-') < 250.0 and invariant_mass(event, 'mu-', 'ta+') < 250.0:
        failed_particles = 1

    if failed_particles == 0:
        return True
    else:
        return False

############################################################################################
# function to determine if an event has passed all cuts
############################################################################################
def cut_pass(event): 
    if cuts3_5_pass(event) and cuts3_7_pass(event) and cuts3_8_pass(event):
        return True
    else:
        return False


############################################################################################
############################################################################################
# Begin analysis of events
############################################################################################
############################################################################################

print('Retrieving events...')

############################################################################################
# Open the .lhe files in python
############################################################################################
file_handles = glob.glob("*.lhe")
if len(file_handles) == 0:
    print("I need at least one .lhe file to parse")
    exit()
############################################################################################
# Open the .lhe file in python, construct an xml tree, and isolate the event generation info
############################################################################################
invariant_mass_list = []
n_of_events = []
cs = []
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
#Call the cut function for each event
############################################################################################
cut_events = 0
kept_events = []

for event in all_events:
    if cut_pass(event):
        kept_events.append(event)
    else:
        cut_events += 1

############################################################################################
# Output : generation information
############################################################################################
print('------------------------------------------------------')
print('Total number of events analyzed:    ', sum(n_of_events))
print('Cross section before cuts:          ', sum([n_of_events[i] * cs[i] for i in range(len(n_of_events))])
     /sum(n_of_events), 'pb')
print('Number of kept events:              ', len(kept_events))
print('Number of cut events:               ', cut_events)
print('Cross section after cuts:           ', (len(kept_events)) / 10000.0 * sum([n_of_events[i] * cs[i] for i in range(len(n_of_events))])
     /sum(n_of_events), 'pb')
print('------------------------------------------------------')
############################################################################################
# Output : plots
############################################################################################
#plt.hist(invariant_mass_list, bins=300)
#plt.title('histogram for ta+ ta- events distribution')
#plt.show()
