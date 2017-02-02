############################################################################################
#        This code parses through .lhe files obtained from madgraph
#        Author :  B. Bhattacharya  (Wayne State University)
#        Date   :  February 2, 2017
############################################################################################
from __future__ import print_function
from bs4 import BeautifulSoup
import numpy
import re
############################################################################################
# pdg particle dictionary
############################################################################################
pdg_sm_particle_dict = {'d': 1.0, 'dbar':-1.0,'u': 2.0, 'ubar':-2.0, 's': 3.0, 'sbar':-3.0, 
'c': 4.0, 'cbar':-4.0, 'b': 5.0, 'bbar':-5.0,'t': 6.0, 'tbar':-6.0, 'e-': 11.0, 'e+':-11.0, 
've': 12.0, 've~':-12.0, 'mu-': 13.0,'mu+':-13.0, 'vm': 14.0, 'vm~':-14.0,'ta-': 15.0, 
'ta+':-15.0, 'vt': 16.0, 'vt~':-16.0, 'g':21.0, 'z':21.0, 'w+':24.0, 'w-':-24.0}
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
# Open the .lhe file in python, construct an xml tree, and isolate the event generation info
############################################################################################
data = BeautifulSoup(open("unweighted_events.lhe", 'r').read(), "lxml")
gen_info = dict([[str_format(y) for y in x.split(':')] for x in 
           data.mggenerationinfo.string.strip().split('\n')])
n_of_events = gen_info['#  Number of Events']
total_cs = gen_info['#  Integrated weight (pb)']
############################################################################################
# Output : generation information
############################################################################################
print('Total number of events analyzed:', n_of_events)
print('Total cross section:', total_cs, 'pb')
print('Reading events ...')
############################################################################################
events = data.find_all('event')
############################################################################################
# function to extract pt of a particle
############################################################################################
def pt(particle,event,order=1):
    events = [[float(y) for y in x.strip().split()] for x in event.text.strip().split('\n')[1:]]
    if particle == 'miss':
        return sum([x**2 for x in miss(event)])
    pt_dict = {}
    for item in events:
        pt_dict[item[0]] = pt_dict.get(item[0], [])
        pt_dict[item[0]].append(numpy.sqrt(item[6]**2 + item[7]**2))
        pt_dict[item[0]].sort(reverse=True)
    try:
        if len(pt_dict[pdg_sm_particle_dict[particle]]) > order:
            return -1.0
        else:
            return pt_dict[pdg_sm_particle_dict[particle]][order-1]
    except:
        return -1.0
############################################################################################
for event in events:
    print(event.text.strip().split('\n')[1:])
