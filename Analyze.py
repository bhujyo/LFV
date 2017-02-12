############################################################################################
#        This code parses through .lhe files obtained from madgraph
#        The code produces a plot of number of events vs invariant mass of two final state
#         particles
#        Author :  B. Bhattacharya  (Wayne State University)
#        Date   :  February 2, 2017
############################################################################################
from LHEP import *
import glob
from bs4 import BeautifulSoup
from numpy import *
import re
import matplotlib.pyplot as plt
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
plot_list = []
for handle in file_handles:
    data = BeautifulSoup(open(handle, 'r').read(), "lxml")
    gen_info = dict([[str_format(y) for y in x.split(':')]
        for x in data.mggenerationinfo.string.strip().split('\n')])
    n_of_events.append(gen_info['#  Number of Events'])
    cs.append(gen_info['#  Integrated weight (pb)'])
    all_events = data.find_all('event')
    for event in all_events:
        plot_list.append(particle_info(event, 'b')[0][3]/particle_info(event, 't')[0][3])      
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
plt.hist(plot_list, bins=300)
plt.title('histogram for t tbar events distribution')
plt.xlim([0.0, 1.1])
plt.show()
