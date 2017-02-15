############################################################################################
#        This code parses through .lhe files obtained from madgraph
#        The code produces a plot of number of events vs invariant mass of two final state
#         particles
#        Author :  B. Bhattacharya  (Wayne State University)
#        Date   :  February 2, 2017
############################################################################################
from __future__ import print_function, division
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
# Open the .lhe file in python, construct an xml tree, and isolate the event generation info
############################################################################################
n_of_events = []
cs = []
plot_list = []
count1 = 0
count2 = 0
for handle in file_handles:
    data = BeautifulSoup(open(handle, 'r').read(), "lxml")
    gen_info = dict([[str_format(y) for y in x.split(':')]
        for x in data.mggenerationinfo.string.strip().split('\n')])
    n_of_events.append(gen_info['#  Number of Events'])
    cs.append(gen_info['#  Integrated weight (pb)'])
    all_events = data.find_all('event')
    flags = []
    for event in all_events:
        binfo = particle_info(event, 'b')
        wpinfo = particle_info(event, 'w+')[0]
        flag = 0
        for info in binfo:
            if centerr(par_in['mt'],5)[0] <= invariant_mass(momentum_sum(info,
                               wpinfo)) <= centerr(par_in['mt'],5)[1]:
                EboEt = info[3] / momentum_sum(info, wpinfo)[3]
                flag += 1
        if flag == 1: plot_list.append(EboEt); count1 += 1
        if flag > 1: count2 += 1
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
print('Total number of events used = ', count1, count2)
############################################################################################
plt.hist(plot_list, bins=40)
plt.title('histogram for missing p_t distribution in ttbar events')
plt.show()
