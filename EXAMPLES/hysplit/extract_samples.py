import numpy as np
import subprocess
import datetime
import os,sys
import re
import shutil
from extract_wind import write_atm
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 

from input_file import *

INDEX, LAT, LON = np.loadtxt('con2stn.inp',skiprows=0, unpack=True)

nsampl = len(INDEX)

block_length = nsampl * npart

con2stn = np.loadtxt('con2stn.txt',skiprows=1)

# print 'size',con2stn.shape

JDAY = con2stn[:,0]

Value = con2stn[:,-1]

# print 'JDAY',JDAY
# print 'Value',Value

con2std_len = len(JDAY)

nblocks = con2std_len/block_length

loading = np.zeros((nblocks,nsampl,npart))

for i in range(nblocks):
    for j in range(nsampl):
        loading[i,j,:] = Value[i*block_length+(j+1)+nsampl*np.arange(npart)-1]


fig = plt.figure()
ax = plt.subplot(111)

width = 0.9/nsampl

legend_strings = []
color=iter(cm.rainbow(np.linspace(0,1,nsampl)))

for j in range(nsampl):
    c=next(color)
    ax.bar(np.array(diam_phi)+(j-0.5*nsampl)*width, loading[nblocks-1,j,:]/np.sum(loading[nblocks-1,j,:])*100,width,color=c)
    stringj = "Loc %s (%s,%s), Loading=%.2e [g/m2]" % (str(j+1), str(LAT[j]), str(LON[j]),np.sum(loading[nblocks-1,j,:])*100 )
    legend_strings.append(stringj)

ax.legend(legend_strings)

ax.set_ylabel('Loading [mass wt%]')
ax.set_title('Loading at sampling locs')

ax.set_xlabel('Diameter [phi]')
ax.set_xticks(np.array(diam_phi))

# ax.legend((rects1[0], rects2[0]), ('Men', 'Women'))

fig.savefig('gsd.pdf')


