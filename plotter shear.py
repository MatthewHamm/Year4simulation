# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 21:07:36 2020

@author: Matthew Hamm
"""

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import csv
from scipy import optimize

def filereader(filename,n):
    count=0
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
    
        results=[np.zeros(n)]

        for row in csv_reader:
            count+=1

            row0=np.array(row).astype(np.float)
               
            results=np.append(results,[row0],axis=0)
            

            

        
        return results[1 :,:]


quick=np.loadtxt('vlin func i int dt 0.01  t 100v 1.000000.csv', delimiter=',',usecols=(3))
l=len(quick[:-3])
t=np.arange(l/100,step=0.01)
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t,quick[:-3])
plt.savefig('v 1',dpi=1000)
plt.show()
slow=np.loadtxt('vlin func i int dt 0.01  t 100v 0.001000.csv', delimiter=',',usecols=(3))
l=len(slow[:-3])
t=np.arange(l/100,step=0.01)
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t[:-175000],slow[:-175003])
plt.savefig('v short 001',dpi=1000)
plt.show()
ave=filereader('ave_posfinal.csv',3)
n=np.argmin(ave[:,1])
print(np. min(ave[:,1]))
print(ave[n,0 ])
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('<x>')
ax.set_xlabel('v')
ax.plot(ave[:,0],ave[:,1])

plt.savefig('average x',dpi=1000)
plt.show()

plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('<y>')
ax.set_xlabel('v')
ax.plot(ave[:,0],ave[:,2])
plt.savefig('average y',dpi=1000)
plt.show()
equ=filereader('equilibrum.csv',3)

print(equ[:,0])



plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('y')
ax.set_xlabel('x')
ax.plot(equ[:,1],equ[:,2])
plt.savefig('equilibrum',dpi=1000)
plt.show()
