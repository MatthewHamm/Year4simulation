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

equ=filereader('equ.csv',5)





plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('y')
ax.set_xlabel('x')
ax.plot(equ[:,1],equ[:,2])
plt.savefig('equilibrum close',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('y')
ax.set_xlabel('x')
ax.plot(equ[:,3],equ[:,4])
plt.savefig('equilibrum far',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('y')
ax.set_xlabel('Shifted x')
ax.plot(equ[:,0],equ[:,3])
plt.savefig('equilibrum far x',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('y')
ax.set_xlabel('Shifted x')
ax.plot(equ[:,0],equ[:,1])
plt.savefig('equilibrum close x',dpi=1000)
plt.show()

quick=np.loadtxt('vlin new func dt 0.01 i t 200v 1.000000.csv',delimiter=',',usecols=(5,6))
l=len(quick[:-3,0])
t=np.arange(l/100,step=0.01)
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t,quick[:-3,0])
plt.savefig('2chan fast v close',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t,quick[:-3,1])
plt.savefig('2chan fast v far',dpi=1000)
plt.show()
slow=np.loadtxt('vlin func dt 0.01 i t 200v 0.001000.csv',delimiter=',',usecols=(5,6))
l=len(slow[:-3,0])
t=np.arange(l/100,step=0.01)
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t[:-175000],slow[:-175003,0])
plt.savefig('2chan 0001000 small v close',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t[:-175000],slow[:-175003,1])
plt.savefig('2chan 0001000 small v far',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t,slow[:-3,0])
plt.savefig('2chan 0001000 v close',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('vx')
ax.set_xlabel('t')
ax.plot(t,slow[:-3,1])
plt.savefig('2chan 0001000 v far',dpi=1000)
plt.show()
ave=filereader('newavepos11.csv',5)

plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('<x>')
ax.set_xlabel('v')
ax.plot(ave[:,0],ave[:,1])

plt.savefig('close x',dpi=1000)
n=np.argmin(ave[:,1])

print(ave[n,1])
print( ave[n,0])
plt.show()

plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('<y>')
ax.set_xlabel('v')
ax.plot(ave[:,0],ave[:,2])
plt.savefig('close y',dpi=1000)
plt.show()
plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('<x>')
ax.set_xlabel('v')
ax.plot(ave[:,0],ave[:,3])
n=np.argmin(ave[:,3])
x=ave[n,3]
print(ave[n,3])
print( ave[n,0])
plt.savefig('far x',dpi=1000)
plt.show()

plt.figure()
ax = plt.axes()
ax.grid()
ax.set_ylabel('<y>')
ax.set_xlabel('v')
ax.plot(ave[:,0],ave[:,4])
plt.savefig('far y',dpi=1000)
plt.show()
