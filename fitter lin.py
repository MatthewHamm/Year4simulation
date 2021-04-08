

import numpy as np
from matplotlib import pyplot as plt

from scipy import optimize
import csv

w=(14/2)*np.sqrt(3)
def test_func1(x, a,b,):
    ans=b*(-np.exp((x-w)/a)+np.exp((w-x)/a))/(2*np.sinh(w/a))

    return ans
x=np.linspace(0.001,0.02,20)
y=np.empty(20)

count=0
for i in x:
    print(i)
    ypos=np.array([[0,i]])
    results=np.loadtxt('vlin13 func dt 0.01 i t 200flow %6f.csv'% i,delimiter=',',usecols=(0,1))
    results[:,0]= results[:,0]+np.sqrt(3)/2
    
    ypos=results
    print(ypos)
    params, params_covariance = optimize.curve_fit(test_func1, ypos[:,0], ypos[:,1],p0=[0.1,i])
    y[count]=params[0]
    print(params[0],params[1])
    count+=1
    plt.figure()
    ax=plt.axes()
    ax.set_xlabel('y')
    ax.set_ylabel('|v|')
    xpos=np.linspace(0,w,100)
    ax.scatter(ypos[:,0], ypos[:,1], label='Data')
    ax.set_ylim(0,0.005)
    ax.plot(xpos, test_func1(xpos,params[0],params[1]),
         label='Fitted function',color='red')

    ax.legend(loc='best')
    plt.savefig('vlin small 13 %s'% count, dpi=1000)
    plt.show()

plt.figure()
ax=plt.axes()
ax.set_ylabel('$\delta$')
ax.set_xlabel('v')
ax.scatter(x,y,color='red')

ax.legend(loc='best')
plt.savefig('vlin small delta 13', dpi=1000)
plt.show()