

import numpy as np
from matplotlib import pyplot as plt

from scipy import optimize
import csv

w=(14/2)*np.sqrt(3)
def test_func1(x, a,b,c):
    R=np.sqrt((1+b)/(2*a))
    I=np.sqrt(b*(1+b)/(2*a))
    ans=(-b)*np.sqrt(np.square(np.sinh(R*(x-w))*np.cos(I*(x-w))*np.sinh(R*(w))*np.cos(I*(w))+np.cosh(R*(x-w))*np.sin(I*(x-w))*np.cosh(R*(w))*np.sin(I*(w)))+np.square(-np.sinh(R*(x-w))*np.cos(I*(x-w))*np.cosh(R*(w))*np.sin(I*(w))+np.cosh(R*(x-w))*np.sin(I*(x-w))*np.sinh(R*(w))*np.cos(I*(w))))
    ans=ans/(np.sinh(R*(-w))*np.cos(I*(-w))*np.sinh(R*(w))*np.cos(I*(w))+np.cosh(R*(-w))*np.sin(I*(-w))*np.cosh(R*(w))*np.sin(I*(w)))
    return ans
x=np.linspace(0.1,1,10)
y=np.empty(10)

count=0
for i in x:
    print(i)
    ypos=np.array([[0,i]])
    results=np.loadtxt('vosc13 func dt 0.01 i t 200flow %6f.csv'% i,delimiter=',',usecols=(0,1))
    results[:,0]= results[:,0]+np.sqrt(3)/2
    ypos=np.append(ypos,results,0)

    print(ypos)
    params, params_covariance = optimize.curve_fit(test_func1, ypos[:,0], ypos[:,1],p0=[0.1,i,0.5])
    y[count]=params[0]
    print(params[0],params[1])
    count+=1
    plt.figure()
    ax=plt.axes()
    ax.set_xlabel('y')
    ax.set_ylabel('|v|')
    ax.set_ylim(0,0.005)
    ypos=results
    ax.scatter(ypos[:,0], ypos[:,1], label='Data')
    xpos=np.linspace(0,w,100)
    ax.plot(xpos, test_func1(xpos,params[0],params[1],params[2]),
         label='Fitted function',color='red')

    ax.legend(loc='best')
    plt.savefig('vosc 13 %s'% count, dpi=1000)
    plt.show()
print(y[9])
plt.figure()
ax=plt.axes()
ax.set_ylabel('$\delta$')
ax.set_xlabel('v')
ax.scatter(x, y,color='red')


plt.savefig('vosc delta 13', dpi=1000)
plt.show()