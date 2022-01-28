#  @author: Yifeng Peng PID: 730366058

import numpy as np
import matplotlib.pyplot as plt


def Dfunc(t,Y,params):
    [I1,I3,M,g,l] = params
    theta = Y[0]
    p_theta = Y[1]
    phi = Y[2]
    p_phi = Y[3]
    psi = Y[4]
    p_psi = Y[5]
    
    def kick(t):
        tperiod = 0.5
        if t/tperiod > 1 and t % tperiod < 0.01:
            kick = 1
        else:
            kick = 1
        return kick
    
    ret = np.zeros(6)
    
    ret[0] = p_theta/I1
    ret[1] = (p_phi-p_psi*np.cos(theta))**2*np.cos(theta)/(I1*np.sin(theta)**3) - p_psi*(p_phi-p_psi*np.cos(theta))/(I1*np.sin(theta)) + M*kick(t)*g*l*np.sin(theta)
    ret[2] = (p_phi-p_psi*np.cos(theta))/(I1*np.sin(theta)**2)
    ret[3] = 0
    ret[4] = p_psi/I3 - (p_phi-p_phi*np.cos(theta))*np.cos(theta)/(I1*np.sin(theta)**2)
    ret[5] = 0
    return ret

def RK4(f,Y0,dt,T,params):
    t = np.arange(0,T,dt)
    iteration = len(t)
    Y = np.zeros([len(Y0),iteration])
    Y[:,0] = Y0[:,0]
    for i in range(iteration-1):
        k1 = f(t[i],Y[:,i],params)
        k2 = f(t[i]+dt/2,Y[:,i]+dt/2*k1,params)
        k3 = f(t[i]+dt/2,Y[:,i]+dt/2*k2,params)
        k4 = f(t[i]+dt,Y[:,i]+dt*k3,params)
        Y[:,i+1] = Y[:,i]+dt/6*(k1+2*k2+2*k3+k4)
    return t,Y