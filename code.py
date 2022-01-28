#  @author: Yifeng Peng PID: 730366058

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
from util import Dfunc, RK4

I1 = 3/4; I3 = 3/10;
M = 1; g = 9.8; l = 1;
params = [I1,I3,M,g,l]

theta0 = np.pi/6; phi0 = -np.pi/2; p_phi0 = 10; psi0 = 0;

# resolution is used to speed up the graphing process, it won't affect the
# quality of the solution, which is related to dt. However, it will affect
# the quality of the time dependent drawing

# #for steady precession
# resolution = 100
# p_theta0 = 0
# dt = 0.001

#for backward swirl precession
resolution = 100
p_theta0 = 1
dt = 0.001

# #for critical nutation
# resolution = 100
# p_theta0 = 0.3
# dt = 0.001

p_psi0 = 0.5*(2*p_phi0*np.cos(theta0) + p_phi0*np.sin(theta0)*np.tan(theta0) + np.sqrt(-4*g*I1*l*M*np.tan(theta0)**-1*np.sin(theta0)**-1+p_phi0**2*np.sin(theta0)**-2)*np.sin(theta0)**2*np.tan(theta0))

Y0 = np.array([[theta0],[p_theta0],[phi0],[p_phi0],[psi0],[p_psi0]])

T = 7

t,Y = RK4(Dfunc,Y0,dt,T,params)
theta = Y[0]
phi = Y[2]
psi = Y[4]
for i in range(len(phi)):
    if abs(phi[i] - 1.5*np.pi) < 1e-2:
        print("the period of precession is about {:.3f}s".format(t[i]))
        break
p_theta = Y[1]
p_psi = Y[5]
p_phi = Y[3]
Etot = p_theta**2/(2*I1)+(p_psi-p_phi*np.cos(theta))**2/(2*I1*np.sin(theta)**2)+p_phi**2/(2*I3)+M*g*l*np.cos(theta)

max_theta = max(theta)/np.pi * 180
min_theta = min(theta)/np.pi * 180
print("the maximum theta is {:.2f} degree and the minimum theta is {:.2f} degree".format(max_theta,min_theta))

x = l*np.sin(theta)*np.cos(phi)
y = l*np.sin(theta)*np.sin(phi)
z = l*np.cos(theta)

dispRange = range(0,len(x),resolution)
trace = np.zeros([3,len(dispRange)])
index = 1
trace[0,0] = x[0]
trace[1,0] = y[0]
trace[2,0] = z[0]
for i in dispRange:
    fig = plt.figure(1)
    if i+resolution < len(x):
        trace[0,index] = x[i+resolution]
        trace[1,index] = y[i+resolution]
        trace[2,index] = z[i+resolution]
    
    #plot axis
    ax = plt.axes(projection='3d')
    ax.plot3D([0,0,0,0,0,0,l,-l],[0,0,0,l,-l,0,0,0],[l,-l,0,0,0,0,0,0],"r-",linewidth=0.5)
    #plot the top and it's coordinate visualizer
    ax.plot3D([x[i],0],[y[i],0],[z[i],0],"b-")
    ax.plot3D([x[i],x[i],x[i],0,0],[0,y[i],y[i],y[i],0],[0,0,z[i],z[i],z[i]],"b--",linewidth=1)
    #plot the trajectory
    ax.plot3D(trace[0,0:index],trace[1,0:index],trace[2,0:index],"k-",linewidth=0.7)
    ax.set_title("t = {:.3f}".format(t[i]))
    index += 1
    
    ax.set_zlim([-l,l])
    ax.set_xlim([-l,l])
    ax.set_ylim([-l,l])
    
    plt.draw()
    plt.pause(0.05)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D([0,0,0,0,0,0,l,-l],[0,0,0,l,-l,0,0,0],[l,-l,0,0,0,0,0,0],"r-",linewidth=0.5)
ax.plot3D(x,y,z,"k-",linewidth=0.7)
ax.plot3D([x[-1],0],[y[-1],0],[z[-1],0],"b-",linewidth=1)
ax.set_title("trajectory")
ax.set_zlim([-l,l])
ax.set_xlim([-l,l])
ax.set_ylim([-l,l])

plt.figure()
plt.plot(t,Etot)
plt.title("energy vs time")
plt.ylim([0,max(Etot)*1.5])

plt.figure()
plt.plot(phi[0:round(len(phi)/4)],theta[0:round(len(phi)/4)])
# plt.plot(phi,theta)
plt.title("theta vs phi")

plt.figure()
plt.plot(phi[0:round(len(phi)/4)],theta[0:round(len(phi)/4)])
# plt.plot(t,theta)
plt.title("theta vs time")