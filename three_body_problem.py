# Python simulation for the 3 body problem
# Different initial conditions give rise to different configurations.
# Author: Hemaditya Malla


import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass #This is a handier version of just using the regular classes
import copy

def computeForce(masses,positions):
    #Computes the gravitational force of attraction by k-th,j-th body on a i-th body at ti
    # and returns a 2 element array [Fx, Fy] (#TODO: can this be generalized?)
    ri, rj, rk = positions
    mj, mk = masses
    dr = [ri-rj, ri-rk]
    Fi = sum([G*mass*(r/np.linalg.norm(r)**3) for mass,r in zip([mj,mk], dr)])
    return Fi


#Initial variables/constants

@dataclass
class Body:
    mass: int
    bidx: int
    position: np.ndarray = np.zeros(2)
    velocity: np.ndarray = np.zeros(2)

#Gravitational constant
G = 4*np.pi**2

#m = np.asarray([3.0028e-6, 1, 0]) #small mass, large mass, zero mass
# Initializing 3 bodies of mass and certain velocities
body = []
#Body 1
body.append(Body(1,1, np.asarray([3.3030197, -0.82771837]), np.asarray([1.587433767, 1.47221479])))
#body.append(Body(1,1, np.asarray([1.0, 0.0]), np.asarray([0.0, 6.1682])))
#Body 2
body.append(Body(1,2, np.asarray([-3.3030197, 0.82771837]), np.asarray([1.587433767, 1.47221479])))
#body.append(Body(1,2, np.asarray([0.0, 0.0]), np.asarray([0.0, 0.0])))
#Body 3
body.append(Body(1,3, np.asarray([0.0, 0.0]), np.asarray([-3.174867535, -2.94442961])))


#Time integration stuff
dt = 0.01
tFinal = 5.0
N = int(tFinal/dt)
t = np.linspace(0, tFinal, N)
body_pos = {"1":[], "2":[], "3":[]}

#Energy variables
Energy = np.zeros([N,3])
Total_energy = np.zeros([N,1])
x1 = []
y1 = []


#Implementation of RK3
ti = 0
while ti < tFinal:
    ti+=dt
    #print(ti)
    #Saving the conditions
    for i in range(3):
        body_pos[str(i+1)].append(body[i].position)

    #Finding evaluation at ti
    k1_v_x = []
    idx_comb = [[0, 1,2], [1,0,2], [2,0,1]]
    for i,j,k in idx_comb:
        masses = [body[j].mass, body[k].mass]
        positions = [body[i].position, body[j].position, body[k].position]
        k1_v_x.append(computeForce(masses, positions))
    k1_v_x.append(body[0].velocity)
    k1_v_x.append(body[1].velocity)
    k1_v_x.append(body[2].velocity)


    #Solving for k2 (ti+dt/3)
    temp_body = copy.deepcopy(body)
    for i in range(3):
        temp_body[i].velocity+=(dt/3)*k1_v_x[i]
    for i in range(3):
        temp_body[i].position+=(dt/3)*k1_v_x[i+3]
    k2_v_x = []
    idx_comb = [[0, 1,2], [1,0,2], [2,0,1]]
    for i,j,k in idx_comb:
        masses = [body[j].mass, body[k].mass]
        positions = [temp_body[i].position, temp_body[j].position, temp_body[k].position]
        k2_v_x.append(computeForce(masses, positions))
    k2_v_x.append(temp_body[0].velocity)
    k2_v_x.append(temp_body[1].velocity)
    k2_v_x.append(temp_body[2].velocity)

    #Solving for k3 (ti+2dt/3)
    temp_body = copy.deepcopy(body)
    for i in range(3):
        temp_body[i].velocity+=(2*dt/3)*k2_v_x[i]
    for i in range(3):
        temp_body[i].position+=(2*dt/3)*k2_v_x[i+3]
    k3_v_x = []
    idx_comb = [[0,1,2], [1,0,2], [2,0,1]]
    for i,j,k in idx_comb:
        masses = [body[j].mass, body[k].mass]
        positions = [temp_body[i].position, temp_body[j].position, temp_body[k].position]
        k3_v_x.append(computeForce(masses, positions))
    k3_v_x.append(temp_body[0].velocity)
    k3_v_x.append(temp_body[1].velocity)
    k3_v_x.append(temp_body[2].velocity)


    #Final step
    for i in range(3):
        body[i].velocity+=(dt/4)*(k1_v_x[i]+3*k2_v_x[i])
    for i in range(3):
        body[i].position+=(dt/4)*(k1_v_x[i+3]+3*k2_v_x[i+3])
    x1.append(body[0].position[0])
    y1.append(body[0].position[1])
    print(body[0].position[0], body[0].position[1])





plt.show()

    


















