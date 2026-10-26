import numpy as np
import matplotlib.pyplot as plt
from voter_functions import time_step, stationary_check

#Set seed
rng = np.random.default_rng(seed=141421)


"""
10x10  threshold = 3 example
"""


#Generate 10x10 example
N = 10
opinion = rng.choice([0,1],(N,N))
copy = opinion.copy()  #Produces copy of generated intial configuration, this is the one we will evolve
stationary = False

i = 0
t = 4  #Time step size for graphics

#Evolves configuration until it reaches a stationary state
while not stationary:
    stationary = time_step(copy, D = ([0,0],[1,0],[0,1],[-1,0],[0,-1]), T = 3, rng = rng)[1]  #stationary becomes true if it samples for too long without finding a vertex to be updated
    if stationary:
        stationary = stationary_check(copy, T = 3)  #checks if it is in fact in a stationary state
        if stationary:  #if so saves a final image and exists
            plt.imshow(copy)
            plt.savefig(f"Figure_T3_{int(i/t +2)}")
            plt.clf()
    if i % t == 0:      #Saves heatmap of current configuration every t steps
        plt.imshow(copy)
        plt.savefig(f"Figure_T3_{int(i/t +1)}")
        plt.clf()
    i += 1





"""
65x65  threshold = 2 example
"""

#Generate 65x65 example
N = 65
opinion = rng.choice([0,1],(N,N))
copy = opinion.copy()
stationary = False

i = 0
t = 7750

while not stationary:
    stationary = time_step(copy, D = ([0,0],[1,0],[0,1],[-1,0],[0,-1]), T = 2, rng = rng)[1]
    if stationary:
        stationary = stationary_check(copy, T = 2)
        if stationary:
            plt.imshow(copy)
            plt.savefig(f"Figure_T2_{int(i/t +2)}")
            plt.clf()
    if i % t == 0:
        plt.imshow(copy)
        plt.savefig(f"Figure_T2_{int(i/t +1)}")
        plt.clf()
    i += 1





"""
65x65  threshold = 1 example
"""

copy = opinion.copy()
for j in range(t*4+1):
    time_step(copy, D = ([0,0],[1,0],[0,1],[-1,0],[0,-1]), T = 1, rng = rng)
    if j % t == 0:
        plt.imshow(copy)
        plt.savefig(f"Figure_T1_{int(j/t +1)}")
        plt.clf()





"""
Scale invariant plot
"""

for N in [10, 25 , 50]:  #Iterates over grid size
    diff = np.zeros(4)
    samples = 20  #Number of initial configurations beinging avereaged over
    for j in range(samples):
        opinion = rng.choice([0,1],(N,N))
        for hold in range(2,6):  #Iterates over threshold sizes
            stationary = False
            copy = opinion.copy()
            while not stationary:  #Evolves until stationary
                stationary = time_step(copy, D = ([0,0],[1,0],[0,1],[-1,0],[0,-1]), T = hold, rng = rng)[1]
                if stationary:
                    stationary = stationary_check(copy, T = hold)
            diff[hold-2] += np.sum(np.abs(opinion - copy))  #Adds total number of changed opinions from initial configuration
    plt.plot(np.array([2,3,4,5]),diff/(samples*(N**2)), label=f"N={N}")
plt.xlabel(r'$T$')
plt.ylabel(r"$\Delta$ opinions / $N^2$")
plt.legend()
plt.savefig("Avg_Change")
plt.clf()
