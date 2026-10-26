import numpy as np
import matplotlib.pyplot as plt
from time_step import time_step, time_step2
from tests import avgdist_l1
from matplotlib.animation import FuncAnimation, PillowWriter

rng = np.random.default_rng(seed=14142)

#Generates 40x40 intial configuration
N = 40
opinions = rng.random((N,N))
logOdds = np.log(opinions/(1-opinions))



"""
Produces pictures for unbiased regime
"""

#Heatmap and histogram at t=0
copy = logOdds.copy()
plt.imshow(copy)
plt.colorbar()
plt.savefig("coda_1")
plt.clf()
counts, bins = np.histogram(copy, bins = 100)
plt.stairs(counts, bins)
plt.savefig("coda_hist_1")
plt.clf()

for j in range(500_000):  #Evolves 500,000 time steps
    time_step(copy, a = 3/4, d = 3/4, b = 3/4, c = 3/4, rng = rng)
    if (j+1) % 100_000 == 0:  #Every 100,000 time steps saves heatmap and histogram
        plt.imshow(copy)
        plt.colorbar()
        plt.savefig(f"coda_{int((j+1)/100000)+1}")
        plt.clf()
        counts, bins = np.histogram(copy, bins = 100)
        plt.stairs(counts, bins)
        plt.savefig(f"coda_hist_{int((j+1)/100000)+1}")
        plt.clf()

"""
Produces pictures for four regimes
"""


#Favouring negative reigeme
copy = logOdds.copy()
for j in range(1_000_000):
    time_step(copy, a = 19/20, d = 11/20, b = 11/20, c = 19/20, rng = rng)
plt.imshow(copy)
plt.colorbar()
plt.savefig("negative")
plt.clf()
counts, bins = np.histogram(copy, bins = 100)
plt.stairs(counts, bins)
plt.savefig("negative_hist")
plt.clf()

#Favouring disagreement reigeme
copy = logOdds.copy()
for j in range(1_000_000):
    time_step(copy, a = 19/20, d = 11/20, b = 19/20, c = 11/20, rng = rng)
plt.imshow(copy)
plt.colorbar()
plt.savefig("disagree")
plt.clf()
counts, bins = np.histogram(copy, bins = 100)
plt.stairs(counts, bins)
plt.savefig("disagree_hist")
plt.clf()

#Favours Positive
copy = logOdds.copy()
for j in range(1_000_000):
    time_step(copy, a = 11/20, d = 19/20, b = 19/20, c = 11/20, rng = rng)
plt.imshow(copy)
plt.colorbar()
plt.savefig("positive")
plt.clf()
counts, bins = np.histogram(copy, bins = 100)
plt.stairs(counts, bins)
plt.savefig("positive_hist")
plt.clf()

#Favours agreement
copy = logOdds.copy()
for j in range(1_000_000):
    time_step(copy, a = 11/20, d = 19/20, b = 11/20, c = 19/20, rng = rng)
plt.imshow(copy)
plt.colorbar()
plt.savefig("agree")
plt.clf()
counts, bins = np.histogram(copy, bins = 100)
plt.stairs(counts, bins)
plt.savefig("agree_hist")
plt.clf()





"""
Produces variant model histograms
"""


#Variant Model
copy = logOdds.copy()
counts, bins = np.histogram(copy, bins = 100)
plt.stairs(counts, bins)
plt.savefig("coda_var_hist_1")
plt.clf()
for j in range(200_000):
    time_step2(copy, rng = rng)
    if (j+1) % 40_000 == 0:  #Every 40,000 time steps records a histogram
        counts, bins = np.histogram(copy, bins = 100)
        plt.stairs(counts, bins)
        plt.savefig(f"coda_var_hist_{int((j+1)/40_000)+1}")
        plt.clf()






"""
Produces average distance to disagreement statistic graph
"""


#Avg dist statistic   
points = 10  #how many values of t will be sampled
dist = np.zeros((points,10))  #Initialise list of average distances
for i in range(10):  #Averaged over 10 initial configurations
    opinions = rng.random((N,N))
    logOdds = np.log(opinions/(1-opinions))
    for t in range(points):
        a = t/(2*points) + 1/2 + (1/2 - (points-1)/(2*points))/2  #Instead of starting at t=0 we start a little bit inside as the time_step function has log(0) errors at this point
        d = 1 - (a - 1/2)
        b = max(d, a)
        c = min(d, a)
        copy = logOdds.copy()
        for j in range(100_000):  
            time_step(copy, a = a, d = d, b = b, c = c, rng = rng)
            if (j+1) % 10_000 == 0:  #Every 10,000 timesteps checks if all vertices have the same sign and breaks early if so
                n, m = copy.shape
                sign = np.sign(copy[n-1,m-1])
                agreement = True
                x = 0
                while agreement and x < n:
                    for y in range(m):
                        if np.sign(copy[x,y]) != sign:
                            agreement = False
                            break
                    x += 1
                if agreement:
                    break
                
        dist[t,i] = avgdist_l1(copy) 

plt.plot(np.arange(points)/points + (1/2 - (points-1)/(2*points))/2, np.sum(dist,1)/10)  #Plots average distance averaged over the 10 intial configurations
plt.xlabel("t")
plt.ylabel("Average Distance")
plt.savefig("avgDist2")
plt.clf()

"""
Produces average distance to disagreement statistic graph for larger time scale as well as a gif of how the final configurations change along the curve
"""

#Avg dist statistic larger timescale
points = 10
dist = np.zeros(points)
states = []  #List of states for gif
for t in range(points):
    a = t/(2*points) + 1/2 + (1/2 - (points-1)/(2*points))/2
    d = 1 - (a - 1/2)
    b = max(d, a)
    c = min(d, a)
    copy = logOdds.copy()
    for j in range(1_000_000):  #Larger number of total timesteps
        time_step(copy, a = a, d = d, b = b, c = c, rng = rng)
        if (j+1) % 100000 == 0:
            n, m = copy.shape
            sign = np.sign(copy[n-1,m-1])
            agreement = True
            x = 0
            while agreement and x < n:
                for y in range(m):
                    if np.sign(copy[x,y]) != sign:
                        agreement = False
                        break
                x += 1
            if agreement:
                break
                        
                
    dist[t] = avgdist_l1(copy)
    states.append(copy.copy())


fig, ax = plt.subplots()

plot = ax.imshow(states[0])

#Creates gif where each frame is the heatmap of a state in the states list
ani = FuncAnimation(fig=fig, func = lambda t : ax.imshow(states[t]), frames=len(states), interval=60)

ani.save("avgdis.gif", writer=PillowWriter(fps=2))
plt.clf()


plt.plot(np.arange(points)/points + (1/2 - (points-1)/(2*points))/2, dist)
plt.xlabel("t")
plt.ylabel("Average Distance")
plt.savefig("avgDist")
plt.clf()

