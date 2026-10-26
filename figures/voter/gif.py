import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from voter_functions import time_step, stationary_check


rng = np.random.default_rng(seed=141421)

"""
Produces gif of 65x65 configuration evolving in time with threshold = 2
"""

N = 65
opinion = rng.choice([0,1],(N,N))
copy = opinion.copy()
stationary = False

i = 1
t = 500

states = [copy.copy()]  #List of states which will be turned into the gif
while not stationary:  #Evolves until stationary
    stationary = time_step(copy, D = ([0,0],[1,0],[0,1],[-1,0],[0,-1]), T = 2, rng = rng)[1]
    if stationary:
        stationary = stationary_check(copy, T = 2)
        if stationary:
            states.append(copy.copy())
    if i % t == 0:  #Every t time steps records current state to states list
        states.append(copy.copy())
    i += 1

fig, ax = plt.subplots()

plot = ax.imshow(states[0])

#Creates gif where each frame is the heatmap of a state in the states list
ani = FuncAnimation(fig=fig, func = lambda t : ax.imshow(states[t]), frames=len(states), interval=60)

ani.save("voter.gif", writer=PillowWriter(fps=20))