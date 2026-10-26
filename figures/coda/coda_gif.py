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
Favouring negative gif
"""


copy = logOdds.copy()
states = [copy.copy()]
for j in range(1_000_000):
    time_step(copy, a = 19/20, d = 11/20, b = 11/20, c = 19/20, rng = rng)
    if (j+1) % 10_000 == 0:
        states.append(copy.copy())

fig, ax = plt.subplots()
counts, bins = np.histogram(states[0], bins = 100)
plot = ax.stairs(counts, bins)

def func(t):
    counts, bins = np.histogram(states[t], bins = 100)
    ax.cla()
    return ax.stairs(counts, bins)

ani = FuncAnimation(fig=fig, func = func, frames=len(states), interval=60)

ani.save("negative_hist.gif", writer=PillowWriter(fps=20))
plt.clf()



"""
Produces favouring diagreement gif
"""


copy = logOdds.copy()
states = [copy.copy()]
for j in range(1_000_000):
    time_step(copy, a = 19/20, d = 11/20, b = 19/20, c = 11/20, rng = rng)
    if (j+1) % 10_000 == 0:  #Every 10,000 time steps adds a configuration to the list of states
        states.append(copy.copy())

fig, ax = plt.subplots()
counts, bins = np.histogram(states[0], bins = 100)
plot = ax.stairs(counts, bins)

def func(t):  #Function which sends the t frame to the plot of the t state
    counts, bins = np.histogram(states[t], bins = 100)
    ax.cla()
    return ax.stairs(counts, bins)

ani = FuncAnimation(fig=fig, func = func, frames=len(states), interval=60)

ani.save("disagreement_hist.gif", writer=PillowWriter(fps=20))
plt.clf()



"""
Produces favouring positive gif
"""


copy = logOdds.copy()
states = [copy.copy()]
for j in range(1_000_000):
    time_step(copy, a = 11/20, d = 19/20, b = 19/20, c = 11/20, rng = rng)
    if (j+1) % 10_000 == 0:
        states.append(copy.copy())

fig, ax = plt.subplots()
counts, bins = np.histogram(states[0], bins = 100)
plot = ax.stairs(counts, bins)

def func(t):
    counts, bins = np.histogram(states[t], bins = 100)
    ax.cla()
    return ax.stairs(counts, bins)

ani = FuncAnimation(fig=fig, func = func, frames=len(states), interval=60)

ani.save("positive_hist.gif", writer=PillowWriter(fps=20))
plt.clf()


"""
Produces favouring agreement gif
"""


copy = logOdds.copy()
states = [copy.copy()]
for j in range(1_000_000):
    time_step(copy, a = 11/20, d = 19/20, b = 11/20, c = 19/20, rng = rng)
    if (j+1) % 10_000 == 0:
        states.append(copy.copy())

fig, ax = plt.subplots()
counts, bins = np.histogram(states[0], bins = 100)
plot = ax.stairs(counts, bins)

def func(t):
    counts, bins = np.histogram(states[t], bins = 100)
    ax.cla()
    return ax.stairs(counts, bins)

ani = FuncAnimation(fig=fig, func = func, frames=len(states), interval=60)

ani.save("agreement_hist.gif", writer=PillowWriter(fps=20))
plt.clf()
