import matplotlib.pyplot as plt
import numpy as np

#Sup metric  (Didnt actually end up using this)
def avgdist(opinions):
    """
    Takes a configuration(opinions) and computes the average shortest distance from a vertex to another vertex which disagrees with it in the "sup-metric" on the integer torus lattice 

    Input:
    opinions - 2 dimensional numpy array

    Output:
    float
    """
    
    n, m = opinions.shape
    cum = 0  #Cummulative shortest distance
    for x in range(n):
        for y in range(m):
            r = 1  #Radius of ball being checked
            found = False  #Boolean for finding disagreement
            while not found:
                for u in range(r+1):  #Iterates over the sphere
                    if np.sign(opinions[x,y]) != np.sign(opinions[(x+r)%n,(y+u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-r)%n,(y+u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x+r)%n,(y-u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-r)%n,(y-u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x+u)%n,(y+r)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-u)%n,(y+r)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x+u)%n,(y-r)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-u)%n,(y-r)%m]):
                        found = True
                        break
                r += 1
                if r == max(n,m)+1:
                    return max(n,m)+1
            cum += r
    return cum/(n*m)


#l1 metric
def avgdist_l1(opinions):
    """
    Takes a configuration(opinions) and computes the average shortest distance from a vertex to another vertex which disagrees with it in the "l1-metric" on the integer torus lattice 

    Input:
    opinions - 2 dimensional numpy array

    Output:
    float
    """
    n, m = opinions.shape
    cum = 0  #Cummulative shortest distance
    for x in range(n):
        for y in range(m):
            r = 1  #Radius of spherical shell being checked
            found = False  #Boolean for finding disagreement
            while not found:
                for u in range(r+1):  #Iterates over the sphere
                    if np.sign(opinions[x,y]) != np.sign(opinions[(x+r-u)%n,(y+u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-r+u)%n,(y+u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x+r-u)%n,(y-u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-r+u)%n,(y-u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x+u)%n,(y+r-u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-u)%n,(y+r-u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x+u)%n,(y-r+u)%m]):
                        found = True
                        break
                    elif np.sign(opinions[x,y]) != np.sign(opinions[(x-u)%n,(y-r+u)%m]):
                        found = True
                        break
                r += 1
                if r == n+m+1:
                    return n+m+1
            cum += r
    return cum/(n*m)
