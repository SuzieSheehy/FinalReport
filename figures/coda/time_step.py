import numpy as np

#Standard CODA model.
def time_step(opinion, D = ([1,0],[0,1],[-1,0],[0,-1]), a = 3/5, b = 3/5, c = 3/5, d = 3/5, rng = np.random.default_rng()):    
    """
    Takes in a configuration(opinion) and evolves it a single time step via the CODA model.
    Assumes that the geometry is that of the integer torus.

    Input:
    opinion - 2 dimensional numpy array
    D - tuple of length 2 lists
    a,b,c,d - floats in the open interval (0,1)
    rng - numpy random generator

    Output:
    2 dimensional numpy array
    """
    
    n, m = opinion.shape
    #Samples a coordinate
    x = rng.choice(n)
    y = rng.choice(m)
    #Samples a neighbour
    w, z = rng.choice(D)
    #Applys update rule
    if opinion[(x+w)%n,(y+z)%m] >= 0 and opinion[x,y] >= 0:
        opinion[x,y] += np.log(a/(1-d))
    elif opinion[(x+w)%n,(y+z)%m] < 0 and opinion[x,y] >= 0:
        opinion[x,y] += -np.log(d/(1-a))
    elif opinion[(x+w)%n,(y+z)%m] < 0 and opinion[x,y] < 0:
        opinion[x,y] += -np.log(b/(1-c))
    elif opinion[(x+w)%n,(y+z)%m] > 0 and opinion[x,y] < 0:
        opinion[x,y] += np.log(c/(1-b))
    return opinion


#Variant CODA model considered in paper
def time_step2(opinion, D = ([1,0],[0,1],[-1,0],[0,-1]), rng = np.random.default_rng()):
    """
    Takes in a configuration(opinion) and evolves it a single time step via the variant CODA model.
    Assumes that the geometry is that of the integer torus.

    Input:
    opinion - 2 dimensional numpy array
    D - tuple of length 2 lists
    rng - numpy random generator

    Output:
    2 dimensional numpy array
    """

    
    n, m = opinion.shape
    #Samples a coordinate
    x = rng.choice(n)
    y = rng.choice(m)
    #Samples neighbour
    w, z = rng.choice(D)
    #Applys update rule
    opinion[x,y] += opinion[(x+w)%n,(y+z)%m]
    return opinion