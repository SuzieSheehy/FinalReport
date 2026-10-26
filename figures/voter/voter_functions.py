import numpy as np

def time_step(opinion, D = ([1,0],[0,1],[-1,0],[0,-1],[0,0]), T = 1, p = 0.90 , rng = np.random.default_rng()):
    """Takes in a configuration(opinion) and repeatedly samples it until it finds a site to update via the voter model.
    It then returns the new configuration.
    It also returns a boolean which will typically be False.
    If the function is unable to find a vertex to update and believes there is a p change of being in a stationary state it returns True.
    Assumes that the geometry of the process is an integer torus, 
    
    Input:
    opinion - 2 dimensional numpy array
    D - tuple of length 2 lists
    T - non-negative integer
    p - float in the open interval (0,1)
    rng - numpy random generator

    Output:
    2 dimensional array, boolean
    
    """
    
    #Records number of samples without finding a site that can be updated
    samples = 0
    while True:
        samples += 1
        n, m = opinion.shape
        #Samples a coordinate
        x = rng.choice(range(n))
        y = rng.choice(range(m))
        #Counts how many neihbours are of opinion 1
        total = 0
        for w, z in D:
            total += opinion[(x+w)%n,(y+z)%m]
        #Checks update conditions
        if opinion[x,y] == 0 and total >= T:
            opinion[x,y] = 1
            return opinion, False
        elif opinion[x,y] == 1 and len(D) - total >= T:
            opinion[x,y] = 0
            return opinion, False
        #Exits early if there is a probability p chance of being in a stationary state
        elif samples >= np.log(1-p)/np.log(1-1/opinion.size):
            return opinion, True

def stationary_check(opinion, D = ([1,0],[0,1],[-1,0],[0,-1],[0,0]), T = 1):
    """
    Checks if a given configuration(opinion) is a stationary configuration.

    Input:
    opionion - 2 dimensional numpy array
    D - tuple of length 2 lists
    T - non-negative integer   

    Output:
    boolean
    """
    n, m = opinion.shape
    for x in range(n):
        for y in range(m):
            total = 0
            for w, z in D:
                total += opinion[(x+w)%n,(y+z)%m]
            if opinion[x,y] == 0 and total >= T:
                return False
            elif opinion[x,y] == 1 and len(D) - total >= T:
                return False
    return True