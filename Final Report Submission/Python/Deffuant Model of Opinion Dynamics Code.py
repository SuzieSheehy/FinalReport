import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats

"""I've added a section at the bottom with a few functions to try so you don't need to read every single function to know how to use them"""

#Generates a lattices if no lattice is specified, returns the lattice
def lat_gen(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5):
    #Allow user to input a lattice of their own, otherwise generate a lattice of given dimension x by y
    if lattice is None:
        lat = np.random.beta(alpha, beta, size=(x,y))
    else:
        lat = lattice
    return lat

#Shows the lattice as a heatmap
def show_lattice(lattice):
    fig, ax = plt.subplots(1, 1)
    plot0 = ax.imshow(lattice)
    fig.colorbar(plot0, ax = ax)
    ax.set_title("Voter Block")
    plt.show()

#show_lattice(lat_gen(alpha = 0.5, beta = 0.5))

#Deffuant Model, picks two sites and updates them according to Deffuant update rule, returns the updated lattice
#Can take a specified lattice or generate one if none is specified
def def_neigh_convert(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.4, r = 0.25):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)

    #Random site
    i = np.random.randint(0, len(lat))
    j = np.random.randint(0, len(lat[0]))

    #Random site to talk to
    k = np.random.randint(0, len(lat))
    l = np.random.randint(0, len(lat[0]))

    if abs(lat[i,j]-lat[k,l]) <= r:
        temp = lat[i,j]
        lat[i,j] = lat[i,j] + mu*(lat[k,l] - lat[i,j])
        lat[k,l] = lat[k,l] + mu*(temp - lat[k,l])
    return lat

#Picks one site at random, picks a neighbour at random and updates sites based on Deffuant model, returns the updated lattice
#Can take lattice or generate one if none is specified
def neigh_conv(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.4, r = 0.25):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)

    i = np.random.randint(0, len(lat))
    j = np.random.randint(0, len(lat[0]))
    neigh = np.random.randint(0,4)

    #Choosing random neighbour for updating
    if neigh == 0:
        conf_neigh = lat[i-1,j]
    elif neigh == 1:
        conf_neigh = lat[(i+1) % len(lat),j]
    elif neigh == 2:
        conf_neigh = lat[i,j-1]
    elif neigh == 3:
        conf_neigh = lat[i-1,(j+1) % len(lat[0])]
    #Debug help
    else:
        print("ERROR IN NEIGH_CONV")

    if abs(lat[i,j]-conf_neigh) <= r:
        temp = lat[i,j]
        lat[i,j] = lat[i,j] + mu*(conf_neigh - lat[i,j])
        if neigh == 0:
            lat[i-1,j] = conf_neigh + mu*(temp - conf_neigh)
        elif neigh == 1:
            lat[(i+1) % len(lat),j] = conf_neigh + mu*(temp - conf_neigh)
        elif neigh == 2:
            lat[i,j-1] = conf_neigh + mu*(temp - conf_neigh)
        elif neigh == 3:
            lat[i-1,(j+1) % len(lat[0])] = conf_neigh + mu*(temp - conf_neigh)
        #Debug help
        else:
            print("ERROR IN NEIGH_CONV 2")


    return lat

"""This is from before the implementation of a noise model so doesn't take means and standard deviations or use the noisy model, also only uses iteration, no tolerance
where it stops the loop once reached. Consider this function as depricated"""
#Runs the updating function the specified amount of times and returns the lattice based on the model (or function) used
def vote_iter(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.4, r = 0.25, iterations = 10000, model = "Deffuant"):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
    i = 0
    if model == "Deffuant":
        while i < iterations:
                lat = def_neigh_convert(lattice = lat, mu = mu, r = r)
                i += 1
    elif model == "close":
        while i < iterations:
            lat = neigh_conv(lattice = lat, mu = mu, r = r)
            i += 1
    return lat

"""This is from before the implementation of a noise model so doesn't take means and standard deviations or use the noisy model, also only uses iteration, no tolerance
where it stops the loop once reached. Consider this function as depricated"""
#Plots the heatmap of the initial lattice and the updated lattice
def plot_heatmap(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.4, r = 0.25, iterations = 10000, model = "Deffuant"):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
    fig, ax = plt.subplots(2, 1, figsize=(11, 10), layout="constrained")
    plot0 = ax[0].imshow(lat)
    fig.colorbar(plot0, ax = ax[0])
    ax[0].set_title("Initial Voter Block")
    lat1 = vote_iter(lattice = lat, alpha = alpha, beta = beta, mu = mu, r = r, iterations = iterations, model = model)
    plot1 = ax[1].imshow(lat1)
    fig.colorbar(plot1, ax = ax[1])
    ax[1].set_title("Voters after " + str(iterations) + " time steps")
    plt.show()

#Measure the average difference between a site and its neighbours (averaged over all sites) and the variance of the sites
#In the report, the correlation was renamed to the "average difference to neighbours (ADN)" to avoid confusion with the statistical definition of correlation
def correlation_measure(lattice):
    #Correlation calculation, creates an array of the same dimension as the lattice and fills spaces with the absolute difference between neighbours and the site
    #corresponding to the site of the lattice
    corr_lat = np.zeros((len(lattice),len(lattice[0])))
    corr = np.zeros((len(lattice),len(lattice[0])))
    for i in range(len(lattice[:,0])):
        for j in range(len(lattice[0,:])):
            corr_lat[i,j] = (lattice[i-1,j] + lattice[(i+1) % len(lattice), j] + lattice[i,j-1] + lattice[i,(j+1) % len(lattice[0])])/4
            corr[i,j] = abs(corr_lat[i,j]-lattice[i,j])
    #Average over all sites
    expCorr = np.mean(corr)

    #Variance calculation
    latSqu = lattice**2
    expLat = np.mean(lattice)
    expLatSqu = np.mean(latSqu)
    
    return expCorr, (expLatSqu - expLat**2)

#Noisy Deffuant model, same as the normal Deffuant moddel but with r being based on a beta distribution so that it is inherently noisy
#Returns updated lattice
def noise_neigh_convert(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.4, mean = 0.25, stand_dev = 0.1):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)

    #Random site
    i = np.random.randint(0, len(lat))
    j = np.random.randint(0, len(lat[0]))

    #Random site to talk to
    k = np.random.randint(0, len(lat))
    l = np.random.randint(0, len(lat[0]))

    #https://stats.stackexchange.com/questions/316086/distribution-that-has-a-range-from-0-to-1-and-with-peak-between-them#
    #beta distribution reparameterised to have a mean so that it can be directly compared to static r
    phi = (mean*((1-mean)/(stand_dev**2)))-1
    r = np.random.beta(phi*mean, phi*(1-mean))
    
    if abs(lat[i,j]-lat[k,l]) <= r:
        temp = lat[i,j]
        lat[i,j] = lat[i,j] + mu*(lat[k,l] - lat[i,j])
        lat[k,l] = lat[k,l] + mu*(temp - lat[k,l])
    return lat

#Update lattice based on model selected (this time with noisy model included), returns the lattice, correlation, variation and the number of timesteps taken as loop can be cut short to save computing time
#The number of times the correlation/variance/timestep is recorded can be increased or reduced to change compute times too
def vote_iter_corr(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.4, r = 0.25, mean = 0.25, stand_dev = 0.05, tol = 0.00001, iterations = 10000, model = "Deffuant", corr_data = 1000):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
    i = 0
    timestep = []
    corr = []
    var = []

    if model == "close":
        while i < iterations:
            lat = neigh_conv(lattice = lat, mu = mu, r = r)

            #Record specified number of measurements on the lattice
            if i % (int(iterations/corr_data)) == 0 and iterations >= corr_data:
                timestep.append(i)
                corr.append(correlation_measure(lat)[0])
                var.append(correlation_measure(lat)[1])
            elif iterations < corr_data:
                print("ERROR, please ensure iterations is greater than corr_data")

            #If the variance difference is less than a threshold, stop the loop to save compute resources
            if  len(var) > 10:
                if abs(var[-1] - var[-5]) <= tol:
                    break
            i += 1

    #Refer to close model above for breakdown of how the below works
    elif model == "Deffuant":
        while i < iterations:
            lat = def_neigh_convert(lattice = lat, mu = mu, r = r)
            if i % (int(iterations/corr_data)) == 0 and iterations >= corr_data:
                timestep.append(i)
                corr.append(correlation_measure(lat)[0])
                var.append(correlation_measure(lat)[1])
            elif iterations < corr_data:
                print("ERROR, please ensure iterations is greater than corr_data")
            if  len(var) > 10:
                if abs(var[-1] - var[-5]) <= tol:
                    break
            i += 1

    #Refer to close model for breakdown of how the below works
    elif model == "Deff_noise":
        while i < iterations:
            lat = noise_neigh_convert(lattice = lat, mu = mu, mean = mean, stand_dev = stand_dev)
            if i % (int(iterations/corr_data)) == 0 and iterations >= corr_data:
                timestep.append(i)
                corr.append(correlation_measure(lat)[0])
                var.append(correlation_measure(lat)[1])
            elif iterations < corr_data:
                print("ERROR, please ensure iterations is greater than corr_data")
            if  len(var) > 10:
                if abs(var[-1] - var[-5]) <= tol:
                    break
            i += 1
    
    return lat, corr, var, timestep

#Plot the heatmap along with graphs on the measurements as the lattice is iterated on
def plot_heatmap_corr(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.4, r = 0.25, mean = 0.25, stand_dev = 0.05, tol = 0.000001, iterations = 10000, model = "Deffuant", corr_data = 1000):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
    fig, ax = plt.subplots(4, 1, figsize=(11, 10), layout="constrained")
    plot0 = ax[0].imshow(lat)
    fig.colorbar(plot0, ax = ax[0])
    ax[0].set_title("Initial Voter Block")
    lat1, corr1, var1, time = vote_iter_corr(lattice = lat, alpha = alpha, beta = beta, mu = mu, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    plot1 = ax[1].imshow(lat1)
    fig.colorbar(plot1, ax = ax[1])
    ax[1].set_title("Voters after " + str(iterations) + " time steps")
    ax[2].plot(time,corr1)
    ax[2].set_title("Average Difference of Neighbours (ADN)")
    ax[2].set_xlabel("Timestep")
    ax[2].set_ylabel("ADN")
    ax[3].plot(time,var1)
    ax[3].set_title("Variation with final value " + str(var1[-1]))
    ax[3].set_xlabel("Timestep")
    ax[3].set_ylabel("Variance")
    plt.show()

#plot_heatmap_corr(x = 50, y = 50, alpha = 1, beta = 1, mu = 0.4, r = 0.25, mean = 0.25, stand_dev = 0.085, model = "Deffuant", tol = 0.0000000000001, iterations = 8000000, corr_data=5000)
#plot_heatmap_corr(alpha = 1, beta = 1, mu = 0.4, r = 0.6, mean = 0.25, stand_dev = 0.085, model = "Deffuant", tol = 0.0000000000001, iterations = 6000000, corr_data=500)
#plot_heatmap_corr(alpha = 5, beta = 1, mu = 0.8, r = 0.25, iterations = 300000, corr_data=5000)
#plot_heatmap_corr(alpha = 2, beta = 5, mu = 0.8, r = 0.25, iterations = 300000, corr_data=5000)

#The below is for plotting the lattice with two different values of mu to show how mu makes the lattice behaviour speed up or slow down (timestep wise not computing time wise)
def plot_heatmap_corr_mu(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, r = 0.25, mean = 0.25, stand_dev = 0.05, tol = 0.001, iterations = 10000, model = "Deffuant", corr_data = 100):
    lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
    fig, ax = plt.subplots(2, 1, figsize=(11, 10), layout="constrained")
    lat1, corr1, var1, time1 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.1, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat1 done")
    lat2, corr2, var2, time2 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.2, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat2 done")
    lat3, corr3, var3, time3 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.3, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat3 done")
    lat4, corr4, var4, time4 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.4, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat4 done")
    lat5, corr5, var5, time5 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.5, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat5 done")
    lat6, corr6, var6, time6 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.6, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat6 done")
    lat7, corr7, var7, time7 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.7, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat7 done")
    lat8, corr8, var8, time8 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.8, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat8 done")
    lat9, corr9, var9, time9 = vote_iter_corr(lattice = lat.copy(), alpha = alpha, beta = beta, mu = 0.9, r = r, mean = mean, stand_dev = stand_dev, tol = tol, iterations = iterations, model = model, corr_data = corr_data)
    print("lat9 done")
    ax[0].plot(time1,corr1, label = "mu = 0.1")
    ax[0].plot(time2,corr2, label = "mu = 0.2")
    ax[0].plot(time3,corr3, label = "mu = 0.3")
    ax[0].plot(time4,corr4, label = "mu = 0.4")
    ax[0].plot(time5,corr5, label = "mu = 0.5")
    ax[0].plot(time6,corr6, label = "mu = 0.6")
    ax[0].plot(time7,corr7, label = "mu = 0.7")
    ax[0].plot(time8,corr8, label = "mu = 0.8")
    ax[0].plot(time9,corr9, label = "mu = 0.9")
    ax[0].set_title("Average Difference to Neighbours (ADN)")
    ax[0].set_xlabel("Timestep")
    ax[0].set_ylabel("ADN")
    ax[0].legend()
    ax[1].plot(time1,var1, label = "mu = 0.1")
    ax[1].plot(time2,var2, label = "mu = 0.2")
    ax[1].plot(time3,var3, label = "mu = 0.3")
    ax[1].plot(time4,var4, label = "mu = 0.4")
    ax[1].plot(time5,var5, label = "mu = 0.5")
    ax[1].plot(time6,var6, label = "mu = 0.6")
    ax[1].plot(time7,var7, label = "mu = 0.7")
    ax[1].plot(time8,var8, label = "mu = 0.8")
    ax[1].plot(time9,var9, label = "mu = 0.9")
    ax[1].set_title("Variation")
    ax[1].set_xlabel("Timestep")
    ax[1].set_ylabel("Variance")
    ax[1].legend()
    plt.show()

"""Defining clusters is not a straight forward task, typical values sites take longterm for Deffuant are ~0.85, ~0.15 and ~0.5 for r = 0.25.
This behaviour originates from the r value used (0.25 by default), meaning after a while no changes can occur as sites are outside eachother's r range
Will use this for defining clusters as being within ranges of something like 0-1/3, 1/3-2/3, 2/3-1
It is worth noting that this very much a manual way of defining things, I might rework this later so that r values define cluster ranges"""

"""I've removed the old implementation that was only useful for r = 0.25, now made it work for r values such that (1/r)-1 ~ integer if r < 1/3
 and just split it into clusters of size 0.2 otherwise as I don't expect clustering for values above r = 0.5 and phase transitions have been occuring for values ~0.38"""

def cluster_lat_dev(lattice, r):

    #Calculate the number of values the lattice will take
    theo_num_clust = round(np.reciprocal(r)) - 1
    clust_arr = np.zeros((len(lattice),len(lattice[0])))
    for i in range(len(lattice[:,0])):
        for j in range(len(lattice[0,:])):

            if r <= 0.34:
                for k in range(theo_num_clust):
                    if (k/theo_num_clust) <= lattice[i,j] < ((k+1)/theo_num_clust):
                        clust_arr[i,j] = k+1
            else:
                for l in range(5):
                    if (l/5) <= lattice[i,j] < ((l+1)/5):
                        clust_arr[i,j] = l+1


    #I found a scipy function, scipy.ndimage.measurements.label, that is able to produce the number of clusters and label them
    #https://stackoverflow.com/questions/25664682/how-to-find-cluster-sizes-in-2d-numpy-array
    #Technically the labelling isn't done above, it simply sorts values into groups that I've defined, hence the use of the above
    #https://github.com/ussserrr/percolation-python/blob/master/README.md
    #The use of ndimages didn't work as it only sorts into two groups, so a different implementation is used
    
    label_clust = np.zeros_like(clust_arr, dtype=int)

    neighbours = [(-1, 0), (1, 0), (0, -1), (0, 1)]

    current_label = 0

    for i in range(len(clust_arr)):
        for j in range(len(clust_arr[0])):

            #If already assigned to a clump, go to next site/cell
            if label_clust[i, j] != 0:
                continue

            #New label for new cluster
            current_label += 1
            value = clust_arr[i, j]

            # Start a new clump
            label_clust[i, j] = current_label
            #List of elements that need investigating
            stack = [(i, j)]

            #while stack is used so that this repeats until stack is empty
            while stack:
                #Remove element from stack once its investigated
                i1, j1 = stack.pop()

                for i_coord, j_coord in neighbours:
                    #Neighbour coordinates
                    neigh_i = (i1 + i_coord) % len(clust_arr)
                    neigh_j = (j1 + j_coord) % len(clust_arr[0])

                    #Check if neighbour cells are already labelled and whether neighbour has same value as original cell
                    if label_clust[neigh_i, neigh_j] == 0 and clust_arr[neigh_i, neigh_j] == value:
                        #Label and append neighbour to stack for later investigation
                        label_clust[neigh_i, neigh_j] = current_label
                        stack.append((neigh_i, neigh_j))
    
    num_clust = current_label
    
    values, area_clust = np.unique(label_clust, return_counts=True)
    
    return label_clust, num_clust, area_clust, clust_arr

def plot_lab_clust(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.8, r = 0.25, mean = 0.25, stand_dev = 0.05, tol = 0.0000001, iterations = 1000, model = "Deffuant", corr_data = 100):
    InitLat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
    TestLat, TestCorr, TestVar, TestTime = vote_iter_corr(lattice = InitLat.copy(), mu = mu, r = r, mean = mean, stand_dev = stand_dev, tol = tol, model = model, iterations = iterations, corr_data = corr_data)
    fig, ax = plt.subplots(1, 1, figsize=(11, 10), layout="constrained")
    ax.set_title("Lattice with Labelled Clusters")
    plot2 = ax.imshow(cluster_lat_dev(lattice = TestLat, r = r)[0])
    fig.colorbar(plot2, ax = ax)
    plt.show()

def plot_clust(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.8, r = 0.25, mean = 0.25, stand_dev = 0.05, tol = 0.0000001, iterations = 1000, model = "Deffuant", corr_data = 100):
    InitLat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
    TestLat, TestCorr, TestVar, TestTime = vote_iter_corr(lattice = InitLat.copy(), mu = mu, r = r, mean = mean, stand_dev = stand_dev, tol = tol, model = model, iterations = iterations, corr_data = corr_data)
    fig, ax = plt.subplots(5, 1, figsize=(11, 10), layout="constrained")
    ax[0].set_title("(a) Initial Lattice")
    plot0 = ax[0].imshow(InitLat)
    fig.colorbar(plot0, ax = ax[0])
    ax[1].set_title("(b) Final Lattice")
    plot1 = ax[1].imshow(TestLat, vmin = 0, vmax = 1)
    fig.colorbar(plot1, ax = ax[1])
    ax[2].set_title("(c) Lattice with Labelled Clusters")
    plot2 = ax[2].imshow(cluster_lat_dev(lattice = TestLat, r = r)[0])
    fig.colorbar(plot2, ax = ax[2])
    ax[3].set_title("(d) Area of Clusters")
    ax[3].set_xlabel("Area")
    ax[3].set_ylabel("Frequency")

    #https://stackoverflow.com/questions/18325706/how-to-plot-a-scatter-plot-using-the-histogram-output-in-matplotlib
    plot1, bins, patches = ax[3].hist(cluster_lat_dev(lattice = TestLat, r = r)[2], bins = np.max(cluster_lat_dev(lattice = TestLat, r = r)[2]))
    bins_mean = [0.5 * (bins[i] + bins[i+1]) for i in range(len(plot1))]
    ax[4].scatter(bins_mean, plot1)
    #https://stackoverflow.com/questions/18773662/scatter-plot-logarithmic-scale
    ax[4].set_yscale('log')
    ax[4].set_xscale('log')
    ax[4].set_xlabel("Area")
    ax[4].set_ylabel("Frequency")
    plt.show()

#Plot the ADN, variance, number of clusters and area of clusters over multiple sims. This function could use multithreading extremely well as each thread could run a few sims
def macro_stats(lattice = None, x = 50, y = 50, alpha = 0.5, beta = 0.5, mu = 0.8, stand_dev = 0.01, tol = 0.00001, iterations = 10000000, corr_data = 1000, sims = 10, r_values = [0.1, 0.2, 0.3, 0.4, 0.5]):
    #Define another function that produces the final ADN, variance, number of clusters and area of clusters
    def corr_var():
        corr_close_arr = []
        var_close_arr = []
        num_clust_close = []
        area_clust_close = []
        corr_deff_arr = []
        var_deff_arr = []
        num_clust_deff = []
        area_clust_deff = []
        corr_noise_arr = []
        var_noise_arr = []
        num_clust_noise = []
        area_clust_noise = []
        lat = lat_gen(lattice = lattice, x = x, y = y, alpha = alpha, beta = beta)
        for r in r_values:
            lat_close, corr_close, var_close, timestep_close = vote_iter_corr(lattice = lat, mu = mu, r = r, tol = tol, iterations = iterations, model = "close", corr_data = corr_data)
            #Due to how the iterations can stop prematurely, I can't vectorise the problem
            corr_close_arr.append(corr_close[-1])
            var_close_arr.append(var_close[-1])
            num_clust_close.append(cluster_lat_dev(lat_close, r)[1])
            area_clust_close.append(np.mean(cluster_lat_dev(lat_close, r)[2]))
            lat_deff, corr_deff, var_deff, timestep_deff = vote_iter_corr(lattice = lat, mu = mu, r = r, tol = tol, iterations = iterations, model = "Deffuant", corr_data = corr_data)
            corr_deff_arr.append(corr_deff[-1])
            var_deff_arr.append(var_deff[-1])
            num_clust_deff.append(cluster_lat_dev(lat_deff, r)[1])
            area_clust_deff.append(np.mean(cluster_lat_dev(lat_deff, r)[2]))
            lat_noise, corr_noise, var_noise, timestep_noise = vote_iter_corr(lattice = lat, mu = mu, mean = r, stand_dev = stand_dev, tol = tol, iterations = iterations, model = "Deff_noise", corr_data = corr_data)
            corr_noise_arr.append(corr_noise[-1])
            var_noise_arr.append(var_noise[-1])
            num_clust_noise.append(cluster_lat_dev(lat_noise, r)[1])
            area_clust_noise.append(np.mean(cluster_lat_dev(lat_noise, r)[2]))
            print("r value " + str(r) + " done")
        return corr_close_arr, var_close_arr, num_clust_close, area_clust_close, corr_deff_arr, var_deff_arr, num_clust_deff, area_clust_deff, corr_noise_arr, var_noise_arr, num_clust_noise, area_clust_noise

    #Run the above function multiple times and average over it to help eliminate variations
    corr_close_arr_res = []
    var_close_arr_res = []
    num_clust_close_res = []
    area_clust_close_res = []
    corr_deff_arr_res = []
    var_deff_arr_res = []
    num_clust_deff_res = []
    area_clust_deff_res = []
    corr_noise_arr_res = []
    var_noise_arr_res = []
    num_clust_noise_res = []
    area_clust_noise_res = []
    for i in range(sims):
        corr_close_arr, var_close_arr, num_clust_close, area_clust_close, corr_deff_arr, var_deff_arr, num_clust_deff, area_clust_deff, corr_noise_arr, var_noise_arr, num_clust_noise, area_clust_noise = corr_var()
        #This can technically be vectorised to reduce compute time but the function works as is
        corr_close_arr_res.append(corr_close_arr)
        var_close_arr_res.append(var_close_arr)
        num_clust_close_res.append(num_clust_close)
        area_clust_close_res.append(area_clust_close)
        corr_deff_arr_res.append(corr_deff_arr)
        var_deff_arr_res.append(var_deff_arr)
        num_clust_deff_res.append(num_clust_deff)
        area_clust_deff_res.append(area_clust_deff)
        corr_noise_arr_res.append(corr_noise_arr)
        var_noise_arr_res.append(var_noise_arr)
        num_clust_noise_res.append(num_clust_noise)
        area_clust_noise_res.append(area_clust_noise)
        print("sim " + str(i) + " done")
    corr_close_arr_av = np.mean(corr_close_arr_res, axis = 0)
    var_close_arr_av = np.mean(var_close_arr_res, axis = 0)
    num_clust_close_av = np.mean(num_clust_close_res, axis = 0)
    area_clust_close_av = np.mean(area_clust_close_res, axis = 0)
    corr_deff_arr_av = np.mean(corr_deff_arr_res, axis = 0)
    var_deff_arr_av = np.mean(var_deff_arr_res, axis = 0)
    num_clust_deff_av = np.mean(num_clust_deff_res, axis = 0)
    area_clust_deff_av = np.mean(area_clust_deff_res, axis = 0)
    corr_noise_arr_av = np.mean(corr_noise_arr_res, axis = 0)
    var_noise_arr_av = np.mean(var_noise_arr_res, axis = 0)
    num_clust_noise_av = np.mean(num_clust_noise_res, axis = 0)
    area_clust_noise_av = np.mean(area_clust_noise_res, axis = 0)

    #Plot the results
    fig, ax = plt.subplots(5, 1, figsize=(11, 10), layout="constrained")
    #https://stackoverflow.com/questions/30577439/scipy-and-the-beta-distribution-is-is-possible-to-change-the-range-of-x-values
    
    x_val = np.linspace(0, 1, 1000)
    y_val = stats.beta.pdf(x_val, alpha, beta)
    ax[0].plot(x_val, y_val, label = "alpha = " + str(alpha) + " and beta = " + str(beta))
    ax[0].legend()
    ax[0].set_title("(a) Beta Distribution for Initial Lattice")
    ax[0].set_xlabel("x")
    ax[0].set_ylabel("Probability Density")

    ax[1].plot(r_values,corr_close_arr_av, label = "Close Model")
    ax[1].plot(r_values,corr_deff_arr_av, label = "Deffuant Model")
    ax[1].plot(r_values,corr_noise_arr_av, label = "Deffuant with Noise")
    ax[1].legend()
    ax[1].set_title("(b) Final ADN for 3 Different Models")
    ax[1].set_xlabel("r Values")
    ax[1].set_ylabel("ADN")
    
    ax[2].plot(r_values,var_close_arr_av, label = "Close Model")
    ax[2].plot(r_values,var_deff_arr_av, label = "Deffuant Model")
    ax[2].plot(r_values,var_noise_arr_av, label = "Deffuant with Noise")
    ax[2].legend()
    ax[2].set_title("(c) Final Variation for 3 Different Models")
    ax[2].set_xlabel("r Values")
    ax[2].set_ylabel("Variance")

    ax[3].plot(r_values,num_clust_close_av, label = "Close Model")
    ax[3].plot(r_values,num_clust_deff_av, label = "Deffuant Model")
    ax[3].plot(r_values,num_clust_noise_av, label = "Deffuant with Noise")
    ax[3].legend()
    ax[3].set_title("(d) Average Number of Clusters for 3 Different Models")
    ax[3].set_xlabel("r Values")
    ax[3].set_ylabel("Number of Clusters")

    ax[4].plot(r_values,area_clust_close_av, label = "Close Model")
    ax[4].plot(r_values,area_clust_deff_av, label = "Deffuant Model")
    ax[4].plot(r_values,area_clust_noise_av, label = "Deffuant with Noise")
    ax[4].legend()
    ax[4].set_title("(e) Average Area of Clusters for 3 Different Models")
    ax[4].set_xlabel("r Values")
    ax[4].set_ylabel("Area of Clusters")
    
    plt.show()

#Plot the beta distribution for different alpha and beta values
def plot_beta_dist():
    fig, ax = plt.subplots(1, 1)
    #https://stackoverflow.com/questions/30577439/scipy-and-the-beta-distribution-is-is-possible-to-change-the-range-of-x-values
        
    x_val = np.linspace(0, 1, 1000)
    y_val1 = stats.beta.pdf(x_val, 0.5, 0.5)
    y_val2 = stats.beta.pdf(x_val, 1, 1)
    y_val3 = stats.beta.pdf(x_val, 5, 1)
    y_val4 = stats.beta.pdf(x_val, 1, 5)
    y_val5 = stats.beta.pdf(x_val, 2, 2)
    ax.plot(x_val, y_val1, label = "alpha = 0.5 and beta = 0.5")
    ax.plot(x_val, y_val2, label = "alpha = 1 and beta = 1")
    ax.plot(x_val, y_val3, label = "alpha = 5 and beta = 1")
    ax.plot(x_val, y_val4, label = "alpha = 1 and beta = 5")
    ax.plot(x_val, y_val5, label = "alpha = 2 and beta = 2")
    ax.legend()
    ax.set_title("Beta Distribution for Initial Lattice")
    ax.set_xlabel("x")
    ax.set_ylabel("Probability Density")

    plt.show()


"""Below are some functions that are ready to try out, just remove the # to run"""

###############################################
###############################################
###############################################

"""Plots the beta distribution for different combinations of alpha and beta. Doesn't take any arguments"""
#plot_beta_dist()

"""Shows a 2D lattice on a heatmap"""
show_lattice(lat_gen(alpha = 0.5, beta = 0.5))

"""Shows the initial and final lattice on a heatmap"""
"""Can take a lattice x and y values, alpha and beta values for the initial lattice generation, mu, r, iterations to run and a model, either "Deffuant" or "close"."""
#plot_heatmap(x = 100, y = 100, r = 0.25, mu = 0.4, iterations = 5000000)

"""Plots the initial/final lattices, ADN and variance over iterations"""
"""Can take a lattice x and y values, alpha and beta values for the initial lattice generation, mu, r, mean if using the noisy model, stand_dev if using noisy model, iterations to run and a model"""
"""Models to choose from are "close", "Deffuant", "Deff_noise"."""
#plot_heatmap_corr(x = 50, y = 50, alpha = 1, beta = 1, mu = 0.4, r = 0.25, mean = 0.25, stand_dev = 0.085, model = "Deffuant", tol = 0.0000000000001, iterations = 8000000, corr_data=5000)
#plot_heatmap_corr(alpha = 1, beta = 1, mu = 0.4, r = 0.6, mean = 0.25, stand_dev = 0.085, model = "Deff_noise", tol = 0.0000000000001, iterations = 6000000, corr_data=5000)

"Plots the variance of multiple lattices as they're iterated on, with multiple lattices with different mu values being plotted on the same plot"
"""Takes lattice, x, y, alpha, beta, r, mean, stand_dev, tol, iterations, model, corr_data"""
"""Models to choose from are "close", "Deffuant", "Deff_noise"."""
#plot_heatmap_corr_mu(alpha = 1, beta = 1, r = 0.25, mean = 0.25, stand_dev = 0.085, model = "Deffuant", tol = 0.000000000000001, iterations = 200000, corr_data=5000)

"""Plots just the labelled-cluster lattice"""
"""Takes lattice, x, y, alpha, beta, r, mean, stand_dev, tol, iterations, model, corr_data"""
"""Models to choose from are "close", "Deffuant", "Deff_noise"."""
#plot_lab_clust(alpha = 1, beta = 1, mu = 0.4, r = 0.25, iterations = 1000000, model = "Deffuant")

"""Plots the lattice and the labelled-cluster lattice, along with plotting the area of the clusters on both a histogram and a log-log scatterplot"""
"""Takes lattice, x, y, alpha, beta, r, mean, stand_dev, tol, iterations, model, corr_data"""
"""Models to choose from are "close", "Deffuant", "Deff_noise"."""
#plot_clust(iterations = 600000, model = "Deffuant")
#plot_clust(alpha = 1, beta = 1, mu = 0.4, r = 0.25, iterations = 1000000, model = "Deffuant")

"""Be warned, this function can take a long time to run, I've lowered some of the parameters so it isn't as bad but for the report different values were used.
The run time is still very long either way though."""
"""Plots the final ADN, variance, number of clusters and area of clusters for different r values and for different models, using the same initial lattice per run
and averaging over the runs (or sims as I've called them) to reach the final plots"""
"""Takes lattice, x, y, alpha, beta, r, mean, stand_dev, tol, iterations, model, corr_data"""
"""Models to choose from are "close", "Deffuant", "Deff_noise"."""
r_elem = np.linspace(0.05,0.95, num = 50)
#macro_stats(x = 30, y = 30, alpha = 1, beta = 1, mu = 0.4, tol = 0.00000000000001, r_values = r_elem, corr_data = 10000, sims = 3)