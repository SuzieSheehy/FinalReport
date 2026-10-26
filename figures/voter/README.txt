Dependencies:
- matplotlib 3.11.2
- numpy 2.5.4

Running voter.py produces all figures in the voter model chapter of the report, figures 2.1-2.4.
Producing figures 2.1,2.2 and 2.4 should not take long to produce however figure 2.3 may take a few minutes.
Additionally gif.py may be run to produce a gif (voter.gif) of a 65x65 configuration evolving till stationary with a threshold of 2 for ones enjoyment.

Extra Notes on Functions:

voter_functions.py contains functions which give a single time step for the voter model processes as well as a checker for it a state is a stationar state.
Our functions allow the specification of a neighbourhood D, which we write as tuple of length 2 lists.
This neighbourhood allows you to edit how you want the integer torus to be connected.
For example one could specify D = ([0,0], [1,0], [0,1], [-1,0], [0,-1], [1,1], [1,-1], [-1,1], [-1,-1]), which would mean that every vertex on the lattice is connected to its nearest neighbours along the vertical and horozontal lines as usual, as well as its neighbours along the diagonal lines and to itself.
We also allow the specification of a probability p for the time step function.
The purpose of this naively the function would run forever if a stationary state is entered, so if enough random samples of the configuration have occured such that we have a p chance of being in a stationary state it breaks out of the function.
