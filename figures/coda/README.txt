Dependencies:
- matplotlib 3.11.2
- numpy 2.5.4


Running codaIm.py will produce all figures in the CODA model chapter of the report, figures 3.1-3.6.
Additional it produces a gif (avgdis.gif) of the final configurations that were computed in the making of figure 3.5. This allows you to see how the final configurations change as you move along the curve we consder in equation (3.3.3).
codaIm may take a few minutes to produce all figures. Figures  3.4-3.5 in particular take a couple minutes to produce each.

Additionally running coda_gif.py produces gifs of the evolution of the histograms of configurations in the various extreme regimes presented in figure 3.3. 
These are written as disagreement_hist.gif, agreement_hist.gif, negative_hist.gif, positive_hist.gif.


Extra Notes on Functions:

time_step.py contains functions which give a single time step for the CODA model processes we consider.
Our functions allow the specification of a neighbourhood D, which we write as tuple of length 2 lists.
This neighbourhood allows you to edit how you want the integer torus to be connected.
For example one could specify D = ([1,0], [0,1], [-1,0], [0,-1], [1,1], [1,-1], [-1,1], [-1,-1]), which would mean that every vertex on the lattice is connected to its nearest neighbours along the vertical and horozontal lines as usual, as well as its neighbours along the diagonal lines.

test.py contains a function which produces the average shortest l1 distance from a vertex to a disagreeing vertex as we used in the creation of Figures 3.4-3.5.
It also contains a similar function for the sup metric, also known as the l infinity metric, which was not used in this report.
