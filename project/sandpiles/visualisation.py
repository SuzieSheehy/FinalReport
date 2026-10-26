import matplotlib.pyplot as plt
import numpy as np
from typing import Callable

# Plot the empirical pmf from a simulation statistic series
def visualise_statistic(
        stat_series: np.ndarray[int],       # series to plot
        include_zero: bool = False,         # whether to include 0 (i.e. do we care if there is no topple)
        log: bool = False,                  # whether to log the counts and stats
        ax = None,                          # plot axes for formatting
        marker = 'o',                       # default to a scatter plot
        pmf: Callable = None,               # pmf for MLE
        mle_params: tuple[float] = None,    # fitted MLE params
        pmf_kwargs: dict = {},              # other pmf arguments
        mle_plot_kwargs: dict = {},         # line plot kwargs
        x_max_MLE: int = None,              # optional x_max for MLE
        **kwargs                            # plot keyword arguments for formatting
):
    # get the frequency by size of event
    stat, counts = np.unique(stat_series, return_counts=True)

    # remove 0 if requested
    if not include_zero and stat[0] == 0:
        stat = stat[1:]
        counts = counts[1:]

    freq = counts / np.sum(counts)

    if ax is None:
        # generate the required axes
        ax = plt.gca()

    if log:
        ax.set_xscale('log')
        ax.set_yscale('log')

    ax.plot(stat, freq, marker=marker, **kwargs)

    if pmf is not None:
        if x_max_MLE is None:
            x_max_MLE = max(stat)+1
        x = np.arange(1, x_max_MLE)
        y = pmf(params=mle_params, x=x, loglik=False, **pmf_kwargs)
        ax.plot(x, y, **mle_plot_kwargs)


# Visualise a sandpile on a lattice
def visualise_sandpile(board: np.ndarray[np.int8, np.int8],     # sandpile
                       ax = None                                # existing plot axes (for manual titles etc)
                       ):

    if ax is None:
        fig, ax = plt.subplots()
    else:
        fig = ax.figure
    
    im = ax.imshow(board, vmin=0, vmax = max(4, np.max(board)))

    r, c = board.shape

    # Show all ticks and label them with the respective list entries
    ax.set_xticks(range(c), labels=range(1, c+1),
                rotation=45, rotation_mode="xtick")
    ax.set_yticks(range(r), labels=range(r, 0, -1))

    # Loop over data dimensions and create text annotations.
    for i in range(r):
        for j in range(c):
            text = ax.text(j, i, str(board[i, j]),
                        ha="center", va="center", color="w")

    fig.tight_layout()


# Plots the avalanche size and area graphs for the report - centralises formatting
def visualise_mle_fit(
        example_subset: list[int],      # subset of gridsizes to plot
        avalanche_data,                 # simulation results
        pmf,                            # pmf to use for the MLE
        statistics,                     # list of statistics to plot
        params,                         # optimisation results
        x_max,                          # hard cutoff values
        colors,                         # statistic color mapping
        names,                          # statistic names
        symbols,                        # statistic symbols
        subplots = None,                # existing fig, ax for formatting
        suptitle = None,                # big title to write at the top
        sup_y = None,                   # y position of the big title
):
    if subplots is None:
        fig, ax = plt.subplots(len(example_subset), len(statistics), 
                               figsize=(5*len(statistics), 10*len(example_subset)/3))
    else:
        fig, ax = subplots

    for j, stat in enumerate(statistics):
        for i, L_example in enumerate(example_subset):
            # plot the data using the previous helper function
            visualise_statistic(
                stat_series = avalanche_data[stat][L_example][avalanche_data[stat][L_example] <= x_max[stat][L_example]],
                log=True,
                ax = ax[i, j],
                pmf = pmf,
                mle_params = params[stat][L_example],
                pmf_kwargs = {"x_max": x_max[stat][L_example]},
                mle_plot_kwargs={"color": "black", 
                                "label": fr"MLE: $\hat{{\tau}}_{symbols[stat]} = {params[stat][L_example][0]: .2f}$", 
                                "marker": "o",
                                "ms": 1.5},
                ms = 1.5,
                label = f"L = {L_example}",
                color = colors[stat],
                linestyle="None",
            )
            ax[i, j].axvline(x_max[stat][L_example], linestyle="dashed", color=colors[stat], label=fr"${symbols[stat]}_{{max}}(L)$")
            ax[i, j].legend()
            ax[i, j].set_xlim(0.5, max([2*x_max[stat][L] for L in example_subset]))
            ax[i, j].set_ylim(10**(-6), 10**0)

            # only put a y label on the rows
            ax[i, 0].set_ylabel("Frequency")


        # only title the top and put an xlabel on the bottom of each plot column
        ax[0, j].set_title(f"Distribution of Avalanche {names[stat]}")
        ax[-1, j].set_xlabel(f"Avalanche {names[stat]}")

    if suptitle is not None:
        fig.suptitle(suptitle, y=sup_y)


def visualise_tau_plots(
    statistics,                                 # list of statistics to plot tau for
    gridsizes,                                  # list of all values of L
    params,                                     # optimisation results
    symbols,                                    # dict connecting statistics to notation letters
    colors,                                     # dict connecting statistics to standard colors
    subplots = None,                            # pass in existing subplot
    coef_index = 0,                             # index of param to plot within param array
    title = r"MLE $\hat{\tau}(L)$ vs $L$",      # plot title
    ylabel = r"Power Law Coefficient $\hat{\tau}$ MLE", # plot ylabel
    latex = r"\hat{{\tau}}",                    # latex formatting character to drop in
    hlines = None,                              # horizontal line values for pooled MLE plots
    **kwargs                                    # other arguments
):
    if subplots is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig, ax = subplots
    
    for stat in statistics:
        ax.plot(gridsizes, 
                [params[stat][L][coef_index] for L in gridsizes], 
                label=fr"${latex}_{symbols[stat]}(L)$", 
                color=colors[stat],
                marker = "o",
                ms = 1.5,
                **kwargs)

        if hlines is not None:
            ax.hlines(
                hlines[stat], 
                xmin=0, 
                xmax=gridsizes[-1], 
                colors=[colors[stat]], 
                linestyles=["dashed"], 
                label=fr"${latex}_{symbols[stat]} = {hlines[stat]: .2f}$")

    ax.set_xlabel(r"L")
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=min(0.5, min(params[s][L][coef_index] - 0.1 for s in statistics for L in gridsizes)))
    ax.set_title(title)
    ax.legend()


def visualise_c_L_plots(
    statistics,                                 # list of statistics to plot tau for
    gridsizes,                                  # list of all values of L
    params,                                     # optimisation results
    symbols,                                    # dict connecting statistics to notation letters
    colors,                                     # dict connecting statistics to standard colors
    subplots = None,                            # pass in existing subplot
    title = r"MLE $\hat{c}(L)$ vs $L$",         # plot title
    k = None,                                   # pooled line intercept
    D = None,                                   # pooled fractal dimension
    **kwargs                                    # other arguments
    ):
    if subplots is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig, ax = subplots

    x = np.arange(gridsizes[0], gridsizes[-1] + 1)

    for stat in statistics:
        ax.plot(gridsizes, 
                [params[stat][L][1] for L in gridsizes], 
                label=fr"$\hat{{c}}_{symbols[stat]}(L)$", 
                marker="o", 
                ms=1.5, 
                color = colors[stat],
                **kwargs)

        if k is not None and D is not None:
            y = k[stat] * np.power(x, D[stat])
            ax.plot(x, y, color=colors[stat], label=fr"$\hat{{D}}_{symbols[stat]} = {D[stat]: .2f}$", linestyle="dashed")
        
        ax.legend()
        ax.set_xlabel(r"L")
        ax.set_ylabel(r"Soft Truncation Parameter $\hat{c}$ MLE")
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_title(title)
