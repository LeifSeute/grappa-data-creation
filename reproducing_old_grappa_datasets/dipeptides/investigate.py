#%%
from grappa.data import MolData
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from tqdm import tqdm
import copy
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %%
ff = 'charmm36'

dspath = Path(__file__).parent / 'data' / f'grappa_dipeptides_300K_{ff}'

# %%
gradients = []
ff_gradients = []
energies = []
ff_energies = []
names = []

files = list([f for f in dspath.iterdir() if not f.is_dir()])
for ds in tqdm(files):
    names.append(ds.stem)
    data = MolData.load(ds)
    gradients.append(data.gradient.flatten())
    ff_gradients.append(data.ff_gradient[ff].flatten())
    energies.append(data.energy - data.energy.mean())
    ff_energies.append(data.ff_energy[ff] - data.ff_energy[ff].mean())

# %%


# Concatenate lists of arrays for scatter plots
gradients_concat = np.concatenate(gradients)
ff_gradients_concat = np.concatenate(ff_gradients)
energies_concat = np.concatenate(energies)
ff_energies_concat = np.concatenate(ff_energies)

# Calculate RMSE for each molecule
rmse_gradients = [np.sqrt(np.mean((g - ffg) ** 2)) for g, ffg in zip(gradients, ff_gradients)]
rmse_energies = [np.sqrt(np.mean((e - ffe) ** 2)) for e, ffe in zip(energies, ff_energies)]

# Calculate maximum absolute force values
max_force_qm = [np.max(np.abs(g)) for g in gradients]

max_energy_qm = [np.max(e) - np.min(e) for e in energies]

k = 5
top_k_max_forces = np.argsort(-np.array(max_force_qm))[:k]
top_k_max_force_rmse = np.argsort(-np.array(rmse_gradients))[:k]

print(f'Top {k} molecules with highest QM forces:\n{" ".join([f"{names[i]} ({round(max_force_qm[i], ndigits=1)}), " for i in top_k_max_forces])}\n')
print(f'Top {k} molecules with highest force RMSE:\n{" ".join([f"{names[i]} ({round(rmse_gradients[i], ndigits=1)}), " for i in top_k_max_force_rmse])}\n')


def calculate_density_scatter(x, y, delta_factor=100, seed=0):
    np.random.seed(seed)
    points = []
    frequencies = []

    # Create a deep copy of the input arrays
    x_copy = copy.deepcopy(x)
    y_copy = copy.deepcopy(y)

    # Calculate delta
    delta = max(max(x) - min(x), max(y) - min(y)) / delta_factor

    while len(x_copy) > 0:
        # Pick a random point
        idx = np.random.randint(len(x_copy))
        point_x, point_y = x_copy[idx], y_copy[idx]

        # Calculate distances to the random point
        distances = np.sqrt((x_copy - point_x)**2 + (y_copy - point_y)**2)

        # Find points within the distance delta
        within_delta = distances < delta

        # Count the number of points within delta
        frequency = np.sum(within_delta)

        # Store the point and its frequency
        points.append((point_x, point_y))
        frequencies.append(frequency)

        # Remove the points within delta from the copied arrays
        x_copy = x_copy[~within_delta]
        y_copy = y_copy[~within_delta]

    return np.array(points), np.array(frequencies)


def scatter_plot(ax, x, y, n_max=None, seed=0, symmetric=False, alpha=1., s=15, num_ticks=None, ax_symmetric=False, cluster=False, delta_factor=100, cmap='viridis', logscale=False, **kwargs):
    """
    Create a scatter plot of two arrays x and y.
    Args:
        ax: Matplotlib axis object
        x: Array of x values
        y: Array of y values
        n_max: Maximum number of points to plot
        seed: Random seed for selecting n_max points
        symmetric: Make the plot symmetric around the origin
        alpha: Transparency of the points
        s: Size of the points
        num_ticks: Number of ticks on the axes
        ax_symmetric: Whether to make the axes symmetric
        cluster: Whether to cluster points
        delta_factor: Factor for clustering points
        cmap: Colormap for clustered scatter plot
        logscale: Whether to use a log scale for the colorbar
        **kwargs: Additional keyword arguments for ax.scatter call
    """
    if n_max is not None and n_max < len(x):
        np.random.seed(seed)
        idxs = np.random.choice(len(x), n_max, replace=False)
        x = x[idxs]
        y = y[idxs]

    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    if ax_symmetric:
        min_val = -max(abs(min_val), abs(max_val))
        max_val = max(abs(min_val), abs(max_val))

    if symmetric:
        val = max(abs(min_val), abs(max_val))
        min_val = -val
        max_val = val

    ax.set_ylim(min_val, max_val)
    ax.set_xlim(min_val, max_val)

    ax.set_aspect('equal', 'box')

    ax.plot([min_val, max_val], [min_val, max_val], color='black', linestyle='-', linewidth=1)

    if num_ticks is not None:
        ax.xaxis.set_major_locator(plt.MaxNLocator(num_ticks))
        ax.yaxis.set_major_locator(plt.MaxNLocator(num_ticks))

    # Make the ticks go into and out of the plot and make them larger
    ax.tick_params(axis='both', which='major', direction='inout', length=10, width=1)

    # Re-set the limits
    ax.set_ylim(min_val, max_val)
    ax.set_xlim(min_val, max_val)

    if cluster:
        points, frequencies = calculate_density_scatter(x, y, delta_factor=delta_factor, seed=seed)
        # revert the order of both to make the high frequency points appear on top:
        points = points[::-1]
        frequencies = frequencies[::-1]
        norm = plt.Normalize(vmin=min(frequencies), vmax=max(frequencies))

        if logscale:
            norm = colors.LogNorm(vmin=min(frequencies), vmax=max(frequencies))
        else:
            norm = plt.Normalize(vmin=min(frequencies), vmax=max(frequencies))
            
        sc = ax.scatter(points[:, 0], points[:, 1], c=frequencies, cmap=cmap, norm=norm, s=s, alpha=alpha, edgecolor='k', linewidths=0.5, **kwargs)

        # Create an axes divider for the colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        
        # Add the colorbar
        cbar = plt.colorbar(sc, cax=cax)
        cbar.set_label('Frequency')

    else:
        ax.scatter(x, y, alpha=alpha, s=s, **kwargs)

    return ax

# %%

FONTSIZE = 15
# FONT = 'Arial'

# plt.rc('font', family=FONT)
plt.rc('xtick', labelsize=FONTSIZE)
plt.rc('ytick', labelsize=FONTSIZE)
plt.rc('axes', labelsize=FONTSIZE, titlesize=FONTSIZE+2)
plt.rc('legend', fontsize=FONTSIZE, title_fontsize=FONTSIZE-2)


# top_k = 3

# Create subplots
fig, axs = plt.subplots(3, 2, figsize=(12, 18))

# Scatter plot of energies
axs[0, 0] = scatter_plot(axs[0, 0], energies_concat, ff_energies_concat, cluster=True, logscale=True)
axs[0, 0].set_xlabel('QM Energy')
axs[0, 0].set_ylabel(f'{ff} Energy')
axs[0, 0].set_title('Energy Comparison')

# Scatter plot of gradients
scatter_plot(axs[0, 1], gradients_concat, ff_gradients_concat, cluster=True, logscale=True, symmetric=True)
axs[0, 1].set_xlabel('QM Gradient')
axs[0, 1].set_ylabel(f'{ff} Gradient')
axs[0, 1].set_title('Gradient Comparison')

# Histograms
axs[1, 0].hist(rmse_energies, bins=30, alpha=0.7)
axs[1, 0].set_xlabel('RMSE Energies')
axs[1, 0].set_ylabel('Frequency')
axs[1, 0].set_title('Histogram of RMSE between Energies')


axs[1, 1].hist(rmse_gradients, bins=30, alpha=0.7)
axs[1, 1].set_xlabel('RMSE Gradients')
axs[1, 1].set_ylabel('Frequency')
axs[1, 1].set_title('Histogram of RMSE between Gradients')

axs[2, 0].hist(max_energy_qm, bins=30, alpha=0.7)
axs[2, 0].set_xlabel('Max Energy QM')
axs[2, 0].set_ylabel('Frequency')
axs[2, 0].set_title('Histogram of Max Energy QM')

axs[2, 1].hist(max_force_qm, bins=30, alpha=0.7)
axs[2, 1].set_xlabel('Max Gradient QM')
axs[2, 1].set_ylabel('Frequency')
axs[2, 1].set_title('Histogram of Max Gradient QM')



# Save and show plot
plt.tight_layout()
plt.savefig(f'analysis_{ff}.png')
plt.show()
# %%
