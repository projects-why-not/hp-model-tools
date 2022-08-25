from .backend import Plotter
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_3d_fit(hp_sequence, core, monomer_coords):

    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(projection="3d")
    Plotter.plot_fit(ax, core, hp_sequence, monomer_coords)

    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    ax.set_zlim(-10,10)

    plt.show()


def plot_2d_fit(hp_sequence, core, monomer_coords):
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    Plotter.plot_fit(ax, core, hp_sequence, monomer_coords)

    ax.axis("equal")

    plt.show()


def plot_fit(hp_sequence, core, monomer_coords, fit_3d):
    if fit_3d:
        plot_3d_fit(hp_sequence, core, monomer_coords)
    else:
        plot_2d_fit(hp_sequence, core, monomer_coords)
