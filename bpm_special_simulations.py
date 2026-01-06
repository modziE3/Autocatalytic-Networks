from maxRAF import *
import random
import numpy as np
import matplotlib.pyplot as plt
from typing import NamedTuple
from binary_polymer_model import BinaryCRSGenerator, number_of_reactions, plot_n_range_varied_mean_catalysts
from special_functions import CAF_existence

def get_RAF_size_span_from_mc_range(n, mc_span, sample_size, t=2, l=2):
    RAF_size_span = []

    generator = BinaryCRSGenerator()
    generator.generate_reactions(n, t, l)

    for i in range(len(mc_span)):
        mc = mc_span[i]
        max_raf_running_size_num = 0
        for j in range(sample_size):
            print(f"n = {n}: Processing mc index={i} out of {len(mc_span)}: {j/sample_size * 100 :.0f}% complete", end='\r')
            generator.catalyze_reactions_level_of_catalysis(mc)
            max_raf_running_size_num += len(phi(generator.CRS.reactions, generator.CRS.food_set))
        RAF_size_span.append(max_raf_running_size_num / sample_size)
    return mc_span, RAF_size_span

def get_CAF_probability_from_mc_range(n, mc_span, sample_size, t=2, l=2):
    RAF_size_span = []

    generator = BinaryCRSGenerator()
    generator.generate_reactions(n, t, l)

    for i in range(len(mc_span)):
        mc = mc_span[i]
        caf_count = 0
        for j in range(sample_size):
            print(f"n = {n}: Processing mc index={i} out of {len(mc_span)}: {j/sample_size * 100 :.0f}% complete", end='\r')
            generator.catalyze_reactions_level_of_catalysis(mc)
            if CAF_existence(generator.CRS): caf_count += 1
        RAF_size_span.append(caf_count / sample_size)
    return mc_span, RAF_size_span

def plot_n_range_special(n_range, func, args: list, name, x_label, y_label, save_to_file = False):
    for n in n_range:
        x, y = func(n, *args)
        plt.plot(x, y, label = f"n = {n}")
    plt.grid(True)
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if save_to_file: plt.savefig(f"{name}.png")
    else: plt.show()

def plot_n_range_RAF_size(n_range, args: list, name, x_label, y_label, save_to_file = False):
    for n in n_range:
        x, y = get_RAF_size_span_from_mc_range(n, *args)
        max_y = number_of_reactions(n)

        (line,) = plt.plot(x, y, label=f"n = {n}")

        plt.axhline(
            y=max_y,
            linestyle=":",
            color=line.get_color(),
            label = f"|R| for n={n}"
        )
    plt.grid(True)
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if save_to_file: plt.savefig(f"({name})(n_range={n_range})(sample={args[1]})(number_of_points={len(args[0])}).png")
    else: plt.show()


if __name__ == "__main__":
    # mc_start = 1.5
    # mc_end = 20
    # number_of_points = 50
    # span = np.linspace(mc_start, mc_end, number_of_points)
    # sample_size = 20
    # plot_n_range_RAF_size(
    #     [4,5,6],
    #     [span, sample_size],
    #     "average_raf_size",
    #     "Level of Catalysis",
    #     "Average RAF Size",
    #     True
    # )

    # mc_start = 0
    # mc_end = 3.5
    # number_of_points = 50
    # span = np.linspace(mc_start, mc_end, number_of_points)
    # n_range = [4,5,6,7]
    # sample_size = 500
    # plot_n_range_varied_mean_catalysts(n_range, span, sample_size, allow_food_catalyst=False)

    mc_start = 0
    mc_end = 4.5
    number_of_points = 50
    span = np.linspace(mc_start, mc_end, number_of_points)
    n_range = [4,5,6,7]
    sample_size = 500
    plot_n_range_special(
        n_range,
        get_CAF_probability_from_mc_range,
        [span, sample_size],
        f"(caf_probalility)(n_range={n_range})(sample={sample_size})(number_of_points={number_of_points})",
        "Level of Catalysis",
        "Probability of CAF",
        save_to_file=True
    )