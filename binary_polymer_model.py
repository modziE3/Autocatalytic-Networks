from maxRAF import *
import random
import numpy as np
import matplotlib.pyplot as plt
from typing import NamedTuple


class CRS(NamedTuple):
    reactions: set[Reaction]
    food_set: set

class BinaryCRSGenerator:
    def __init__(self):
        self.CRS = None
        self.elements = []

    def generate_reactions(self, n, t=2, l=2):
        elements = set(['0','1','2','3','4','5','6','7','8','9'][:l])
        new_elements = set()
        reactions = set()
        changed = True
        while changed:
            changed = False
            for a in elements:
                for b in elements:
                    c = a+b
                    if len(c) <= n and not contains_reaction(reactions, Reaction('', [a, b], [], [a+b])):
                        changed = True
                        new_elements.add(c) # new element
                        reactions.add(Reaction(f"r{len(reactions)}", [a, b], [], [a+b])) #concat
                        reactions.add(Reaction(f"r{len(reactions)}", [a+b], [], [a, b])) #cutting
            elements.update(new_elements)
            new_elements = set()
        self.elements = elements
        self.CRS = CRS(reactions, {element for element in elements if len(element) <= t})

    def catalyze_reactions(self, p, allow_food_catalyst = True):
        for reaction in self.CRS.reactions:
            reaction.catalyst_sets = []
            for element in self.elements:
                if not allow_food_catalyst and element in self.CRS.food_set: continue
                if random.random() <= p:
                    reaction.catalyst_sets.append({element})

    def catalyze_reactions_mean_catalysts_per_reaction(self, mean_catalysts, allow_food_catalyst = True):
        self.catalyze_reactions(mean_catalysts / len(self.elements), allow_food_catalyst)

    def catalyze_reactions_level_of_catalysis(self, mean_catalysts, allow_food_catalyst = True):
        self.catalyze_reactions(mean_catalysts / len(self.CRS.reactions), allow_food_catalyst)


def contains_reaction(reaction_set, reaction):
    for r in reaction_set:
        if r.reactants == reaction.reactants and r.products == reaction.products:
            return True
    return False

def number_of_reactions(n, l=2):
    return 2*sum((l**k) * (k-1) for k in range(n+1)[1:])

def get_probability_span_from_mc_range(n, mc_span, sample_size, t=2, l=2, allow_food_catalyst = True):
    probablity_span = []

    generator = BinaryCRSGenerator()
    generator.generate_reactions(n, t, l)

    for i in range(len(mc_span)):
        mc = mc_span[i]
        max_raf_count = 0
        for j in range(sample_size):
            print(f"n = {n}: Processing mc index={i} out of {len(mc_span)}: {j/sample_size * 100 :.0f}% complete", end='\r')
            generator.catalyze_reactions_level_of_catalysis(mc, allow_food_catalyst)
            if phi(generator.CRS.reactions, generator.CRS.food_set) != set(): max_raf_count += 1
        probablity_span.append(max_raf_count / sample_size)
    return probablity_span

def plot_varied_mean_catalysts(n, mc_span, sample_size, t=2, l=2):
    plt.plot(mc_span, get_probability_span_from_mc_range(n, mc_span, sample_size, t, l))
    plt.grid(True)
    plt.savefig(f"(n={n})(sample_size={sample_size})(number_of_points={len(mc_span)}).png")
    # plt.show() #optional show graph

def plot_n_range_varied_mean_catalysts(n_range, mc_span, sample_size, t=2, l=2, allow_food_catalyst = True):
    for n in n_range:
        plt.plot(mc_span, get_probability_span_from_mc_range(n, mc_span, sample_size, t, l, allow_food_catalyst), label = f"n = {n}")
    plt.grid(True)
    plt.legend()
    plt.xlabel("Level of Catalysis")
    plt.ylabel("Probability of a RAF")
    plt.savefig(f"(n_range={n_range})(sample_size={sample_size})(number_of_points={len(mc_span)})(food_catalysts={allow_food_catalyst}).png")
    # plt.show() #optional show graph

if __name__ == "__main__":
    mc_start = 0
    mc_end = 3.5
    number_of_points = 50
    span = np.linspace(mc_start, mc_end, number_of_points)

    # n = 8
    # l = 2
    # t = 2
    # sample = 1
    # plot_varied_mean_catalysts(n, span, sample, t, l)

    n_range = [4, 5, 6, 7, 8]
    l = 2
    t = 2
    sample = 500
    plot_n_range_varied_mean_catalysts(n_range, span, sample, t, l)