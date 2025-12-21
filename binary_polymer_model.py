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

    def catalyze_reactions(self, p):
        for reaction in self.CRS.reactions:
            reaction.catalyst_sets = []
            for element in self.elements:
                if random.random() <= p:
                    reaction.catalyst_sets.append({element})

    def catalyze_reactions_mean_catalysts_per_reaction(self, mean_catalysts):
        self.catalyze_reactions(mean_catalysts / len(self.elements))

    def catalyze_reactions_level_of_catalysis(self, mean_catalysts):
        self.catalyze_reactions(mean_catalysts / len(self.CRS.reactions))


def contains_reaction(reaction_set, reaction):
    for r in reaction_set:
        if r.reactants == reaction.reactants and r.products == reaction.products:
            return True
    return False

def number_of_reactions(n, l=2):
    return 2*sum((l**k) * (k-1) for k in range(n+1)[1:])

def plot_varied_mean_catalysts(n, mc_span, sample_size, t=2, l=2):
    probablity_span = []

    generator = BinaryCRSGenerator()
    generator.generate_reactions(n, t, l)

    for i in range(len(mc_span)):
        mc = mc_span[i]
        max_raf_count = 0
        for j in range(sample_size):
            print(f"Processing mc index={i} out of {len(mc_span)}: {j/sample * 100 :.0f}% complete", end='\r')
            generator.catalyze_reactions_level_of_catalysis(mc)
            if phi(generator.CRS.reactions, generator.CRS.food_set) != set(): max_raf_count += 1
        probablity_span.append(max_raf_count / sample_size)

    plt.plot(mc_span, probablity_span)
    plt.grid(True)
    # plt.show()

    plt.savefig(f"(n={n})(sample_size={sample_size})(number_of_points={len(mc_span)}).png")


if __name__ == "__main__":
    mc_start = 1
    mc_end = 3
    number_of_points = 50
    span = np.linspace(mc_start, mc_end, number_of_points)

    n = 4
    l = 2
    t = 2
    sample = 2000
    plot_varied_mean_catalysts(n, span, sample, t, l)