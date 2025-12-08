from maxRAF import *
import random


def generate_reactions(n, p, t=2, l=2):
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
    for reaction in reactions:
        for element in elements:
            if random.random() <= p:
                reaction.catalyst_sets.append({element})
    return reactions
                    

def contains_reaction(reaction_set, reaction):
    for r in reaction_set:
        if r.reactants == reaction.reactants and r.products == reaction.products:
            return True
    return False


if __name__ == "__main__":
    for r in generate_reactions(8,0.001):
        print(r)