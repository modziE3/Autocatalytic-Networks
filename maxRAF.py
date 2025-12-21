"""
Author: Luke Burton
Date: 24/11/2025

Used to compute maxRAFs. Function phi(R: set[Reaction], F: set[str]) -> set[Reaction] uses 
Algorithm from section 3. Mathematical aspects of RAFs:
Huson D, Xavier JC, Steel M.
 2024 Self-generating autocatalytic networks:
 structural results, algorithms and their relevance
 to early biochemistry. J. R. Soc. Interface 21:
 20230732.

Examples from HusonLab: https://github.com/husonlab/catrenet/tree/master/examples
"""

import re
from typing import List, Set


class Reaction:
    def __init__(self, label: str, reactants: list[str], catalyst_sets: list[set[str]], products: list[str]):
        self.label = label
        self.reactants = reactants
        self.catalyst_sets = catalyst_sets
        self.products = products

    def is_satisfied(self, available_agents: set) -> bool:
        return all(agent in available_agents for agent in self.reactants)
    
    def needs_change(self, available_agents: set) -> bool:
        return any(agent not in available_agents for agent in self.products)
    
    def __str__(self):
        return f"{self.label}: {' + '.join(self.reactants)} " \
           f"[{', '.join('{' + ', '.join(s) + '}' for s in self.catalyst_sets)}] " \
           f"-> {' + '.join(self.products)}"
    
    def non_complex_str(self, empty_catalyst_placeholder: bool = False):
        catalysts = [{"NOT_CATALYSED"}] if empty_catalyst_placeholder and self.catalyst_sets == [] else self.catalyst_sets
        return f"{self.label}: {'+'.join(self.reactants)} " \
           f"[{', '.join(','.join(s) for s in catalysts)}] " \
           f"-> {'+'.join(self.products)}"
    
    def __repr__(self):
        return self.label
    
    def is_catalyzed(self, available_agents: set) -> bool:
        return any(all(u in available_agents for u in U) for U in self.catalyst_sets)
    
    def is_strictly_autocatalyzed(self, available_agents, food_set: set) -> bool:
        return any(all(u in available_agents for u in U) and not U.issubset(food_set) for U in self.catalyst_sets)
    
    def rho(self):
        return set(self.reactants)

    def pi(self):
        return set(self.products)


def reaction_str_to_class(reaction_str: str) -> Reaction:
    pattern = re.compile(
        r'^(?P<head>\w+):\s+'
        r'(?P<reactants>[\w\s\+]+)\s+'
        r'(?P<catalysts>\[.*\])\s*->\s*'
        r'(?P<products>[\w\s\+]+)$'
    )

    m = pattern.fullmatch(reaction_str)
    if not m:
        raise ValueError(f"Invalid reaction string: {reaction_str}")
    head = m.group("head")
    reactants = [x.strip() for x in m.group("reactants").split("+")]
    products = [x.strip() for x in m.group("products").split("+")]
    catalysts_raw = m.group("catalysts").strip("[]").strip()
    catalysts: List[Set[str]] = []
    if catalysts_raw == "":  
        catalysts = []
    else:
        parts = [p.strip() for p in re.split(r',(?![^{]*\})', catalysts_raw)]
        for part in parts:
            if part.startswith("{") and part.endswith("}"):
                inner = part.strip("{}").strip()
                if inner == "":
                    catalysts.append(set())
                else:
                    catalysts.append(set(x.strip() for x in inner.split(",") if x.strip()))
            else:
                catalysts.append({part})
    return Reaction(head, reactants, catalysts, products)

def closure(reactions: set[Reaction], food_set: set[str]) -> set[str]:
    """Computes the closure of a set of reactions under a given food set. 
       Returns the computed closure set of reactions. 
    """
    availabe_agents = set()
    for f in food_set: availabe_agents.add(f)
    changed = True
    while changed:
        changed = False
        for r in reactions:
            if r.is_satisfied(availabe_agents) and r.needs_change(availabe_agents):
                changed = True
                availabe_agents.update(r.products)
    return availabe_agents
    
def phi(R: set[Reaction], F: set[str]) -> set[Reaction]:
    Rk = set(R)
    while len(Rk) > 0:
        Rk_plus_one = set()
        closure_Rk = closure(Rk, F)
        for r in Rk:
            if r.rho().issubset(closure_Rk) and r.is_catalyzed(closure_Rk):
                Rk_plus_one.add(r)
        if Rk == Rk_plus_one: break
        Rk = Rk.intersection(Rk_plus_one)
    return Rk

def strictly_autocatalytic_RAF(R: set[Reaction], F: set[str]) -> set[Reaction]:
    Rk = set(R)
    while len(Rk) > 0:
        Rk_plus_one = set()
        closure_Rk = closure(Rk, F)
        for r in Rk:
            if r.rho().issubset(closure_Rk) and r.is_strictly_autocatalyzed(closure_Rk, F):
                Rk_plus_one.add(r)
        if Rk == Rk_plus_one: break
        Rk = Rk.intersection(Rk_plus_one)
    return Rk

def R_Q_poly(R: set[Reaction], F: set[str]) -> set[Reaction]:
    if phi(R,F)==set(): return set()
    Rk = R
    changed = True
    while changed:
        changed = False
        discard_set = set()
        for reaction in Rk:
            Rk_changed = {r for r in Rk if r != reaction}
            if phi(Rk_changed, F) != set():
                changed = True
                discard_set.add(reaction)
        Rk = {r for r in Rk if r not in discard_set}
        if len(Rk) == 0:
            break
    return Rk

def R_Q_poly2(R: set[Reaction], F: set[str]) -> set[Reaction]:
    if phi(R,F)==set(): return set()
    return {r for r in R if set() == phi(F=F, R={
        r_prime for r_prime in R if r != r_prime
    })}

def R_Q_exp(R: set[Reaction], F: set[str]) -> set[Reaction]:
    try:
        return set.intersection(*all_rafs(R, F))
    except:
        return set()

def all_rafs(R: set[Reaction], F: set[str]) -> set[set[Reaction]]:

    def all_sub_rafs(R: set[Reaction], F: set[str]) -> set[set[Reaction]]:
        rafs = {frozenset(phi(R, F))}
        for reaction in R:
            R_changed = frozenset({r for r in R if r != reaction})
            rafs.update(all_sub_rafs(frozenset(phi(R_changed, F)), F))
        try:
            rafs.remove(frozenset())
        except KeyError:
            pass
        return rafs

    return [set(raf) for raf in all_sub_rafs(R, F)]
        
reaction_str_set = {
    'r1: a + a [{c,d},e] -> c',
    'r2: b + c [{}] -> d',
    'r3: b + b [] -> e',
    'r4: a + e [{a}] -> b',
    'r5: c+d [{d}]->g+g'
}

reaction_str_set2 = {
    'r1: a+b [{a,d}] -> e',      
    'r2: b+c [{a,b},{e}] -> d',   
    'r3: d [{a,b}] -> c'        
}

# maxCAF none, maxRAF size 3, max-pRAF size 6
example_0 = {
    "food_set": {'a','b','c'},
    "reaction_set": {
        reaction_str_to_class(r_str) for r_str in {
            'r1: a+b [{e},{f}] -> d',
            'r2: b+c [{f}] -> e',
            'r3: e+c [{d}] -> b+f',
            'r4: c+h [{g}] -> g',
            'r5: i [{h}] -> h',
            'r6: h [{d}] -> i',
        }
    }
}

# maxCAF none, maxRAF size 6, max-pRAF size 6
example_1 = {
    "food_set": {"a1","a2","a3","ap1","ap2","ap3","b1","b2","b3","bp1","bp2","bp3"},
    "reaction_set": {
        reaction_str_to_class(r_str) for r_str in {
            "r1: a1+b1 [{c3}] -> c1",
            "r2: a2+b2 [{c1}] -> c2",
            "r3: a3+b3 [{c2}] -> c3",
            "rp1: ap1+bp1 [{c3}] -> c1",
            "rp2: ap2+bp2 [{c1}] -> c2",
            "rp3: ap3+bp3 [{c2}] -> c3",
        }
    }
}

# Has a maxRAF of size 4 (which is an iRAF),  a pseudoRAF of size 7 (everything) and no CAF.
example_9 = {
    "food_set": {'f1','f2','f3','f4','f5','f6',},
    "reaction_set": {
        reaction_str_to_class(r_str) for r_str in {
            'r1: f1 + f2 [{g}] -> a',
            'r2: a + f3 [{c}] -> b',
            'r3: b + f4 [{a}] -> c',
            'r4: c + f5 [{b}] -> a + g',
            'r5: b + c [{h}] -> d',
            'r6: e + f6 [{b}] -> h',
            'r7: f5 + h [{c}] -> e',
        }
    }
}

example_custom_0 = {
    "food_set": {'0','1','00','01','10','11'},
    "reaction_set": {
        reaction_str_to_class(r_str) for r_str in {
            'r1: 10 + 0 [{01100}] -> 100',
            'r2: 01 + 100 [{0}] -> 01100',
            'r3: 10 + 1 [{0}] -> 101',
            'r4: 11 + 10 [{101}] -> 1110',
            'r5: 1110 + 0 [{101}] -> 11100',
        }
    }
}

example_custom_1 = {
    "food_set": {'f'},
    "reaction_set": {
        reaction_str_to_class(r_str) for r_str in {
            'r1: f [{c2}] -> c1',
            'r2: f [{c1}, {c3}] -> c2',
            'r3: f [{c2}] -> c3',
        }
    }
}

example_custom_2 = {
    "food_set": {'f'},
    "reaction_set": {
        reaction_str_to_class(r_str) for r_str in {
            'r1: f [{c3}] -> c1',
            'r2: f [{c1}, {c3}] -> c2',
            'r3: f [{c2}, {c4}] -> c3',
            'r4: f [{c2}, {c3}] -> c4',
        }
    }
}


example_custom_3 = {
    "food_set": {'f1', 'f2'},
    "reaction_set": {
        reaction_str_to_class(r_str) for r_str in {
            'r11: f1+f2 [{UA},{UC},{UG},{UU}] -> AA',
            'r12: f1+f2 [{GA},{GC},{GG},{GU}] -> AC',
            'r13: f1+f2 [{CA},{CC},{CG},{CU}] -> AG',
            'r14: f1+f2 [{AA},{AC},{AG},{AU}] -> AU',
            'r21: f1+f2 [{UA},{UC},{UG},{UU}] -> CA',
            'r22: f1+f2 [{GA},{GC},{GG},{GU}] -> CC',
            'r23: f1+f2 [{CA},{CC},{CG},{CU}] -> CG',
            'r24: f1+f2 [{AA},{AC},{AG},{AU}] -> CU',
            'r31: f1+f2 [{UA},{UC},{UG},{UU}] -> GA',
            'r32: f1+f2 [{GA},{GC},{GG},{GU}] -> GC',
            'r33: f1+f2 [{CA},{CC},{CG},{CU}] -> GG',
            'r34: f1+f2 [{AA},{AC},{AG},{AU}] -> GU',
            'r41: f1+f2 [{UA},{UC},{UG},{UU}] -> UA',
            'r42: f1+f2 [{GA},{GC},{GG},{GU}] -> UC',
            'r43: f1+f2 [{CA},{CC},{CG},{CU}] -> UG',
            'r44: f1+f2 [{AA},{AC},{AG},{AU}] -> UU',
        }
    }
}



def print_maxRAF(example):
    reactions = example["reaction_set"]
    food_set = example["food_set"]
    print(f"maxRAF(R) = {phi(reactions, food_set)}\n")

if __name__ == "__main__":
    # print_maxRAF(example_0)
    # print_maxRAF(example_1)
    # print_maxRAF(example_9)

    # reactions = {reaction_str_to_class(r_str) for r_str in reaction_str_set}
    # food_set = {'a', 'b'}
    # print(f"maxRAF(R) = {phi(reactions, food_set)}\n")
    # for raf in all_rafs(reactions, food_set):
    #     print(raf)

    example = example_custom_2

    print("All RAFs:")
    for raf in all_rafs(example["reaction_set"], example["food_set"]):
        print(raf)

    print("\nPersistent reactions:")
    temp = R_Q_exp(example["reaction_set"], example["food_set"])
    print(f"exp: {temp}")
    temp = R_Q_poly2(example["reaction_set"], example["food_set"])
    print(f"poly2: {temp}")


    # reactions = {reaction_str_to_class(r_str) for r_str in reaction_str_set2}
    # food_set = {'a', 'b', 'c'}
    # print(f"strictly autocatalytic RAF = {strictly_autocatalytic_RAF(reactions, food_set)}\n")