import math
import matplotlib.pyplot as plt
import numpy as np
import igraph as ig
import binary_polymer_model as bpm

def get_digraph_cycle_probability_bounds(f, sample=100):
    upper_bound = sum(((f**(i+2)) / (i+2)) for i in range(sample))
    lower_bound = 1 - math.exp(-(f**2)/2)
    return max(lower_bound, 0), min(upper_bound, 1)

def get_digraph_cycle_probability_bounds_span(f_span, sample=100):
    lower_span = []
    upper_span = []
    for f in f_span:
        lower, upper = get_digraph_cycle_probability_bounds(f, sample)
        lower_span.append(lower)
        upper_span.append(upper)
    return lower_span, upper_span

def plot_digraph_cycle_probability(min_f = 0, max_f = 3.5, num = 100):
    f_span = np.linspace(min_f, max_f, num)
    lower_p, upper_p = get_digraph_cycle_probability_bounds_span(f_span, num)
    plt.plot(f_span, lower_p)
    plt.plot(f_span, upper_p)
    plt.grid(True)
    plt.show()

def crs_digraph_has_directed_cycle(crs: bpm.CRS) -> bool:
    reactions = list(crs.reactions)
    id_of = {r: i for i, r in enumerate(reactions)}
    pi_of = {r: r.pi() for r in reactions}

    edges = []

    for r_from in reactions:
        produced = pi_of[r_from]
        i = id_of[r_from]

        for r_dest in reactions:
            if any(cat.issubset(produced) for cat in r_dest.catalyst_sets):
                j = id_of[r_dest]
                edges.append((i, j))

    g = ig.Graph(n=len(reactions), edges=edges, directed=True)

    return not g.is_dag()

def crs_digraph_has_RAF(crs: bpm.CRS) -> bool:
    from maxRAF import phi
    return phi(crs.reactions, crs.food_set) != set()

if __name__ == "__main__":
    # plot_digraph_cycle_probability(0, 3.5, 3000)
    # from maxRAF import reaction_str_to_class
    # crs = bpm.CRS(
    #     {reaction_str_to_class(r_str) for r_str in {
    #         'r1: f [{c3}] -> c1',
    #         'r2: f [{c1}, {c3}] -> c2',
    #         'r3: f [{c2}, {c4}] -> c3',
    #         'r4: f [{c2}, {c3}] -> c4',
    #     }},
    #     {'f'},
    # )
    import crs_file_read_write

    # gen = bpm.BinaryCRSGenerator()
    # gen.generate_reactions(3, 1, 2)

    # for _ in range(10000):
    #     gen.catalyze_reactions(0.01, False)
    #     if crs_digraph_has_directed_cycle(gen.CRS) != crs_digraph_has_RAF(gen.CRS):
    #         print("found!")
    #         crs_file_read_write.export_crs(gen.CRS, "degen_case.crs", True)
    #         break

    crs = crs_file_read_write.import_crs("degen_case.crs")
    print(crs_digraph_has_directed_cycle(crs))
    print(crs_digraph_has_RAF(crs))
     