from maxRAF import Reaction, reaction_str_to_class
from binary_polymer_model import CRS
import os


def export_crs(crs: CRS, write_filename: str, overwrite_existing: bool = False, empty_catalyst_placeholder: bool = False) -> None:
    if (not os.path.exists(write_filename)) or overwrite_existing:
        with open(write_filename, 'w') as f:
            f.write("# " + write_filename + "\n\n")
            f.write("Food: " + str.join(",", crs.food_set) + "\n\n")
            for reaction in crs.reactions:
                f.write(reaction.non_complex_str(empty_catalyst_placeholder) + "\n")
    else:
        print(f"The file {write_filename} already exists. Not writing.")

def import_crs(read_filename: str) -> CRS:
    try:
        with open(read_filename, 'r') as file:
            content = file.read()

            reactions = set()
            
            for line in content.splitlines():
                if line.isspace() or line == "" or line[0] == '#':
                    continue

                if line[:6] == "Food: ":
                    food_set = set(line[6:].split(", "))
                    continue

                reaction = reaction_str_to_class(line)
                if reaction.catalyst_sets == [{"NOT_CATALYSED"}]: reaction.catalyst_sets == []
                reactions.add(reaction)

        return CRS(reactions, food_set)

    except FileNotFoundError:
        print("The file was not found")


if __name__ == "__main__":
    import binary_polymer_model as bpm
    gen = bpm.BinaryCRSGenerator()
    gen.generate_reactions(n=3)
    gen.catalyze_reactions_level_of_catalysis(mean_catalysts=2)
    export_crs(gen.CRS, "my_reactions_1_will_overwrite.crs", True, False)
    crs2 = import_crs("my_reactions_1_will_overwrite.crs")
    export_crs(crs2, "my_reactions_2_will_overwrite.crs", True, True)

    from maxRAF import phi
    print(len(raf:=phi(gen.CRS.reactions, gen.CRS.food_set)))
    print(raf)
