from maxRAF import Reaction
from binary_polymer_model import CRS, BinaryCRSGenerator

def CAF_existence(crs: CRS) -> bool:
    """ Determines if a given CRS contains a Constructive Autocatalytic 
        Food-generated set of reactions. This is done in linear time with |R|.
    """
    for reaction in crs.reactions:
        food_constructed = all(reactant in crs.food_set for reactant in reaction.rho())
        food_catalysed = any(catalyst.issubset(crs.food_set) for catalyst in reaction.catalyst_sets)
        if food_constructed and food_catalysed:
            return True
    return False


if __name__ == "__main__":
    r1 = Reaction(
        label= 'r1',
        reactants= ['a','b'], 
        catalyst_sets= [{'c'}, {'d', 'e'}], 
        products= ['f']
    )
    r2 = Reaction(
        label= 'r2',
        reactants= ['x'], 
        catalyst_sets= [{'y'}], 
        products= ['z']
    )
    crs1 = CRS({r1, r2}, {'a', 'b', 'c'})
    crs2 = CRS({r1, r2}, {'a', 'b', 'e'})
    print(CAF_existence(crs1)) # true
    print(CAF_existence(crs2)) # false