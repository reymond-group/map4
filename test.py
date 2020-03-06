from rdkit import Chem
import tmap as tm
from map4 import MAP4Calculator

dim = 1024

MAP4 = MAP4Calculator(dimensions=dim)
ENC = tm.Minhash(dim)

smiles_a = 'c1ccccc1'
mol_a = Chem.MolFromSmiles(smiles_a)
map4_a = MAP4.calculate(mol_a)


smiles_b = 'c1cccc(N)c1'
mol_b = Chem.MolFromSmiles(smiles_b)
map4_b = MAP4.calculate(mol_b)

# or use parallelized version:
fps = MAP4.calculate_many([mol_a, mol_b])


print(ENC.get_distance(map4_a, map4_b))

print(ENC.get_distance(fps[0], fps[1]))