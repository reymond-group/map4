#!/usr/bin/env python

import argparse
import itertools
from collections import defaultdict

import tmap as tm
from mhfp.encoder import MHFPEncoder
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdmolops import GetDistanceMatrix


def to_smiles(mol):
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)


class MAP4Calculator:

    def __init__(self, dimensions=1024, radius=2, is_counted=False, is_folded=False, return_strings=False):
        """
        MAP4 calculator class
        """
        self.dimensions = dimensions
        self.radius = radius
        self.is_counted = is_counted
        self.is_folded = is_folded
        self.return_strings = return_strings

        if self.is_folded:
            self.encoder = MHFPEncoder(dimensions)
        else:
            self.encoder = tm.Minhash(dimensions)

    def calculate(self, mol):
        """Calculates the atom pair minhashed fingerprint

        Arguments:
            mol -- rdkit mol object

        Returns:
            tmap VectorUint -- minhashed fingerprint
        """
        
        atom_env_pairs = self._calculate(mol)
        if self.is_folded:
            return self._fold(atom_env_pairs)
        elif self.return_strings:
            return atom_env_pairs
        return self.encoder.from_string_array(atom_env_pairs)

    def calculate_many(self, mols):
        """ Calculates the atom pair minhashed fingerprint

        Arguments:
            mols -- list of mols

        Returns:
            list of tmap VectorUint -- minhashed fingerprints list
        """

        atom_env_pairs_list = [self._calculate(mol) for mol in mols]
        if self.is_folded:
            return [self._fold(pairs) for pairs in atom_env_pairs_list]
        elif self.return_strings:
            return atom_env_pairs_list
        return self.encoder.batch_from_string_array(atom_env_pairs_list)

    def _calculate(self, mol):
        return self._all_pairs(mol, self._get_atom_envs(mol))

    def _fold(self, pairs):
        fp_hash = self.encoder.hash(set(pairs))
        return self.encoder.fold(fp_hash, self.dimensions)

    def _get_atom_envs(self, mol):
        atoms_env = {}
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            for radius in range(1, self.radius + 1):
                if idx not in atoms_env:
                    atoms_env[idx] = []
                atoms_env[idx].append(MAP4Calculator._find_env(mol, idx, radius))
        return atoms_env

    @classmethod
    def _find_env(cls, mol, idx, radius):
        env = rdmolops.FindAtomEnvironmentOfRadiusN(mol, radius, idx)
        atom_map = {}

        submol = Chem.PathToSubmol(mol, env, atomMap=atom_map)
        if idx in atom_map:
            smiles = Chem.MolToSmiles(submol, rootedAtAtom=atom_map[idx], canonical=True, isomericSmiles=False)
            return smiles
        return ''

    def _all_pairs(self, mol, atoms_env):
        atom_pairs = []
        distance_matrix = GetDistanceMatrix(mol)
        num_atoms = mol.GetNumAtoms()
        shingle_dict = defaultdict(int)
        for idx1, idx2 in itertools.combinations(range(num_atoms), 2):
            dist = str(int(distance_matrix[idx1][idx2]))

            for i in range(self.radius):
                env_a = atoms_env[idx1][i]
                env_b = atoms_env[idx2][i]

                ordered = sorted([env_a, env_b])

                shingle = '{}|{}|{}'.format(ordered[0], dist, ordered[1])

                if self.is_counted:
                    shingle_dict[shingle] += 1
                    shingle += '|' + str(shingle_dict[shingle])

                atom_pairs.append(shingle.encode('utf-8'))
        return list(set(atom_pairs))


def main():
    args = parse_args()

    def _parse_line(line):
        line = line.strip()
        fields = line.split(args.delimiter)
        mol = Chem.MolFromSmiles(fields[0])
        if mol:
            if args.clean_mols:
                mol = sorted(Chem.GetMolFrags(mol, asMols=True),
                             key=lambda mol: mol.GetNumHeavyAtoms(), reverse=True)[0]
                mol = Chem.MolFromSmiles(to_smiles(mol))
            return (line, mol)
        else:
            return None

    calculator = MAP4Calculator(args.dimensions, args.radius, args.is_counted, args.is_folded)


    def process(batch, output_file):
        parsed_lines = [_parse_line(line) for line in batch]
        parsed_lines = [_tuple for _tuple in parsed_lines if _tuple is not None] #remove all lines with unreadable mols
        lines, mols = zip(*parsed_lines)
        fingerprints = calculator.calculate_many(mols)
        for line, mol, fingerprint in zip(lines, mols, fingerprints):
            if len(fingerprint):
                fp_str = args.fp_delimiter.join(str(v) for v in fingerprint)
                output_file.write("{}{}{}{}{}\n".format(line, args.delimiter, to_smiles(mol), args.delimiter, fp_str))

    with open(args.input_path, "r") as input_file:
        with open(args.output_path, "w+") as output_file:
            batch = []
            for line in input_file:
                batch.append(line)
                if len(batch)>=args.batch_size:
                    process(batch, output_file)
                    batch=[]
            process(batch, output_file)

def parse_args():
    parser = argparse.ArgumentParser(description="MAP4 calculator")
    parser.add_argument("--input-path", "-i", help="", type=str, required=True)
    parser.add_argument("--output-path", "-o", help="", type=str, required=True)
    parser.add_argument("--dimensions", "-d", help="Number of dimensions of the MinHashed fingerprint [DEFAULT: 1024]",
                        type=int, default=1024, choices = [128, 512, 1024, 2048])
    parser.add_argument("--radius", "-r", help="Radius of the fingerprint [DEFAULT: 2]",
                        type=int, default=2)
    parser.add_argument("--is-counted", help="The fingerprint stores all shingles.",
                        action="store_true", default=False)
    parser.add_argument("--is-folded", help="The fingerprint is folded with modulo (instead of MinHash).",
                        action="store_true", default=False)
    parser.add_argument("--clean-mols", help="Molecules will be canonicalized, cleaned, and chirality information will be removed, \
    NECESSARY FOR FINGERPRINT CONSISTENCY ACROSS DIFFERENT SMILES INPUT [DEFAULT: True].",
                        type=lambda x: (str(x).lower() == "true"), default="True", metavar = "True/False")
    parser.add_argument("--delimiter", help="Delimiter used for both the input and output files [DEFAULT: \\t]", type=str, default="\t")
    parser.add_argument("--fp-delimiter", help="Delimiter used between the numbers in the fingerprint output [DEFAULT: ;]", type=str, default=";")
    parser.add_argument("--batch-size", "-b", help="Numbers of molecules to process in a batch [DEFAULT: 500]", type=int, default=500)
    return parser.parse_args()


if __name__ == "__main__":
    main()
