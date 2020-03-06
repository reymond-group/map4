from rdkit import Chem
from rdkit.Chem import rdchem, rdmolfiles, rdmolops

# pop_size, mut_rate, gen_gap, query, sim treshold


def attach_capping(mol1, mol2):
    """it is connecting all Nterminals with the desired capping

    Arguments:
        mol1 {rdKit mol object} -- first molecule to be connected
        mol2 {rdKit mol object} -- second molecule to be connected - chosen N-capping

    Returns:
        rdKit mol object -- mol1 updated (connected with mol2, one or more)
    """

    count = 0

    # detects all the N terminals in mol1
    for atom in mol1.GetAtoms():
        atom.SetProp('Cterm', 'False')
        if atom.GetSmarts() == '[N:2]' or atom.GetSmarts() == '[NH2:2]' or atom.GetSmarts() == '[NH:2]':
            count += 1
            atom.SetProp('Nterm', 'True')
        else:
            atom.SetProp('Nterm', 'False')

    # detects all the C terminals in mol2 (it should be one)
    for atom in mol2.GetAtoms():
        atom.SetProp('Nterm', 'False')
        if atom.GetSmarts() == '[C:1]' or atom.GetSmarts() == '[CH:1]':
            atom.SetProp('Cterm', 'True')
        else:
            atom.SetProp('Cterm', 'False')

    # mol2 is addes to all the N terminal of mol1
    for i in range(count):
        combo = rdmolops.CombineMols(mol1, mol2)
        Nterm = []
        Cterm = []

        # saves in two different lists the index of the atoms which has to be connected
        for atom in combo.GetAtoms():
            if atom.GetProp('Nterm') == 'True':
                Nterm.append(atom.GetIdx())
            if atom.GetProp('Cterm') == 'True':
                Cterm.append(atom.GetIdx())

        # creates the amide bond
        edcombo = rdchem.EditableMol(combo)
        edcombo.AddBond(Nterm[0], Cterm[0], order=Chem.rdchem.BondType.SINGLE)
        clippedMol = edcombo.GetMol()

        # removes tags and lables form the atoms which reacted
        clippedMol.GetAtomWithIdx(Nterm[0]).SetProp('Nterm', 'False')
        clippedMol.GetAtomWithIdx(Cterm[0]).SetProp('Cterm', 'False')
        clippedMol.GetAtomWithIdx(Nterm[0]).SetAtomMapNum(0)
        clippedMol.GetAtomWithIdx(Cterm[0]).SetAtomMapNum(0)
        # uptades the 'core' molecule
        mol1 = clippedMol

    return mol1


class PDGA:

    # list of possible aminoacids
    AA = ['R', 'H', 'K', 'E', 'S', 'T', 'N', 'Q', 'G', 'P', 'A', 'V', 'I', 'L', 'F', 'Y', 'W', 'C', 'D', 'M', 'Z', 'O',
          '!', '?', '=', '%', '$', '@', '#']
    # list of possible branching units (1=Dap, 2=Dab, 3=Orn, 4=Lys)
    B = ['1', '2', '3', '4']
    # list of possible C-terminals
    CT = ['+']
    # list of possible N-capping
    NT = ['&']

    # variables for SMILES generation
    B_SMILES = {'1': '[N:2][C@@H](C[N:2])[C:1](O)=O', '2': '[N:2][C@@H](CC[N:2])[C:1](O)=O',
                '3': '[N:2][C@@H](CCC[N:2])[C:1](O)=O', '4': '[N:2][C@@H](CCCC[N:2])[C:1](O)=O',
                '5': '[N:2][C@H](C[N:2])[C:1](O)=O', '6': '[N:2][C@H](CC[N:2])[C:1](O)=O',
                '7': '[N:2][C@H](CCC[N:2])[C:1](O)=O', '8': '[N:2][C@H](CCCC[N:2])[C:1](O)=O'}

    AA_SMILES = {'A': '[N:2][C@@H](C)[C:1](O)=O', 'R': '[N:2][C@@H](CCCNC(N)=N)[C:1](O)=O',
                 'N': '[N:2][C@@H](CC(N)=O)[C:1](O)=O', 'D': '[N:2][C@@H](CC(O)=O)[C:1](O)=O',
                 'C': '[N:2][C@@H](CS)[C:1](O)=O', 'Q': '[N:2][C@@H](CCC(N)=O)[C:1](O)=O',
                 'E': '[N:2][C@@H](CCC(O)=O)[C:1](O)=O', 'G': '[N:2]C[C:1](O)=O',
                 'H': '[N:2][C@@H](CC1=CNC=N1)[C:1](O)=O', 'I': '[N:2][C@@H]([C@@H](C)CC)[C:1](O)=O',
                 'K': '[N:2][C@@H](CCCCN)[C:1](O)=O', 'L': '[N:2][C@@H](CC(C)C)[C:1](O)=O',
                 'M': '[N:2][C@@H](CCSC)[C:1](O)=O', 'F': '[N:2][C@@H](CC1=CC=CC=C1)[C:1](O)=O',
                 'P': 'C1CC[N:2][C@@H]1[C:1](O)=O', 'S': '[N:2][C@@H](CO)[C:1](O)=O',
                 'T': '[N:2][C@@H]([C@H](O)C)[C:1](O)=O', 'W': '[N:2][C@@H](CC1=CNC2=CC=CC=C12)[C:1](O)=O',
                 'Y': '[N:2][C@@H](CC1=CC=C(C=C1)O)[C:1](O)=O', 'V': '[N:2][C@@H](C(C)C)[C:1](O)=O'}

    T_SMILES = {'+': '[N:2]'}

    C_SMILES = {'&': 'C[C:1](=O)'}

    # GA class var
    mut_n = 1
    b_insert_rate = 0.1
    selec_strategy = 'Elitist'
    rndm_newgen_fract = 10
    fitness = 'MXFP'

    # initiatization of class variables updated by the GA
    dist_dict_old = {}
    gen_n = 0
    found_identity = 0
    steady_min = 0
    timelimit_seconds = None
    cbd_av = None
    cbd_min = None
    dist_dict = None
    surv_dict = None
    time = 0
    min_dict = {}
    # used internally to recognize a methylated aa:
    metbond = False
    # can be set with exclude or allow methylation,
    # it refers to the possibility of having methylation in the entire GA:
    methyl = True

    # debug
    verbose = False

    def __init__(self):
        pass

    def split_seq_components(self, seq):
        """split seq in generations and branching units

        Arguments:
            seq {string} -- dendrimer sequence

        Returns:
            lists -- generations(gs, from 0 to..), branching units, terminal and capping
        """

        g = []
        gs = []
        bs = []
        t = []
        c = []

        for ix, i in enumerate(seq):
            if i not in ['1', '2', '3', '4', '5', '6', '7', '8']:
                if i in self.CT:
                    t.append(i)
                elif i in self.NT:
                    c.append(i)
                elif i == 'X':
                    continue
                elif i == '-':
                    if seq[ix - 1] in ['1', '2', '3', '4', '5', '6', '7', '8']:
                        bs.append(i)
                    else:
                        g.append(i)
                else:
                    g.append(i)
            else:
                gs.append(g[::-1])
                bs.append(i)
                g = []
        gs.append(g[::-1])
        gs = gs[::-1]
        bs = bs[::-1]

        return gs, bs, t, c

    def connect_mol(self, mol1, mol2):
        """it is connecting all Nterminals of mol1 with the Cterminal 
        of the maximum possible number of mol2s

        Arguments:
            mol1 {rdKit mol object} -- first molecule to be connected
            mol2 {rdKit mol object} -- second molecule to be connected

        Returns:
            rdKit mol object -- mol1 updated (connected with mol2, one or more)
        """
        count = 0

        # detects all the N terminals in mol1
        for atom in mol1.GetAtoms():
            atom.SetProp('Cterm', 'False')
            atom.SetProp('methyl', 'False')
            if atom.GetSmarts() == '[N:2]' or atom.GetSmarts() == '[NH2:2]' or atom.GetSmarts() == '[NH:2]':
                count += 1
                atom.SetProp('Nterm', 'True')
            else:
                atom.SetProp('Nterm', 'False')

        # detects all the C terminals in mol2 (it should be one)
        for atom in mol2.GetAtoms():
            atom.SetProp('Nterm', 'False')
            atom.SetProp('methyl', 'False')
            if atom.GetSmarts() == '[C:1]' or atom.GetSmarts() == '[CH:1]':
                atom.SetProp('Cterm', 'True')
            else:
                atom.SetProp('Cterm', 'False')

        # mol2 is addes to all the N terminal of mol1
        for i in range(count):
            combo = rdmolops.CombineMols(mol1, mol2)
            Nterm = []
            Cterm = []

            # saves in two different lists the index of the atoms which has to be connected
            for atom in combo.GetAtoms():
                if atom.GetProp('Nterm') == 'True':
                    Nterm.append(atom.GetIdx())
                if atom.GetProp('Cterm') == 'True':
                    Cterm.append(atom.GetIdx())

            # creates the amide bond
            edcombo = rdchem.EditableMol(combo)
            edcombo.AddBond(Nterm[0], Cterm[0],
                            order=Chem.rdchem.BondType.SINGLE)
            edcombo.RemoveAtom(Cterm[0] + 1)
            clippedMol = edcombo.GetMol()

            # removes tags and lables form c term atoms which reacted
            clippedMol.GetAtomWithIdx(Cterm[0]).SetProp('Cterm', 'False')
            clippedMol.GetAtomWithIdx(Cterm[0]).SetAtomMapNum(0)

            # methylates amide bond
            if self.metbond == True and self.methyl == True:
                Nterm = []
                Met = []
                methyl = rdmolfiles.MolFromSmiles('[C:4]')
                for atom in methyl.GetAtoms():
                    atom.SetProp('methyl', 'True')
                    atom.SetProp('Nterm', 'False')
                    atom.SetProp('Cterm', 'False')
                metcombo = rdmolops.CombineMols(clippedMol, methyl)
                for atom in metcombo.GetAtoms():
                    if atom.GetProp('Nterm') == 'True':
                        Nterm.append(atom.GetIdx())
                    if atom.GetProp('methyl') == 'True':
                        Met.append(atom.GetIdx())
                metedcombo = rdchem.EditableMol(metcombo)
                metedcombo.AddBond(
                    Nterm[0], Met[0], order=Chem.rdchem.BondType.SINGLE)
                clippedMol = metedcombo.GetMol()
                clippedMol.GetAtomWithIdx(Met[0]).SetProp('methyl', 'False')
                clippedMol.GetAtomWithIdx(Met[0]).SetAtomMapNum(0)

            # removes tags and lables form the atoms which reacted
            clippedMol.GetAtomWithIdx(Nterm[0]).SetProp('Nterm', 'False')
            clippedMol.GetAtomWithIdx(Nterm[0]).SetAtomMapNum(0)

            # uptades the 'core' molecule
            mol1 = clippedMol
        self.metbond = False
        return mol1

    def smiles_from_seq(self, seq):
        """Calculates the smiles of a given peptide dendrimer sequence

        Arguments:
            seq {string} -- peptide dendrimer sequence
        Returns:
            string -- molecule_smile - SMILES of the peptide
        """

        gs, bs, terminal, capping = self.split_seq_components(seq)

        # modifies the Cterminal
        if terminal:
            molecule = rdmolfiles.MolFromSmiles(self.T_SMILES[terminal[0]])
        else:
            molecule = ''

        # creates the dendrimer structure
        for gen in gs:
            for aa in gen:
                if aa == '-':
                    self.metbond = True
                    continue
                if molecule == '':
                    molecule = rdmolfiles.MolFromSmiles(self.AA_SMILES[aa])
                else:
                    molecule = self.connect_mol(
                        molecule, rdmolfiles.MolFromSmiles(self.AA_SMILES[aa]))

            if bs:
                if bs[0] == '-':
                    self.metbond = True
                    bs.pop(0)
                if molecule == '':
                    molecule = rdmolfiles.MolFromSmiles(self.B_SMILES[bs[0]])
                else:
                    molecule = self.connect_mol(
                        molecule, rdmolfiles.MolFromSmiles(self.B_SMILES[bs[0]]))
                bs.pop(0)

        # adds capping to the N-terminal (the called clip function is different, cause the listed smiles
        # for the capping are already without OH, it is not necessary removing any atom after foming the new bond)
        if capping:
            molecule = attach_capping(
                molecule, rdmolfiles.MolFromSmiles(self.C_SMILES[capping[0]]))

        # clean the smile from all the tags
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)

        molecule_smile = rdmolfiles.MolToSmiles(
            molecule, isomericSmiles=True).replace('[N]', 'N').replace('[C]', 'C')
        return molecule_smile
