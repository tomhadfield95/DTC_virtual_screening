from collections import OrderedDict
from operator import itemgetter
import os
import pickle
import random
import tempfile

#import deepchem as dc
import numpy as np
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
import torch
from torch.autograd import Variable
import torch.nn as nn
from copy import deepcopy

ATOM_SYMBOLS = ['C', 'N', 'O', 'S', 'F', 'P', 'Cl', 'Br', 'D']

def create_var(tensor, requires_grad=False): 
    """\
    create_var(...) -> torch.autograd.Variable

    Wrap a torch.Tensor object by torch.autograd.Variable.
    """
    return Variable(tensor, requires_grad=requires_grad)

def one_of_k_encoding_unk(x, allowable_set):
    #"""Maps inputs not in the allowable set to the last element."""
    """\
    one_of_k_encoding_unk(...) -> list[int]

    One-hot encode `x` based on `allowable_set`.
    Return None if `x not in allowable_set`.
    """
    if x not in allowable_set:
        #x = allowable_set[-2]
        return None
    return list(map(lambda s: int(x == s), allowable_set))

def atom_features(atom, include_extra = False):
    """\
    atom_features(...) -> list[int]

    One-hot encode `atom` w/ or w/o extra concatenation.
    """
    retval  = one_of_k_encoding_unk(atom.GetSymbol(), ATOM_SYMBOLS)
    if include_extra:
        retval += [#atom.GetDegree(),
                   atom.GetFormalCharge()]
    return retval

def bond_features(bond, include_extra = False):
    """\
    bond_features(...) -> list[int]

    One-hot encode `bond` w/ or w/o extra concatenation.

    Parameters
    ----------
    bond: rdkit.Chem.Bond
    include_extra: bool
    """
    bt = bond.GetBondType()  # rdkit.Chem.BondType
    retval = [
      bt == Chem.rdchem.BondType.SINGLE, bt == Chem.rdchem.BondType.DOUBLE,
      bt == Chem.rdchem.BondType.TRIPLE, bt == Chem.rdchem.BondType.AROMATIC,
      0 # no bond
      #bond.GetIsConjugated(),
      #bond.IsInRing()
      ]
    if include_extra:
        bs = bond.GetStereo()
        retval += [bs == Chem.rdchem.BondStereo.STEREONONE,
                   bs == Chem.rdchem.BondStereo.STEREOANY,
                   bs == Chem.rdchem.BondStereo.STEREOZ,
                   bs == Chem.rdchem.BondStereo.STEREOE,
                   bs == Chem.rdchem.BondStereo.STEREOCIS,
                   bs == Chem.rdchem.BondStereo.STEREOTRANS
                  ]
    return np.array(retval)


def make_graph(smiles, extra_atom_feature = False, extra_bond_feature = False):
    """\
    make_graph(...) -> edge_dict, node_dict

    Construct a graph representation of a given SMILES.

    The nodes and edges are represented as one-hot vectors,
    which are respectively gathered into a node dict and an edge dict
    having atom indices as their keys.
    See below for the exact return structures.
    This graph representation is constantly used in the methods of `ggm.ggm`.

    Returns
    -------
    g: OrderedDict[int, list[tuple[torch.autogra.Variable, int]]] | None
        Edge (bond) dictionary, which looks like:
            { atom_idx: [ (one_hot_vector, partner_atom_idx), ... ], ... }
    h: OrderedDict[int, torch.autograd.Variable] | None
        Node (atom) dictionary, which looks like:
            { atom_idx: one_hot_vector, ... }
    If untreatable atoms are present, return (None, None).
    """
    g = OrderedDict({})
    h = OrderedDict({})
    if type(smiles) is str or type(smiles) is np.str_:
        molecule = Chem.MolFromSmiles(smiles)
    else:
        molecule = smiles
    Chem.Kekulize(molecule)

    chiral_list = Chem.FindMolChiralCenters(molecule)
    chiral_index = [c[0] for c in chiral_list]
    chiral_tag = [c[1] for c in chiral_list]
    #Chem.Kekulize(molecule)
    #Chem.Kekulize(molecule, clearAromaticFlags=False)
    for i in range(0, molecule.GetNumAtoms()):
        atom = molecule.GetAtomWithIdx(i)  # rdkit.Chem.Atom
        if atom.GetSymbol() not in ATOM_SYMBOLS:
            return None, None
        atom_i = atom_features(atom, extra_atom_feature)  # One-hot vector
        if extra_atom_feature:
            if i in chiral_index:
                if chiral_tag[chiral_index.index(i)]=='R':
                    atom_i += [0, 1, 0]
                else:
                    atom_i += [0, 0, 1]
            else:
                atom_i += [1, 0, 0]
            atom_i.append(atom.GetIsAromatic())

        h[i] = create_var(torch.FloatTensor(atom_i), False).view(1, -1)
        for j in range(0, molecule.GetNumAtoms()):
            e_ij = molecule.GetBondBetweenAtoms(i, j)  # rdkit.Chem.Bond
            if e_ij != None:
                e_ij = list(map(lambda x: 1 if x == True else 0, bond_features(e_ij, extra_bond_feature))) # ADDED edge feat; one-hot vector
                e_ij = create_var(torch.FloatTensor(e_ij).view(1, -1), False)
                atom_j = molecule.GetAtomWithIdx(j)
                if i not in g:
                    g[i] = []
                g[i].append( (e_ij, j) )
    return g, h

def sum_node_state(h):   
    """Return the element-wise sum of the node vectors."""
    retval = create_var(torch.zeros(h[0].size()))
    for i in list(h.keys()):
        retval+=h[i]
    return retval

def average_node_state(h):   
    """Return the element-wise mean of the node vectors."""
    #retval = create_var(torch.zeros(1,9))
    retval = create_var(torch.zeros(h[0].size()))
    for i in list(h.keys()):
        retval+=h[i]
    if len(h)>0:
        retval=retval/len(h)
    return retval

def collect_node_state(h, except_last=False):   
    """Return a matrix made by concatenating the node vectors.
    Return shape -> (len(h), node_vector_length)    # if not except_last
                    (len(h)-1, node_vector_length)  # if except_last
    """
    retval = []
    for i in list(h.keys())[:-1]:
        retval.append(h[i])
    if except_last==False:
        retval.append(h[list(h.keys())[-1]])
    return torch.cat(retval, 0)

def cal_formal_charge(atomic_symbol, bonds) -> int:
    """\
    Compute the formal charge of `atomic_symbol`
    based on its partner atoms and bond orders.

    Parameters
    ----------
    atomic_symbol: str
    bonds: list[tuple[str, int]]
        [ (atom_symbol, bond_order), ... ]
    """
    if atomic_symbol=='N':
        #if sorted(bonds, key=lambda x: (x[0], x[1])) == [('C', 1), ('O', 1), ('O', 2)]:
        #    return 1
        if sum(j for i, j in bonds)==4:
            return 1
    #if atomic_symbol=='O':
    #    if sorted(bonds, key=lambda x: (x[0], x[1])) == [('N', 1)]:
    #        return -1

    return 0

def graph_to_smiles(g, h) -> str:
    """Prepare atom symbols, bond orders and formal charges
    and call `self.BO_to_smiles` to return a SMILES str."""
    # Determine atom symbols by argmax of each node vector.
    atomic_symbols = OrderedDict({})
    for i in h.keys():
        atomic_symbols[i] = ATOM_SYMBOLS[np.argmax(h[i].data.cpu().numpy())]

    idx_to_atom = {i:k for i,k in enumerate(h.keys())}
    atom_to_idx = {k:i for i,k in enumerate(h.keys())}

    # Determine bond orders by argmax of each edge vector.
    BO = np.zeros((len(atomic_symbols), len(atomic_symbols)))
    for i in h.keys():
        for j in range(len(g[i])):
            BO[atom_to_idx[i],atom_to_idx[g[i][j][1]]] = \
                np.argmax(g[i][j][0].data.cpu().numpy())+1
    if not np.allclose(BO, BO.T, atol=1e-8):
        print ('BO is not symmetry')
        exit(-1)
    fc_list = []
    for i in range(len(atomic_symbols)):
        for j in g[idx_to_atom[i]]:
            bond = [(atomic_symbols[j[1]], BO[i][atom_to_idx[j[1]]]) \
                for j in g[idx_to_atom[i]]]
        fc = cal_formal_charge(atomic_symbols[idx_to_atom[i]], bond)
        if fc!=0:
            fc_list.append([i+1, fc])
    #print (fc_list)
    if len(fc_list)==0:
        fc_list = None
    smiles = BO_to_smiles(atomic_symbols, BO, fc_list)
    return smiles

def BO_to_smiles(atomic_symbols, BO, fc_list=None) -> str:
    """\
    Obtain a SMILES str corresponding to the given atom and bond types.

    During the routine, a temporary SDF file is written
    and read by RDKit to get a SMILES.

    Parameters
    ----------
    atomic_symbols: list[str]
    BO: 2D-ndarray of float or int
    fc_list: None | list[int]

    Returns
    -------
    smiles: str | None
        If a valid SMILES cannot be made, None is returned.
    """
    natoms = len(atomic_symbols)
    nbonds = int(np.count_nonzero(BO)/2)
    # Temporary file descriptor and path to write SDF
    sdf_fd, sdf_path = tempfile.mkstemp(prefix='GGM_tmp', dir=os.getcwd(), text=True)
    with open(sdf_fd, 'w') as w:
        w.write('\n')
        w.write('     GGM\n')
        w.write('\n')
        w.write(str(natoms).rjust(3,' ')+str(nbonds).rjust(3, ' ') + '  0  0  0  0  0  0  0  0999 V2000\n')
        for s in atomic_symbols.values():
            w.write('    0.0000    0.0000    0.0000 '+s+'   0  0  0  0  0  0  0  0  0  0  0  0\n')
        for i in range (int(natoms)):
            for j in range(0,i):
                bond_order = BO[i,j]
                if bond_order!=0:
                    #if BO[i,j]==4:
                    #    BO[i,j]=2
                    w.write(str(i+1).rjust(3, ' ') + str(j+1).rjust(3, ' ') + str(int(bond_order)).rjust(3, ' ') + '0'.rjust(3, ' ') + '\n')
        if fc_list is not None:
            w.write('M  CHG  '+str(len(fc_list)))
            for fc in fc_list:
                w.write(str(fc[0]).rjust(4, ' ')+str(fc[1]).rjust(4, ' '))
            w.write('\n')
        w.write('M  END\n')
        w.write('$$$$')
    # Rewrite the SDF using `babel`.
    #os.system(f'babel -isdf {sdf_path} -osdf {sdf_path} 2> {os.devnull}')
    # Get a SMILES if valid.
    try:
        m = Chem.SDMolSupplier(sdf_path)[0]
        s = Chem.MolToSmiles(m)
        # Final validation.
        s = Chem.MolToSmiles(Chem.MolFromSmiles(s))
    finally:
        #pass
        os.unlink(sdf_path)
    return s

def one_hot(tensor, depth):
    """\
    Return an one-hot vector given an index and length.

    Parameters
    ----------
    tensor: torch.FloatTensor of shape (1,)
        A 0-D tensor containing only an index.
    depth: int
        The length of the resulting one-hot.

    Returns
    -------
    torch.FloatTensor of shape (1, depth)
    """
    ones = torch.sparse.torch.eye(depth)
    return ones.index_select(0,tensor.long())


def is_equal_node_type(h1, h2):
    if len(h1)!=len(h2):
        return False
    for i in h1:
        if not np.array_equal(h1[i].data.cpu().numpy(), h2[i].data.cpu().numpy()):
            return False
    return True

def is_equal_edge_type(g1, g2):
    if len(g1)!=len(g2):
        return False
    for i in g1:
        if len(g1[i])!=len(g2[i]):
            return False
        sorted1 = sorted(g1[i], key=itemgetter(1))
        sorted2 = sorted(g2[i], key=itemgetter(1))
        for j in range(len(sorted1)):
            if not np.array_equal(sorted1[j][0].data.cpu().numpy(), sorted2[j][0].data.cpu().numpy()):
                print (i, j, sorted1[j][0].data.cpu().numpy())
                print (i, j, sorted2[j][0].data.cpu().numpy())
                return False
            if sorted1[j][1]!=sorted2[j][1]:
                print (i, j, sorted1[j][1])
                print (i, j, sorted2[j][1])
                return False
    return True         


def ensure_shared_grads(model, shared_model, gpu=False):
    """Make a shared model have the gradients of a worker model.
    Used in `vaetrain.train` for multiprocessed optimization."""
    for param, shared_param in zip(model.parameters(), shared_model.parameters()):
        if shared_param.grad is not None and not gpu:
            return
        elif not gpu:
            shared_param._grad = param.grad
        else:
            if param.grad is None:
                continue
            shared_param._grad = param.grad.cpu()

def probability_to_one_hot(tensor, stochastic = False):
    """\
    Convert a vector to one-hot of the same size.

    If not stochastic, the one-hot index is argmax(tensor).
    If stochastic, an index is randomly chosen
    according to a probability proportional to each element.

    Parameters
    ----------
    tensor: torch.autograd.Variable
    stochastic: bool
    """
    if stochastic:
        # Index-selection probability proportional to each element
        prob = tensor.data.cpu().numpy().ravel().astype(np.float64)
        prob = prob/np.sum(prob)
        norm = np.sum(prob)
        prob = [prob[i]/norm for i in range(len(prob))]
        idx = int(np.random.choice(len(prob), 1, p = prob))
    else:
        idx = int(np.argmax(tensor.data.cpu().numpy()))
    return create_var(one_hot(torch.FloatTensor([idx]), list(tensor.size())[-1] ))

def make_graphs(s1, s2, extra_atom_feature = False, extra_bond_feature = False):
    """\
    make_graphs(...) -> edge_dict_1, node_dict_1, edge_dict_2, node_dict_2

    Make graphs of `s1` and `s2` using `make_graph`.

    `s2` must be a substructure of `s1`; otherwise None's are returned.
    A graph of `s1` is first made, from which that of `s2` is extracted
    so that the indices of the shared atoms are the same.
    """
    molecule1 = Chem.MolFromSmiles(s1)
    molecule2 = Chem.MolFromSmiles(s2)
    #Chem.Kekulize(molecule1, clearAromaticFlags=False)
    #Chem.Kekulize(molecule2, clearAromaticFlags=False)
    g1, h1 = make_graph(Chem.Mol(molecule1), extra_atom_feature, extra_bond_feature)
    if g1 is None or h1 is None:
        return None, None, None, None

    scaffold_index = molecule1.GetSubstructMatches(molecule2)
    if len(scaffold_index)==0 :
        return None, None, None, None
    scaffold_index = list(scaffold_index[0])

    g2 = OrderedDict({})
    h2 = OrderedDict({})

    for i in h1:
        if i in scaffold_index:
            h2[i] = deepcopy(h1[i])

    for i in g1:
        if i not in scaffold_index :
            continue
        else:
            g2[i] = []
        for edge in g1[i]:
            if edge[1] in scaffold_index:
                g2[i].append(deepcopy(edge))

    return g1, h1, g2, h2

def index_rearrange(molecule1, molecule2, g, h):
    """\
    index_rearrange(...) -> edge_dict, node_dict

    By comparing `molecule1` and a substructure `molecule2`,
    make the overlapping atoms of `molecule1` have the same indices as in `molecule2`:

        molecule1   molecule2      molecule1
        H - N - C     N - C    ->  H - N - C
        0   1   2     1   0        2   1   0

    and apply the rearrangement to `g` and `h`.
    Note that the order of `g.values()` and `h.values()` is preserved.

    Parameters
    ----------
    molecule1: rdkit.Chem.Mol
    molecule2: rdkit.Chem.Mol
    g: edge dict of molecule1
    h: node dict of molecule1
    """
    #Chem.Kekulize(molecule1)
    #Chem.Kekulize(molecule2)
    # The indices of `molecule1` atoms that overlap `molecule2` atoms.
    # The returned ordering corresponds to the atom ordering of `molecule2`.
    scaffold_index = list(molecule1.GetSubstructMatches(molecule2)[0])
    new_index = OrderedDict({})  # Does `new_index` have to be an "ordered" dict?
    for idx,i in enumerate(scaffold_index):
        new_index[i]=idx  # new_index[index_in_molecule1] -> index_in_molecule2
    # Shift to the end the indices of the left atoms in `molecule1`.
    idx = len(scaffold_index)
    for i in range(len(h)):
        if i not in scaffold_index:
            new_index[i] = idx
            idx+=1
    g, h = index_change(g, h, new_index)
    return g, h
        
def index_change(g, h, new_index):
    """\
    index_change(...) -> edge_dict, node_dict

    Rearrange the node numbering of `g` and `h` according to the mapping by `new_index`.

    Parameters
    ----------
    new_index: dict[int, int]
    """
    new_h = OrderedDict({})
    new_g = OrderedDict({})
    for i in h.keys():
        new_h[new_index[i]]=h[i]
    for i in g.keys():
        new_g[new_index[i]]=[]
        for j in g[i]:
            new_g[new_index[i]].append((j[0], new_index[j[1]]))
    return new_g, new_h            

def enumerate_molecule(s: str):
    """Return a list of all the isomer SMILESs of a given SMILES `s`."""
    m = Chem.MolFromSmiles(s) 
    opts = StereoEnumerationOptions(unique=True, onlyUnassigned=False)
    #opts = StereoEnumerationOptions(tryEmbedding=True, unique=True, onlyUnassigned=False)
    isomers = tuple(EnumerateStereoisomers(m, options=opts)) 
    retval = []
    for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):  
        retval.append(smi)
    return retval

def initialize_model(model, load_save_file=False):
    """\
    Parameters
    ----------
    model: ggm.ggm
    load_save_file: str
        File path of the save model.
    """
    if load_save_file:
        model.load_state_dict(torch.load(load_save_file))
    else:
        for param in model.parameters():
            if param.dim() == 1:
                continue
                nn.init.constant(param, 0)
            else:
                #nn.init.normal(param, 0.0, 0.15)
                nn.init.xavier_normal(param)
    return model    

def dict_from_txt(path, dtype=None):
    """\
    Generate a dict from a text file.

    The structure of `path` should be

        key1  value1_1  value1_2  ...
        key2  value2_1  value2_2  ...
        ...

    and then the returned dict will be like

        { key1:[value1_1, value1_2, ...], ... }

    NOTE that the dict value will always be a list,
    even if the number of elements is less than 2.

    Parameters
    ----------
    path: str
        A data text path.
    dtype: type | None
        The type of values (None means str).
        Keys will always be str type.

    Returns
    -------
    out_dict: dict[str, list[dtype]]
    """
    out_dict = {}
    # np.genfromtxt or np.loadtxt are slower than pure Python!
    with open(path) as f:
        for line in f:
            row = line.split()
            if dtype is None:
                out_dict[row[0]] = row[1:]
            else:
                out_dict[row[0]] = [dtype(value) for value in row[1:]]
    return out_dict

def load_data(smiles_path, *data_paths):
    """\
    Read data text files.

    The structure of `smiles_path` should be 

        #    whole      scaffold
        ID1  SMILES1_1  SMILES1_2
        ID2  SMILES2_1  SMILES2_2
        ...

    And the structure of EACH path in `data_paths` should be

        #    whole               scaffold
        ID1  value_of_SMILES1_1  value_of_SMILES1_2
        ID2  value_of_SMILES2_1  value_of_SMILES2_2
        ...

    Parameters
    ----------
    smiles_path: str
        A data text path of IDs and SMILESs.
    data_paths: iterable of str
        Each is a data text path of IDs and property values.

    Returns
    -------
    id_to_smiles: dict[str, list[str]]
    id_to_whole_conditions: dict[str, list[float]]
        { ID1: [whole_value_of_property1, whole_value_of_property2, ...], ... }
    id_to_scaffold_conditions: dict[str, list[float]]
        { ID1: [scaffold_value_of_property1, scaffold_value_of_property2, ...], ... }
    """
    id_to_smiles = dict_from_txt(smiles_path)  # {id:[whole_smiles, scaffold_smiles], ...}
    data_dicts = [dict_from_txt(path, float) for path in data_paths]

    # Only collect molecule IDs having every label.
    common_keys = set(id_to_smiles.keys()).intersection(
        *(data_dict.keys() for data_dict in data_dicts)
    )

    # Collect condition values of multiple properties.
    id_to_whole_conditions = {}
    id_to_scaffold_conditions = {}
    for key in common_keys:
        id_to_whole_conditions[key] = [data_dict[key][0] for data_dict in data_dicts]
        id_to_scaffold_conditions[key] = [data_dict[key][1] for data_dict in data_dicts]
    return id_to_smiles, id_to_whole_conditions, id_to_scaffold_conditions

def divide_data(id_to_conditions, boundaries=.5):
    """\
    Divide data by boundary values.

    Used inside `sample_data`.
    If `id_to_conditions.values()` are all empty, ([], []) is returned.

    Parameters
    ----------
    id_to_conditions: dict[str, list[float]]
        { ID1: [value_of_property1, value_of_property2, ...], ... }
    boundaries: float | list[float]
        (A list of) the boundary value of EACH property.
        The value needs not be in [0, 1] depending on its range.

    Returns
    -------
    high_keys: list[ list[str] ]
        [ [keys_of_high_property1_values...], [keys_of_high_property2_values...] ]
    low_keys: list[ list[str] ]
        [ [keys_of_low_property1_values...], [keys_of_low_property2_values...] ]
    """
    # Preprocess
    num_properties = len(next(iter(id_to_conditions.values())))  # Get one value.
    if isinstance(boundaries, float):
        boundaries = [boundaries for _ in range(num_properties)]

    high_keys = [[]] *num_properties
    low_keys = [[]] *num_properties
    for key, values in id_to_conditions.items():
        for i in range(num_properties):
            if values[i] > boundaries[i]:
                high_keys[i].append(key)
            else:
                low_keys[i].append(key)
    return high_keys, low_keys

def sample_data(id_to_conditions, size, ratios=.5, boundaries=.5):
    """\
    Sample high-valued and low-valued keys of multiple properties by some ratios.

    If `id_to_conditions.values()` are all empty, random sampling is done.

    Parameters
    ----------
    id_to_conditions: dict[str, list[float]]
        { ID1: [value_of_property1, value_of_property2, ...], ... }
    size: int
        The number of samples.
    ratio: flaot | list[float]
        (A list of) the ratio of high-valued samples of EACH property.
    boundaries: float | list[float]
        (A list of) the boundary value of EACH property.
        The value needs not be in [0, 1] depending on its range.

    Returns
    -------
    keys: list[str]
        A list of sampled keys.
    """
    # Preprocess
    random.seed()
    high_keys, low_keys = divide_data(id_to_conditions, boundaries)
    num_properties = len(next(iter(id_to_conditions.values())))  # Get one value.
    # Return random samples if no property is present.
    if not num_properties:
        return random.sample(list(id_to_conditions.keys()), k=size)
    num_each = size // num_properties
    if isinstance(ratios, float):
        ratios = [ratios for _ in range(num_properties)]

    # Sample.
    key_pool = np.array(high_keys[0] + low_keys[0])
    keys = []
    for i in range(num_properties):
        sample_idxs = random.choices(
            np.nonzero(np.in1d(key_pool, high_keys[i]))[0],
            k = int(num_each * ratios[i])
        )
        sample_idxs += random.choices(
            np.nonzero(np.in1d(key_pool, low_keys[i]))[0],
            k = num_each - len(sample_idxs)
        )
        keys.extend(key_pool[sample_idxs])
        key_pool = np.delete(key_pool, sample_idxs)

    # Additional sampling
    keys.extend(random.choices(key_pool, k=size-len(keys)))
    random.shuffle(keys)
    return keys
