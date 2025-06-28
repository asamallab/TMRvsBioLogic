import numpy as np
import pickle as pkl
import pandas as pd

def int_to_DNF(k, BF_int, node_list, simplify = False):
    '''
    k: number of inputs to a node
    BF_int: Boolean function represented as a decimal number
    node_list: the list of nodes in Boolean Network
    returns a DNF expression of the BF
    '''
    f = np.binary_repr(BF_int, 2**k)
    indices = [bin(i)[2:].zfill(k) for i, bit in enumerate(f) if bit == '1']
    dnf = ''
    for ele in indices:
        s = ''
        for i in range(k):
            if ele[i] == '1':
                s += node_list[i]+' & '
            else:
                s += '!'+ node_list[i] +' & '
        s = s[:-3]
        dnf += '(' + s + ')' + ' | '
    dnf = dnf[:-3]

    if simplify == False:
        return dnf
    else:
        new_dnf = dnf.replace('!', '~')
        simplify_dnf = str(simplify_logic(parse_expr(new_dnf), 'dnf'))
        return simplify_dnf.replace('~', '!')

def ints_to_BNet(fname, node_num, in_edges, BF_int_list, simplify=False):
    '''
    path: Location to which the bnet file is to be saved
    Converts a list of boolean functions into a bnet file
    '''
    bnet = ['targets' + ',\t' + 'factors' + '\n']
    node_num = {node_num[ele]:ele for ele in node_num}
    for num in node_num:
        if len(in_edges[num]) == 0:
            bnet += [node_num[num] + ',\t' + f'({node_num[num]})' + '\n']
        else:
            inp_genes = [node_num[index] for index in in_edges[num]]
            bnet += [node_num[num] + ',\t' + int_to_DNF(len(in_edges[num]), BF_int_list[num], inp_genes, simplify) + '\n']
    with open(fname, 'w') as file:
        file.writelines(bnet)
    return f'BoolNet file saved at {fname}...'