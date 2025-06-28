import numpy as np
import networkx as nx
from sympy import simplify_logic, parse_expr

class BN:
    '''
    This code should give all the basic information about the Boolean network:
    -names of the nodes,
    -edges,
    -attractors names and states,
    -basin sizes,
    -inedgree per node
    -out degree per each node
    -Convert our file to BoolNet
    '''
    def __init__(self, inedges, BF_list, bio_fp):
        self.inedges = inedges
        self.BF_list = BF_list
        self.bio_fp = bio_fp
        self.n = len(self.inedges)
    
    def IntsToBitsFast (ints, size):
        '''
        return an array of bits of length 'size' for the set of integers ints such that the MSB is at column 0
        source: https://stackoverflow.com/questions/
        ints: numpy array of integers
        size: required size of the array of bits
        '''
        return (((ints[:,None] & (1 << np.arange(size, dtype = 'int64'))[::-1])) > 0).astype(int)
        

    def IntsToBitsFast_large (int_num, size):
        '''
        returns an array of bits of length 'size' such that the MSB is at column 0
        source: https://stackoverflow.com/questions/
        int_num : integer
        size: required size of the array of bits
        '''
        b = [int_num >> i & 1 for i in range(int_num.bit_length() - 1,-1,-1)]
        if len(b) != size:
            b = (size-len(b))*[0] + b
        return np.array([b])

    def BitsToIntsFast(bits):
        '''
        returns an array of integers
        bits: a matrix of bits whose rows represent the integer of interest
        from: https://stackoverflow.com/questions/41069825/convert-binary-01-numpy-to-integer-or-binary-string
        '''
        m,n = bits.shape # number of columns is needed, not bits.size
        a = 2**np.arange(n)[::-1]  # -1 reverses array of powers of 2 of same length as bits
        return bits @ a  # this matmult is the key line of code

    def STV (self):
        '''
        returns the state transition vector
        PROVIDE example
        '''
        x = BN.IntsToBitsFast (np.arange(2**self.n), self.n)
        new_mat = np.zeros([2**self.n, self.n], dtype = int)
        for node in self.inedges:
            q = self.inedges[node]
            ind = BN.BitsToIntsFast(x[:,q])
            f = BN.IntsToBitsFast_large(self.BF_list[node], 2**len(q))
            new_mat[:,node] = f[:,ind]
        stv = BN.BitsToIntsFast(new_mat)
        return stv

    def fxd_pts (self, stv):
        '''

        '''
        return np.where(np.arange(len(stv)) == stv)[0]

    def basin_fp (self, fp, stv):
        '''

        '''
        stg = nx.DiGraph()
        stg.add_edges_from (list(zip(np.arange(len(stv)), stv)))
        return list(nx.ancestors(stg, fp)) + [fp]

    def basin_sizes_fps (self, stv):
        '''

        '''
        stg = nx.DiGraph()
        stg.add_edges_from (list(zip(np.arange(len(stv)), stv)))
        boa = {}
        for fp in self.bio_fp:
            boa[fp] = len(nx.ancestors(stg, fp)) + 1
        D = {k: v for k, v in sorted(boa.items(), key=lambda item:item[1])}
        return D

    def func_to_dnf (k, func, genes, simp = False):
        '''
        k: number of inputs to a node
        function: Boolean function represented as a decimal number
        genes: the list of genes in the network
        returns a DNF expression of the BF
        '''
        f = np.binary_repr(func, 2**k)
        indices = [bin(i)[2:].zfill(k) for i, bit in enumerate(f) if bit == '1']
        dnf = ''
        for ele in indices:
            s = ''
            for i in range(k):
                if ele[i] == '1':
                    s += genes[i]+' & '
                else:
                    s += '!'+ genes[i] +' & '
            s = s[:-3]
            dnf += '(' + s + ')' + ' | '
        dnf = dnf[:-3]

        if simp == False:
            return dnf
        else:
            new_dnf = dnf.replace('!', '~')
            simp_dnf = str(simplify_logic(parse_expr(new_dnf), 'dnf'))
            return simp_dnf.replace('~', '!')

    def to_BoolNet (path, node_num, inedges, BF_list, mod_num, simp = False):
        '''
        path: Location to which the bnet file is to be saved
        Converts a list of boolean functions into a BoolNet file
        '''
        bnet = ['targets' + ',\t' + 'factors' + '\n']
        node_num = {node_num[ele]:ele for ele in node_num}
        for num in node_num:
            inp_genes = [node_num[index] for index in inedges[num]]
            bnet += [node_num[num] + ',\t' + BN.func_to_dnf(len(inedges[num]), BF_list[num], inp_genes, simp) + '\n']
        with open(path + 'model' + str(mod_num) +'.bnet', 'w') as file:
            file.writelines(bnet)
        return 'BoolNet file saved at', path +'model'+ str(mod_num)+ '.bnet'

class RS(BN):
    '''
    '''
    def __init__(self, inedges, BF_list, bio_fp):
        super().__init__(inedges, BF_list, bio_fp)