import numpy as np
import pandas as pd
from mod_sel_BN_for_Derrida import *
from BF_properties import *
import pickle as pkl
import ast
import multiprocessing
import time

def IntsToBitsFast_lr (ints, size):
    '''
    return an array of bits of length 'size' for the set of integers ints such that the MSB is at column 0
    source: https://stackoverflow.com/questions/
    ints: numpy array of integers
    size: required size of the array of bits
    '''
    return (((ints[:,None] & (1 << np.arange(size, dtype = 'int64'))[::-1])) > 0).astype(int)

def model_data_random(edges, all_nodes, node_num):
    '''
    edges: It is a dataframe with three columns: 'network_number', 'fro' and 'to' where 'fro' and 'to' corresponds 
    to vector of nodes from which an edge emanates and terminates respectively.
    all_nodes: This is a list of nodes ([gene0, gene1, gene2,...])
    node_num: This is a dictionary that links genes with an index ({'gene0': 0, 'gene1': 1, 'gene2': 2, 'gene3': 3,
    'gene4': 4, 'gene5': 5, ...})
    returns information about the interaction in the network
        -in_edges:  a dictionary where each gene is associated with its
        regulatory inputs ({0: [11, 1], 1: [8, 6], 2: [10, 11], 3: [0, 11], 4: [5, 11], ....})
    '''
    in_edges = dict()
    for node in all_nodes:
        inter = edges.loc[edges["to"] == node, ["fro"]]
        in_edges[node_num[node]] = [node_num[ele] for ele in inter['fro']]
        if in_edges[node_num[node]] == []:
            in_edges[node_num[node]] = [node_num[node]]
            print (f'Node {node} does not have any inputs. Self loop assigned!')
    return in_edges



def filter_by_zero_bits(pairs, msb_indices, width):
   
    mask = 0
    for msb_i in msb_indices:
        lsb_pos = width - 1 - msb_i
        mask |= (1 << lsb_pos)
    return [(x, y) for x, y in pairs if (x & mask) == 0 and (y & mask) == 0]


# Network sensitivity or Derrida coefficient for a single model
def network_sensitivity(inedges, BF_list):
    N = len(inedges)
    tot_net_sen = 0
    for i in range(N):
        inp_ith_node = len(inedges[i])
        bf_str = bin(BF_list[i])[2:].zfill(2**inp_ith_node)
        avg_sen = bf(inp_ith_node, bf_str).avg_sensitivity()
        tot_net_sen += avg_sen
    ave_net_sen = tot_net_sen/N
    return ave_net_sen



def final_Hamming_distance(inedges,BF_list,t_0, t_inf, one_neighbors_list):
    bio_fp=[]
    N = len(inedges)
    net_rs = RS(inedges, BF_list, bio_fp)
    stv = net_rs.STV()
    list_after_evol = one_neighbors_list
    sum_ham_array = np.zeros(len(one_neighbors_list))
    for step in range(t_inf):
        list_after_evol = stv[list_after_evol]
        if step >= t_0:
            bitwise_xor_result = np.bitwise_xor.reduce(list_after_evol, axis = 1)
            ham_dist_list = np.vectorize(lambda x: bin(x).count('1'))(bitwise_xor_result)
            sum_ham_array += ham_dist_list
    ave_ham_array = sum_ham_array/(t_inf-t_0)
    final_Hamming_dist = np.mean(ave_ham_array)
    return final_Hamming_dist


def get_derrida_coefficient(inedges,BF_list, one_neighbors_list):
    N = len(inedges)
    bio_fp=[]
    net_rs = RS(inedges, BF_list, bio_fp)
    stv = net_rs.STV()
    list_after_evol = stv[one_neighbors_list]
    bitwise_xor_result = np.bitwise_xor.reduce(list_after_evol, axis = 1)
    ham_dist_list = np.vectorize(lambda x: bin(x).count('1'))(bitwise_xor_result)
    derrida_coeff = np.mean(ham_dist_list)
    return derrida_coeff




def Fragility(inedges,BF_list,t_0, t_inf, one_neighbors_list):
    bio_fp=[]
    N = len(inedges)
    net_rs = RS(inedges, BF_list, bio_fp)
    stv = net_rs.STV()
    x_t_vect = IntsToBitsFast_lr (np.zeros(len(one_neighbors_list),dtype=int), N)
    x_t_perturb_vect = IntsToBitsFast_lr (np.zeros(len(one_neighbors_list),dtype=int), N)
    list_after_evol = list(map(list, zip(*one_neighbors_list)))
    for step in range(t_inf):
        list_after_evol = stv[list_after_evol]
        if step >= t_0:
            x_t_vect += IntsToBitsFast_lr (list_after_evol[0], N)
            x_t_perturb_vect += IntsToBitsFast_lr (list_after_evol[1], N)
    ave_x_t_vect = x_t_vect/(t_inf-t_0)
    ave_x_t_perturb_vect = x_t_perturb_vect/(t_inf-t_0)
    fragility = np.sum(np.abs(ave_x_t_vect - ave_x_t_perturb_vect), axis=1).mean()
    return fragility



def append_rows_to_tsv(networks,models,func_types,derrida_coeffs, net_sens, final_hammings, fragilities, file_path, columns):
    rows = list(zip(networks,models,func_types,derrida_coeffs, net_sens, final_hammings, fragilities))
    df = pd.DataFrame(rows, columns=columns)
    df.to_csv(file_path, sep='\t', index=False, mode='a', header=False)



def Dataframe_of_three_measures(n,k,func_type,t_0,t_inf,net_init,net_final, model_per_net):
    columns = ['Network_no', 'Model_no', 'func_type', 'derrida_coeff', 'net_sen', 'final_hamming', 'fragility']
    df = pd.DataFrame(columns=columns)
    file_path = f'../output/{n}_{k}_{func_type}_{net_init}_{net_final}_derrida_PP.tsv'
    df.to_csv(file_path, sep='\t', index=False)
    edges_list = pd.read_csv(f'../../RBNs_and_the_models/networks//RN_{n}_{k}_10000_PP.tsv', sep = '\t') # Contains the interactions for all the networks generated in R (BoolNet)
    model_df = pd.read_csv(f'../../RBNs_and_the_models/models_scNCF_IMR_and_BMR/{n}_{k}_10000_models_NCF_IMR_BMR_PP.tsv', sep = '\t')
    with open(f'../input/one_ham_neib_list_for_{n}_nodes_no_repeat.pkl', 'rb') as file:
        one_neighbors_list_original = pkl.load(file)

    edges_list['fro'] = edges_list['fro'].apply(lambda x: f'gene{x-1}') # We replace the numbers by name: x -> gene{x-1}
    edges_list['to'] = edges_list['to'].apply(lambda x: f'gene{x-1}') # for both the 'fro' and 'to' column
    edges_list['network_number'] = edges_list['network_number'] - 1
    node_num = {'gene'+str(i): i for i in range(n)}
    all_nodes = list(node_num.keys())
    model_df_relev = model_df[model_df['func_type']==func_type]
    net_indices = range(net_init, net_final)
    for network_number in net_indices:
            networks = [network_number] * model_per_net
            models = list(range(1))
            func_types = [func_type] * model_per_net
            derrida_coeffs, net_sens, final_hammings, fragilities = [], [], [], []
            edges = edges_list[edges_list['network_number'] == network_number] # Edges gives the interactions df for a single network
            in_edges = model_data_random(edges, all_nodes, node_num)
            
            models_for_network_df = model_df_relev[model_df_relev['Network_no'] == network_number]
            models_for_network_df.loc[:, 'model'] = models_for_network_df['model'].apply(ast.literal_eval)
            
            if func_type == 'scNCF':
                for j in range(model_per_net):
                    BF_list = models_for_network_df['model'].iloc[j]
                    derrida_coeffs.append(get_derrida_coefficient(in_edges, BF_list,one_neighbors_list_original))
                    net_sens.append(network_sensitivity(in_edges, BF_list))
                    final_hammings.append(final_Hamming_distance(in_edges,BF_list,t_0, t_inf, one_neighbors_list_original))
                    fragilities.append(Fragility(in_edges,BF_list,t_0, t_inf, one_neighbors_list_original))
                    
            
            
            if func_type == 'IMR':
                for k1 in in_edges:
                    if len(in_edges[k1]) % 2 == 0 and k1 not in in_edges[k1]:
                        in_edges[k1].append(k1)
                        
                for j in range(model_per_net):
                    BF_list = models_for_network_df['model'].iloc[j]
                    derrida_coeffs.append(get_derrida_coefficient(in_edges, BF_list, one_neighbors_list_original))
                    net_sens.append(network_sensitivity(in_edges, BF_list))
                    final_hammings.append(final_Hamming_distance(in_edges,BF_list,t_0, t_inf, one_neighbors_list_original))
                    fragilities.append(Fragility(in_edges,BF_list,t_0, t_inf, one_neighbors_list_original))
                
                
                
            if func_type == 'BMR':
                for k1 in in_edges:
                    if k1 not in in_edges[k1]:
                        in_edges[k1].append(k1)
                        
                for j in range(model_per_net):
                    BF_list = models_for_network_df['model'].iloc[j]
                    msb_indices = [i for i, val in enumerate(BF_list) if val == 0]
                    one_neighbors_list_filtered = filter_by_zero_bits(one_neighbors_list_original, msb_indices, n)
                    derrida_coeffs.append(get_derrida_coefficient(in_edges, BF_list, one_neighbors_list_filtered))
                    net_sens.append(network_sensitivity(in_edges, BF_list))
                    final_hammings.append(final_Hamming_distance(in_edges,BF_list,t_0, t_inf, one_neighbors_list_filtered))
                    fragilities.append(Fragility(in_edges,BF_list,t_0, t_inf, one_neighbors_list_filtered))
                             
            append_rows_to_tsv(networks,models,func_types,derrida_coeffs, net_sens, final_hammings, fragilities, file_path, columns)
    


def batch(total_networks, total_cores):
    networks_per_core = total_networks // total_cores
    remaining_networks = total_networks % total_cores
    allocations = [networks_per_core] * total_cores
    allocations[:remaining_networks] = [x + 1 for x in allocations[:remaining_networks]]
    network_init_list = [sum(allocations[:i]) for i in range(len(allocations))]
    network_final_list = [sum(allocations[:i]) for i in range(1, len(allocations) + 1)]
    return network_init_list, network_final_list

if __name__ == '__main__':
    n = 12
    total_networks = 10
    total_cores = 2
    for k in range(2,3):
        t_0 = 500
        t_inf = 700
        func_types = ['scNCF', 'IMR', 'BMR']
        for func_type in func_types:
            model_per_net = 1
            net_init_list, net_final_list = batch(total_networks, total_cores)
            processes = []
            for net_init, net_final in zip(net_init_list, net_final_list):
                processes.append(multiprocessing.Process(target=Dataframe_of_three_measures,
                                                              args=(n,k,func_type,t_0,t_inf,net_init,net_final, model_per_net)))
            for process in processes:
                process.start()
            for process in processes:
                process.join()