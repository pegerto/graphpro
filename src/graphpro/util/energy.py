import numpy as np

_residue_order = ['CYS','PHE','LEU','TRP','VAL','ILE','MET','HIS','TYR','ALA',
            'GLY','PRO','ASN','THR','SER','ARG','GLN','ASP','LYS','GLU']

_bt_data = [
    [1.34, -0.53, -0.5, -0.74, -0.51, -0.48, -0.49, -0.19, -0.16, -0.26, -0.09, -0.18, 0.28, 0, 0.09, 0.32, 0.04, 0.38, 0.35, 0.46],
    [-0.53, -0.82, -0.78, -0.78, -0.67, -0.65, -0.89, -0.19, -0.49, -0.33, 0.11, -0.19, 0.29, 0, 0.1, 0.08, -0.04, 0.48, 0.11, 0.34],
    [-0.5, -0.78, -0.81, -0.7, -0.8, -0.79, -0.68, 0.1, -0.44, -0.37, 0.14, -0.08, 0.36, 0, 0.26, 0.09, 0.08, 0.62, 0.16, 0.37],
    [-0.74, -0.78, -0.7, -0.74, -0.62, -0.65, -0.94, -0.46, -0.55, -0.4, -0.24, -0.73, -0.09, 0, 0.07, -0.41, -0.11, 0.06, -0.28, -0.15],
    [-0.51, -0.67, -0.8, -0.62, -0.72, -0.68, -0.47, 0.18, -0.27, -0.38, 0.04, -0.08, 0.39, 0, 0.25, 0.17, 0.17, 0.66, 0.16, 0.41],
    [-0.48, -0.65, -0.79, -0.65, -0.68, -0.6, -0.6, 0.19, -0.33, -0.35, 0.21, 0.05, 0.55, 0, 0.35, 0.18, 0.14, 0.54, 0.21, 0.38],
    [-0.49, -0.89, -0.68, -0.94, -0.47, -0.6, -0.56, -0.17, -0.51, -0.23, 0.08, -0.16, 0.32, 0, 0.32, 0.17, -0.01, 0.62, 0.22, 0.24],
    [-0.19, -0.19, 0.1, -0.46, 0.18, 0.19, -0.17, -0.33, -0.21, 0.21, 0.23, -0.05, 0.1, 0, 0.15, 0.04, 0.22, -0.22, 0.26, -0.11],
    [-0.16, -0.49, -0.44, -0.55, -0.27, -0.33, -0.51, -0.21, -0.27, -0.15, -0.04, -0.4, 0.01, 0, 0.07, -0.37, -0.18, -0.07, -0.4, -0.16],
    [-0.26, -0.33, -0.37, -0.4, -0.38, -0.35, -0.23, 0.21, -0.15, -0.2, -0.03, 0.07, 0.24, 0, 0.15, 0.27, 0.21, 0.3, 0.2, 0.43],
    [-0.09, 0.11, 0.14, -0.24, 0.04, 0.21, 0.08, 0.23, -0.04, -0.03, -0.2, -0.01, 0.1, 0, 0.1, 0.14, 0.2, 0.17, 0.12, 0.48],
    [-0.18, -0.19, -0.08, -0.73, -0.08, 0.05, -0.16, -0.05, -0.4, 0.07, -0.01, -0.07, 0.13, 0, 0.17, -0.02, -0.05, 0.25, 0.12, 0.26],
    [0.28, 0.29, 0.36, -0.09, 0.39, 0.55, 0.32, 0.1, 0.01, 0.24, 0.1, 0.13, -0.04, 0, 0.14, 0.02, -0.05, -0.12, -0.14, -0.01],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0.09, 0.1, 0.26, 0.07, 0.25, 0.35, 0.32, 0.15, 0.07, 0.15, 0.1, 0.17, 0.14, 0, 0.13, 0.12, 0.25, 0.01, 0.1, 0.1],
    [0.32, 0.08, 0.09, -0.41, 0.17, 0.18, 0.17, 0.04, -0.37, 0.27, 0.14, -0.02, 0.02, 0, 0.12, 0.13, -0.12, -0.71, 0.5, -0.75],
    [0.04, -0.04, 0.08, -0.11, 0.17, 0.14, -0.01, 0.22, -0.18, 0.21, 0.2, -0.05, -0.05, 0, 0.25, -0.12, 0.14, 0.12, -0.2, 0.1],
    [0.38, 0.48, 0.62, 0.06, 0.66, 0.54, 0.62, -0.22, -0.07, 0.3, 0.17, 0.25, -0.12, 0, 0.01, -0.71, 0.12, 0.27, -0.69, 0.4],
    [0.35, 0.11, 0.16, -0.28, 0.16, 0.21, 0.22, 0.26, -0.4, 0.2, 0.12, 0.12, -0.14, 0, 0.1, 0.5, -0.2, -0.69, 0.38, -0.87],
    [0.46, 0.34, 0.37, -0.15, 0.41, 0.38, 0.24, -0.11, -0.16, 0.43, 0.48, 0.26, -0.01, 0, 0.1, -0.75, 0.1, 0.4, -0.87, 0.45]
]

BT_potential = {}

for i, res1 in enumerate(_residue_order):
    for j, res2 in enumerate(_residue_order):
        BT_potential[(res1, res2)] = _bt_data[i][j]

def bt_potential(res1, res2):
    def normalize_res_name(res_name):
        if res_name in ('HSD', 'HSE', 'HSP'):
            return 'HIS'
        return res_name
    return BT_potential[(normalize_res_name(res1), normalize_res_name(res2))]
    
    
def compute_bt_potential(atom_group, chain, cutoff=6, epsilon=1):
    from scipy.spatial import distance
    from numpy import linalg as LA
    
    ca_position = atom_group.c_alphas_positions(chain)
    residues = atom_group.c_alphas_residues(chain)
    dist = distance.squareform(distance.pdist(ca_position))
    potential = np.zeros((len(dist), len(dist)))
    res_ids = [res['resid'] for res in residues]
    
    for i in range(len(dist)):
        for j in range(i + 1, len(dist)):
            resname_i = residues[i]['resname']
            resname_j = residues[j]['resname']
            
            V_ij = bt_potential(resname_i, resname_j)
            r_ij = dist[i,j]
            
            if r_ij < cutoff:
                # Lennard-Jones weight on distance
                energy = epsilon * V_ij * ((cutoff / r_ij) ** 6 - (cutoff / r_ij) ** 12)
                potential[i,j] = energy
                potential[j,i] = energy
    
    eigen_value, _ = LA.eig(potential)
    return res_ids, eigen_value


def compute_eigen_centrality(atom_group, chain, cutoff=6, max_iter=500):
    from scipy.spatial import distance
    import networkx as nx
    
    ca_position = atom_group.c_alphas_positions(chain)
    residues = atom_group.c_alphas_residues(chain)
    dist = distance.squareform(distance.pdist(ca_position))
    potential = np.zeros((len(dist), len(dist)))
    res_ids = [res['resid'] for res in residues]
    
    for i in range(len(dist)):
        for j in range(i + 1, len(dist)):
            resname_i = residues[i]['resname']
            resname_j = residues[j]['resname']
            
            V_ij = bt_potential(resname_i, resname_j)
            r_ij = dist[i,j]
            
            if r_ij < cutoff:
                energy = np.exp(-V_ij)
                potential[i,j] = energy
                potential[j,i] = energy
    G = nx.from_numpy_array(potential)
    eigen_value = nx.eigenvector_centrality(G,  max_iter=max_iter, weight='weight')
    return res_ids, eigen_value