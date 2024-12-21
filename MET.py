import numpy as np
import UniFrac
import preprocessing
import math
from collections import defaultdict
import mosek
from mosek.fusion import Model, Domain, Expr, ObjectiveSense
from mosek.fusion import OptimizeError, SolutionError, ProblemStatus

"""
Returns the condition for P_e and Q_e to be close to each other
"""
def p_condition(P_e, Q_e, n, c1=1.8):
    diff = np.abs(P_e - Q_e)
    cond1 = np.sqrt(1.1*c1*(P_e + Q_e) * np.ln(n) / n)
    cond2 = np.abs(P_e + Q_e)
    return diff < np.min(cond1, cond2)

"""
Given distributions P and Q, define sets E_r and E_c such that:
- for e in E_r, |P_e - Q_e| is close to 0
- for e in E_r, |P_e - Q_e| is not zero with high probability 

Parameters:
P : a list of proportions of subtree under edge e (length |E|)
G : a list of proportions of subtree under edge e (length |E|)
n : number of sequence reads in expectation
c1 : constant, default

Returns:
E_r : list of indices of edges in E_r
E_c : list of indices of edges in E_c
"""
def split_er_ec(P, Q, n, c1=1.8):
    E_r = []
    E_c = []
    for i in range(len(P)):
        P_e = P[i]
        Q_e = Q[i]
        if p_condition(P_e, Q_e, n, c1):
            E_r.append(i)
        else:
            E_c.append(i)
    return np.array(E_r), np.array(E_c)

def partition_er(P_er, Q_er, n, c1=1.8):
    J = np.floor(np.log2(n/(c1 * np.log(n))))
    part_ranges = []
    sums = P_er + Q_er
    for j in range(1, J + 1):
        upper = 2**(-j-1)
        lower = 2**(-j)
        partition = np.where((sums > lower) & (sums <= upper))[0]
        part_ranges.append(partition)
    E_0 = [np.where(sums <= 2**(-J))[0]]  
    return E_0 + part_ranges

def unbiased_H(P, k, n):
    ls = np.array([P - m / n for m in range(0, k)])
    return np.prod(ls) 

def unbiased_G(P, Q, k, n):
    g = np.zeros(k + 1)
    for l in range(k + 1):
        g[l] = math.comb(k, l) * (-1)**l * unbiased_H(P, l, n) * unbiased_H(Q, k - l, n)
    return np.sum(g)

def decompose_disjoint_paths(E_w, E_j):
    """
    Decomposes a subset of edges into disjoint paths.

    Parameters:
    - E_w (list of tuples): Subset of edges to be decomposed (u, v).
    - E_j (list of tuples): Subset of E_w that determine final paths

    Returns:
    -  S_j (int): Number of disjoint paths using E_w that contain at least one edge in E_j
    """
    # Build the adjacency list for the subset of edges
    adj = defaultdict(list)
    for u, v in E_w:
        adj[u].append(v)
    
    visited_edges = set()
    disjoint_paths = [] # paths will be list of edges
    # Function to find paths using DFS
    def find_path(node, path, path_verts):
        for neighbor in adj[node]:
            edge = (node, neighbor)
            if edge not in visited_edges and node not in path_verts:
                #path_verts.append(neighbor)
                path_verts.append(node)
                visited_edges.add(edge)
                path.append(edge)
                path, path_verts = find_path(neighbor, path, path_verts)
        return path, path_verts
       
    for u, v in E_w:
        edge = (u, v)
        if edge not in visited_edges:
            visited_edges.add(edge)
            path_u, pu = find_path(u, [], path_verts=[])
            path_v, pv = find_path(v, [], path_verts=[])
            if v not in pu:
                disjoint_paths.append(path_u)
                disjoint_paths.append([edge] + path_v)
            else:
                path = path_u + [edge] + path_v
                disjoint_paths.append(path)
            
    required = {(u, v) for u, v in E_j}
    S_j = sum([1 for path in disjoint_paths if set(path).intersection(required)])
    return S_j

def deviation_R(j, Sj, k, M, d, n, c1=1.8):
    """
    Constraint upper bound for deviation of estimator from G for E1..EJ
    """
    C = 6 * d * M**2 (48*c1)**(k/2)
    R = C * d * np.sqrt(Sj * np.log(n))*(np.log(n) / (2**j * n))**(k/2)
    return R

def deviation_R0(k1, k2, M, d, n, c1=1.8):
    """
    Constraint upper bound for deviation of estimator from G for E0
    """
    C = 4* d * M**2 * (76 * c1)**(k1 + k2)
    R0 = C * d * np.sqrt(n*np.log(n)**2) * (np.log(n) / n)**(k1 + k2)

def get_depth(tree):
    # TODO
    return

def moment_screening_estimator(tree, P, Q, L, n, c1=1.8, c2=1):
    """
    Given a tree, empirical distributions P and Q, and edge lengths L, compute the MET of the Wasserstein distance
    
    Parameters:
    -----------
    tree : Phylo.BaseTree
        Phylogenetic tree
    P : dictionary over edges (indexed by int)
        Relative abundances under each edge of sample P
    Q : dictionary over edges (indexed by int)
        Relative abundances under each edge of sample Q
    L : dictionary over edges (indexed by int)
        Lengths of each edge
    n : int
        Number of sequence reads in expectation
    """
    d = get_depth(tree)
    J = np.floor(np.log2(n/(c1 * np.log(n))))
    # split edges into E_r and E_c
    inds_er, inds_ec = split_er_ec(P, Q, n, c1)
    # partition E_r into E_0, E_1, ..., E_J
    # taje a list of P and Q?
    partitions = partition_er(P[inds_er], Q[inds_ec], n, c1)
    # to compute Sj partition all edged into Ew
    Ew_partitions = partition_er(P, Q, n, c1)
    # compute sj for each partition
    Sjs = []
    for part in partitions:
        # get the paired version of edged
        Sj = decompose_disjoint_paths(inds_er[part], inds_ec)
        Sjs.append(Sj)
    # compute the MET
    w0 = 1/(2*n)
    h1h2_range = np.arange(1, np.ceil(c1 * np.log(n) / (w0 * n)))
    E0 = partitions[0]
    W0 = {}
    P_e0 = P[E0]
    Q_e0 = Q[E0]
    L_e0 = L[E0]
    for h1 in h1h2_range:
        for h2 in h1h2_range:
            H_h1h2 = [
                np.arange((h1 - 1) * w0, h1 * w0), 
                np.arange((h2 - 1) * w0, h2 * w0)
            ]
 
            e_inds = np.where(np.isin(P_e0, H_h1h2[0]) & np.isin(Q_e0, H_h1h2[1]))[0]
            W0[(h1, h2)] = np.sum(L_e0[e_inds])

    R0 = {}
    K = math.floor(c2 * np.log(n))
    for k1 in range(K + 1):
        for k2 in range(K + 1):
            R0[(k1, k2)] = deviation_R0(k1, k2, np.max(L), d, n, c1)

    h_inds = [(h1, h2) for h1 in h1h2_range for h2 in h1h2_range]
    k_inds = [(k1, k2) for k1 in range(K + 1) for k2 in range(K + 1)]           
    with Model("LP") as M:
        # dec vars
        W = M.variable("W", len(h_inds), Domain.unbounded())
        Z = M.variable("Z", len(h_inds), Domain.greaterThan(0.0))
        U = M.variable("U", len(k_inds), Domain.greaterThan(0.0))
        M.objective("MinimizeDeviation", ObjectiveSense.Minimize, Expr.sum(Z))

        # add constraints
        # |W - W_0| <= Z -> to linearize abs val
        for i, (h1, h2) in enumerate(h_inds):
            M.constraint(Expr.sub(W.index(i), W0[(h1, h2)]), Domain.lessThan(Z.index(i)))
            M.constraint(Expr.sub(W0[(h1, h2)], W.index(i)), Domain.lessThan(Z.index(i)))
        # constraints on error terms
        for i, (k1, k2) in enumerate(k_inds):
            lhs = Expr.dot(h_inds, [(h1 * w0)** k1 * (h2 *w0) ** k2 for (h1, h2) in h_inds])
            rhs = sum(L[e] * unbiased_H(P[e], k1, n) * unbiased_H(Q[e], k2, n) for e in E0)
            M.constraint(Expr.sub(lhs, rhs), Domain.lessThan(U.index(i)))
            M.constraint(Expr.sub(rhs, lhs), Domain.lessThan(U.index(i)))
            M.constraint(U.index(i), Domain.lessThan(R0[(k1, k2)]))
        try:

            M.solve()
            M.acceptedSolutionStatus(mosek.fusion.AccSolutionStatus.Optimal)
            
            W_opt = W.level()
            Z_opt = Z.level()
            U_opt = U.level()
            W_0_opt = {h: val for h, val in zip(h_inds, W_opt)}
            Z_dict = {h: val for h, val in zip(h_inds, Z_opt)}
            U_dict = {k: val for k, val in zip(k_inds, U_opt)}

        except mosek.Error as e: 
            print("Response code {0}\nMessage       {1}".format(e.errno, e.msg))
            W_0_opt = W0
        
    # compute D0
    D0 = np.zeros(len(h_inds))
    for i (h1, h2) in enumerate(h_inds):
        D0[i] = W_0_opt[(h1, h2)] * np.abs((h1 - h2) * w0)
    D0 = np.sum(D0)

    for j in range(1, J + 1):
        part_inds = partitions[j]
        ej_inds = inds_er[part_inds]
        Sj = Sjs[j]
        Rj = [deviation_R(j, Sj, k, np.max(L), d, n, c1) for k in range(K + 1)]
        wj = 1/np.sqrt(2**j * n)
        
        with Model("LP") as M:
            # dec vars
            W = M.variable("W", len(h_inds), Domain.unbounded())
            Z = M.variable("Z", len(h_inds), Domain.greaterThan(0.0))
            U = M.variable("U", len(k_inds), Domain.greaterThan(0.0))
            M.objective("MinimizeDeviation", ObjectiveSense.Minimize, Expr.sum(Z))

            # add constraints
            # |W - W_0| <= Z -> to linearize abs val
            for i, (h1, h2) in enumerate(h_inds):
                M.constraint(Expr.sub(W.index(i), W0[(h1, h2)]), Domain.lessThan(Z.index(i)))
                M.constraint(Expr.sub(W0[(h1, h2)], W.index(i)), Domain.lessThan(Z.index(i)))
            # constraints on error terms
            for i, (k1, k2) in enumerate(k_inds):
                lhs = Expr.dot(h_inds, [(h1 * w0)** k1 * (h2 *w0) ** k2 for (h1, h2) in h_inds])
                rhs = sum(L[e] * unbiased_H(P[e], k1, n) * unbiased_H(Q[e], k2, n) for e in ej_inds)
                M.constraint(Expr.sub(lhs, rhs), Domain.lessThan(U.index(i)))
                M.constraint(Expr.sub(rhs, lhs), Domain.lessThan(U.index(i)))
                M.constraint(U.index(i), Domain.lessThan(Rj[(k1, k2)]))
            try:

                M.solve()
                M.acceptedSolutionStatus(mosek.fusion.AccSolutionStatus.Optimal)
                
                W_opt = W.level()
                Z_opt = Z.level()
                U_opt = U.level()
                W_0_opt = {h: val for h, val in zip(h_inds, W_opt)}
                Z_dict = {h: val for h, val in zip(h_inds, Z_opt)}
                U_dict = {k: val for k, val in zip(k_inds, U_opt)}

            except mosek.Error as e: 
                print("Response code {0}\nMessage       {1}".format(e.errno, e.msg))
                W_0_opt = W0


    

