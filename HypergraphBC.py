#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:44:57 2021

@author: Hugh Warden
"""

import networkx as nx
import numpy as np
import scipy.sparse as ss
import scipy.sparse.linalg as ssalg
from scipy.sparse import csr_matrix


def get_hyperedge_connectivity_adj_mat(inc_mat):
    m = inc_mat.shape[0]
    n = inc_mat.shape[1]
    adj_mat = np.zeros((n,n))
    for i in range(m):
        for j in range(n-1):
            for k in range(j+1, n):
                if (inc_mat[i,j] == 1) & (inc_mat[i,k] == 1):
                    adj_mat[j,k] = 1
                    adj_mat[k,j] = 1
    return adj_mat

def get_hyperedge_connectivity_adj_mat_fast(inc_mat):
    inc_mat = csr_matrix(inc_mat)
    tinc_mat = ss.csr_matrix.transpose(inc_mat)
    ret_mat = tinc_mat @ inc_mat
    return np.matrix(ret_mat)

def get_hyperedge_connectivity_adj_mat_fast_2(inc_mat):
    tinc_mat = np.transpose(inc_mat)
    ret_mat = tinc_mat @ inc_mat
    return np.matrix(ret_mat)

#TODO test for disconnected hypergraphs
def get_hypergraph_betweenness(adj_mat, inc_mat):
    g = nx.convert_matrix.from_numpy_matrix(adj_mat)
    n_vert = adj_mat.shape[0]
    m_vert = inc_mat.shape[0]
    bc_all = np.zeros(m_vert)
    bc_nt = np.zeros(m_vert)
    bc_act = np.zeros(m_vert)
    
    for i in range(m_vert-1):
        for j in range(i+1, m_vert):
            temp_bc_all = np.zeros(m_vert)
            temp_bc_nt = np.zeros(m_vert)
            temp_bc_act = np.zeros(m_vert)
            #OOP problem, fix having to call this every iteration
            temp = g.copy()
            temp.add_node(n_vert+1)
            for k in range(n_vert):
                if inc_mat[i,k] == 1:
                    temp.add_edge(n_vert+1, k)
            temp.add_node(n_vert+2)
            for k in range(n_vert):
                if inc_mat[j,k] == 1:
                    temp.add_edge(k,n_vert+2)
            
            count = 0
            for p in nx.all_shortest_paths(temp, n_vert+1, n_vert+2):
                temp_inc_mat = inc_mat[:,p[1:(len(p)-1)]]
                r_sum = np.squeeze(np.asarray(temp_inc_mat.sum(1)))
                temp_bc_all[np.where(r_sum > 0)] += 1
                temp_bc_act[np.where(r_sum > 1)] += 1
                if len(p) > 4:
                    temp_inc_mat = inc_mat[:,p[2:(len(p)-2)]]
                    r_sum = np.squeeze(np.asarray(temp_inc_mat.sum(1)))
                    temp_bc_nt[np.where(r_sum > 0)] += 1
                else:
                    temp_bc_nt[np.where(r_sum > 0)] += 1
                    
                count += 1
            
            if count > 0:
                temp_bc_all = temp_bc_all / count
                temp_bc_nt = temp_bc_nt / count
                temp_bc_act = temp_bc_act / count
                
                temp_bc_all[i] = 0
                temp_bc_all[j] = 0
                temp_bc_nt[i] = 0
                temp_bc_nt[j] = 0
                temp_bc_act[i] = 0
                temp_bc_act[j] = 0
            
                bc_all += temp_bc_all
                bc_nt += temp_bc_nt
                bc_act += temp_bc_act
            
    return bc_all, bc_nt, bc_act

'''
inc_mat1 = np.matrix([[1,0,0,0],[1,1,0,0],[0,1,1,0],[0,1,1,0],[0,0,1,1],[0,0,0,1]])
adj_mat1 = np.matrix([[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]])

print(get_hyperedge_connectivity_adj_mat_fast(inc_mat1))
print(adj_mat1)

print(get_hypergraph_betweenness(adj_mat1, inc_mat1))
print(get_hypergraph_betweenness(get_hyperedge_connectivity_adj_mat(inc_mat1), inc_mat1))
print(get_hypergraph_betweenness(get_hyperedge_connectivity_adj_mat_fast(inc_mat1), inc_mat1))

inc_mat2 = np.matrix([[1,0,0,0],[1,1,0,0],[0,1,1,0],[0,1,1,0],[0,0,1,1],[0,0,0,1],[0,0,1,0]])
adj_mat2 = np.matrix([[0,1,0,0],[1,0,1,0],[0,1,0,1],[0,0,1,0]])

print(get_hypergraph_betweenness(adj_mat2, inc_mat2))
'''
