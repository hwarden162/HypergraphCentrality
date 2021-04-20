#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 09:44:57 2021

@author: Hugh Warden
"""

import networkx as nx
import numpy as np

def get_hyperedge_dual_adj_mat(inc_mat):
    tinc_mat = np.transpose(inc_mat)
    ret_mat = tinc_mat @ inc_mat
    return np.matrix(ret_mat)

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

