#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate the expected time of having a double-mutation cell
and the expected number of single-mutation cells at this time.
For three different mechanisms (models).

@author: yuewang
"""
import numpy as np
n = 100 # number of total cells, fixed in this process
model = 1 # 1 means different birth rates
    # 2 means different mutation rates
    # 3 means j can induce t, i.e., mut_jt > mut_ot
if model == 1:
    birthrate = [1.0, 2.0, 4.0, 4.0, 4.0] # birth rate for different types 
        # of cells: non-mutant, J-only, T-only, JT, TJ
    deathrate = [1.0, 1.0, 1.0, 1.0, 1.0] # death rate for different types 
        # of cells: non-mutant, J-only, T-only, JT, TJ
    mut_oj = 0.1 # mutation rate from non-mutant to j only
    mut_ot = 0.1 # mutation rate from non-mutant to t only
    mut_jt = 0.1 # mutation rate from j only to jt
    mut_tj = 0.1 # mutation rate from t only to tj
if model == 2:
    birthrate = [1.0, 2.0, 2.0, 2.0, 2.0] # birth rate for different types 
        # of cells: non-mutant, J-only, T-only, JT, TJ
    deathrate = [1.0, 1.0, 1.0, 1.0, 1.0] # death rate for different types 
        # of cells: non-mutant, J-only, T-only, JT, TJ
    mut_oj = 0.1 # mutation rate from non-mutant to j only
    mut_ot = 0.2 # mutation rate from non-mutant to t only
    mut_jt = 0.2 # mutation rate from j only to jt
    mut_tj = 0.1 # mutation rate from t only to tj
if model == 3:
    birthrate = [1.0, 2.0, 2.0, 2.0, 2.0] # birth rate for different types 
        # of cells: non-mutant, J-only, T-only, JT, TJ
    deathrate = [1.0, 1.0, 1.0, 1.0, 1.0] # death rate for different types 
        # of cells: non-mutant, J-only, T-only, JT, TJ
    mut_oj = 0.1 # mutation rate from non-mutant to j only
    mut_ot = 0.1 # mutation rate from non-mutant to t only
    mut_jt = 0.2 # mutation rate from j only to jt
    mut_tj = 0.1 # mutation rate from t only to tj
    
curr_dist = [[0.0] * (n + 1) for i in range(n+1)]
# probability of having i JAK2-only cells and j TET2-only cells at the 
# current time (no double-mutation cell yet)
curr_dist[0][0] = 1.0
j_prob = [0.0] * n 
# the probability of having i JAK2-only cells when 
# the first JAK2-TET2 cell appears, and no TET2-only cell exists
j_cumu_time = [0.0] * n 
# the time weighted by probability of having i JAK2-only cells when 
# the first JAK2-TET2 cell appears, and no TET2-only cell exists
t_prob = [0.0] * n
# the probability of having i TET2-only cells when 
# the first TET2-JAK2 cell appears, and no JAK2-only cell exists
t_cumu_time = [0.0] * n
# the time weighted by probability of having i TET2-only cells when 
# the first TET2-JAK2 cell appears, and no JAK2-only cell exists
remaining_prob = 1.0
# the probability that no double-mutation cell appears yet
failed_prob = 0.0 
# the probability that one double-mutation cell appears, but there are both
# JAK2-only cells and TET2-only cells
time = 0 # current time
while remaining_prob > 1e-8 and time <= 400:
    time += 1
    prev_dist = curr_dist
    temp_dist = [[0.0] * (n + 1) for i in range(n+1)] 
    # the probabilities of having i JAK2-only cells and j TET2-only cells 
    # after choosing one cell to die
    curr_dist = [[0.0] * (n + 1) for i in range(n+1)]
    for i in range(n+1):
        for j in range(n+1-i):
            if i > 0:
                temp_dist[i-1][j] += prev_dist[i][j] * (i * deathrate[1])\
                    / (i * deathrate[1] + j * deathrate[2] + (n - i - j) \
                       * deathrate[0])
            if j > 0:
                temp_dist[i][j-1] += prev_dist[i][j] * (j * deathrate[2])\
                    / (i * deathrate[1] + j * deathrate[2] + (n - i - j) \
                       * deathrate[0])
            temp_dist[i][j] += prev_dist[i][j] * ((n - i - j) * deathrate[0])\
                / (i * deathrate[1] + j * deathrate[2] + (n - i - j) \
                   * deathrate[0])
    # for cell death
    
    for i in range(n):
        for j in range(n-i):
            # if a JAK2-only cell divides
            curr_dist[i+1][j] += temp_dist[i][j] * (i * birthrate[1])\
                / (i * birthrate[1] + j * birthrate[2] + (n - 1 - i - j) \
                   * birthrate[0]) * (1 - mut_jt)
            lost_prob = temp_dist[i][j] * (i * birthrate[1])\
                / (i * birthrate[1] + j * birthrate[2] + (n - 1 - i - j) \
                   * birthrate[0]) * mut_jt 
                    # probability of having a JT cell
            if j == 0: # no TET2-only cell
                j_prob[i] += lost_prob
                j_cumu_time[i] += time * lost_prob
            else: # has TET2-only cell and the process fails
                failed_prob += lost_prob
            remaining_prob -= lost_prob
            
            # if a TET2-only cell divides
            curr_dist[i][j+1] += temp_dist[i][j] * (j * birthrate[2])\
                / (i * birthrate[1] + j * birthrate[2] + (n - 1 - i - j) \
                   * birthrate[0]) * (1 - mut_tj)
            lost_prob = temp_dist[i][j] * (j * birthrate[2])\
                / (i * birthrate[1] + j * birthrate[2] + (n - 1 - i - j) \
                   * birthrate[0]) * mut_tj
            if i == 0:
                t_prob[j] += lost_prob
                t_cumu_time[j] += time * lost_prob
            else:
                failed_prob += lost_prob
            remaining_prob -= lost_prob
            
            # if a wild type cell divides
            curr_dist[i][j] += temp_dist[i][j] * ((n - 1 - i - j) * birthrate[0])\
                / (i * birthrate[1] + j * birthrate[2] + (n - 1 - i - j) \
                   * birthrate[0]) * (1 - mut_oj - mut_ot)
            curr_dist[i+1][j] += temp_dist[i][j] * ((n - 1 - i - j) * birthrate[0])\
                / (i * birthrate[1] + j * birthrate[2] + (n - 1 - i - j) \
                   * birthrate[0]) * mut_oj
            curr_dist[i][j+1] += temp_dist[i][j] * ((n - 1 - i - j) * birthrate[0])\
                / (i * birthrate[1] + j * birthrate[2] + (n - 1 - i - j) \
                   * birthrate[0]) * mut_ot
    total_prob = np.sum(np.sum(curr_dist)) + sum(j_prob) + sum(t_prob) + failed_prob
    #print(total_prob)
print('model: ', model)
print('mean jt time ', sum(j_cumu_time) / sum(j_prob))
print('mean tj time ', sum(t_cumu_time) / sum(t_prob))
print('mean j only number ', sum([j_prob[i] * i for i in range(n)]) / sum(j_prob))
print('mean t only number ', sum([t_prob[i] * i for i in range(n)]) / sum(t_prob))







