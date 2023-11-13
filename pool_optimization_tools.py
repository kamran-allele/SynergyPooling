from random import *
import numpy as np
from scipy.special import comb
import itertools


## For compound library of specific size, will generate a pool set solution (5 compounds per pool)
def pool_set_generator():    
    
    ## change compound library size here
    Lib_Size = 73

    combinations = {}
    singles_count = [0] * Lib_Size    
    [combinations, singles_count] = init_pools(Lib_Size)
    pool_set_out = all_pools(Lib_Size,combinations,singles_count)
    num_of_pools_in_set = len(pool_set_out)

    return(num_of_pools_in_set)


######################################################################

# function for getting max values of dictionary.
def keys_withmaxval(dict_to_max):
    values = np.array(list(dict_to_max.values()))
    k = np.array(list(dict_to_max.keys()))

    searchval = max(values)
    ii = np.where(values == searchval)[0]
    return(list(k[ii]))


def keys_withminval(dict_to_min):
    values = np.array(list(dict_to_min.values()))
    k = np.array(list(dict_to_min.keys()))

    searchval = min(values)
    ii = np.where(values == searchval)[0]
    return(list(k[ii]))


# scorer(combos_to_score):
# Takes in a list of combinations to score (each element in list is an array representing 1 combination)
# Determines number of new combinations in list
# Returns that number (or score)
def scorer(combos_to_score,combinations_s):
    score = len([1 for p in combos_to_score if combinations_s[p] == 0])
    penalty = 0

    ## penalty formula here is proportional to rep_freq squared 
    #penalty = sum([-(combinations_s[p]**2) for p in combos_to_score])

    ## penalty formula here is proportional to rep_freq
    #penalty = sum([-(combinations_s[p]) for p in combos_to_score])
    
    score = score + penalty

    return score

# comb_generator(a,n):
# Takes in a pool (a) and a new compound (n),
# Returns a list of combinations generated when
# adding the new compound to pool (each comb in list is an array)
def comb_generator(a, n):
    comb_gen = [tuple(sorted((i, n))) for i in a]
    return comb_gen


def add_comb_scorer(pool_current, compounds_to_score,combinations_a):
    scorer2 = scorer
    new_comb_scores1 = {}
    new_comb_scores1 = {i: scorer2(comb_generator(pool_current, i),combinations_a) for i in compounds_to_score}
    return new_comb_scores1

# generates an optimal pool given a starting compound seed
def pool_gen2(seed,Lib_Size,combinations_p,singles_count_p):
    pool = [seed]
    new_pool_score = 0
    
    while len(pool) < 5:

        Lib_List = list(range(Lib_Size))
        new_comb_scores = {}

        for i in pool:
            if i in Lib_List:
                Lib_List.remove(i)

        new_comb_scores = add_comb_scorer(pool, Lib_List,combinations_p)

        max_combs_set = keys_withmaxval(new_comb_scores)
        min_single = {}
        for i in max_combs_set:
            min_single[i] = singles_count_p[i]
        new_member2 = sample(keys_withminval(min_single), 1)

        set1 = comb_generator(pool, new_member2[0])
        for d in set1:
            if combinations_p[d] == 0:
                singles_count_p[d[0]] = singles_count_p[d[0]] + 1
                singles_count_p[d[1]] = singles_count_p[d[1]] + 1
                new_pool_score = new_pool_score + 1
            combinations_p[d] = combinations_p[d] + 1

        pool.append(new_member2[0])
    
    return(pool,combinations_p,singles_count_p,new_pool_score)


def set_to_mat(pool_s):
    p_c = np.zeros((len(pool_s),73))
    
    row_index =0
    
    for p in pool_s:
        for c in p:
            p_c[row_index][c] = 1
        row_index = row_index + 1
        
    return(p_c)

def prune_mat(p_c):
    
    prev_mat = p_c[0:p_c.shape[0]-2] 
    final_vector = p_c[p_c.shape[0]-1].transpose()
    pool_overlap = prev_mat @ final_vector
    
    pool_index_to_eliminate = pool_overlap[0:round(len(pool_overlap)*.9)-1].argmax()

    return(pool_index_to_eliminate)


def seed_selection(singles_count_ss):
    
    if sum(singles_count_ss) == 0:
        #sub_max_list_ss = [0] 
        sub_max_list_ss = [ind for ind, sing in enumerate(singles_count_ss) if sing <= min(singles_count_ss)]
    
    else:
        if len(singles_count_ss) == singles_count_ss.count(singles_count_ss[0]):
            sub_max_list_ss = [ind for ind, sing in enumerate(singles_count_ss) if sing > 0]

        else:
            sub_max_list_ss = [ind for ind, sing in enumerate(singles_count_ss) if sing <= min(singles_count_ss)]
    
    return(sub_max_list_ss)


def remove_pool_from_mem(pool_to_del,combinations,singles_count):
    
    pair_sets = list(itertools.combinations(pool_to_del,2))
                
    for pair_drop in pair_sets:
        s_p_d = tuple(sorted(pair_drop))
        combinations[s_p_d] = combinations[s_p_d] - 1
                    
        if combinations[s_p_d] == 0:
            singles_count[pair_drop[0]] = singles_count[pair_drop[0]] - 1
            singles_count[pair_drop[1]] = singles_count[pair_drop[1]] - 1
    
    return()


def all_pools(Lib_Size,combinations,singles_count):
    pool_list = []
    total_combos_possible = comb(Lib_Size, 2)
    pool_gen3 = pool_gen2
    old_p_score = 0
    new_p_score = 0
    deleted_count = 0
    del_mem = 0
    
    while (sum(singles_count) / 2) < total_combos_possible:
        pop = []
        sub_max_list = seed_selection(singles_count)
        start = choice(sub_max_list)
        
        seed_db = {}
        if len(sub_max_list) < 65 and len(pool_list)>200:
            for seed_choice in sub_max_list:
                [pop,combinations,singles_count,new_p_score] = pool_gen3(seed_choice,Lib_Size,combinations,singles_count)
            
                ## add overlap evaluation here as well
                seed_db[seed_choice] = new_p_score
                remove_pool_from_mem(pop,combinations,singles_count)
        
            start = sample(keys_withmaxval(seed_db),1)[0]
        
        
        [pop,combinations,singles_count,new_p_score] = pool_gen3(start,Lib_Size,combinations,singles_count)
                
        pool_list.append(pop)
        
        if len(pool_list)>80:
            #print(del_mem,deleted_count)
            
            if del_mem !=0:
                del_mem = del_mem + 1
            if del_mem >2:
                del_mem = 0
                        
            if deleted_count < 600 and del_mem == 0:

                del_mem = del_mem + 1
                deleted_count = deleted_count + 1
                
                mat_set = set_to_mat(pool_list)
                removed_pool_ind = prune_mat(mat_set)
                                
                ex_pool = pool_list.pop(removed_pool_ind)
                remove_pool_from_mem(ex_pool,combinations,singles_count)

                                    
    return(pool_list)

def init_pools(Lib_Size):
    combinations1 = {}
    singles_count1 = [0] * Lib_Size
    i = 0

    # creates a dictionary (combinations) of all possible combinations and how often they've appeared
    while i < Lib_Size:
        j = i + 1
        while j < Lib_Size:
            pair_set = (i, j)
            combinations1[pair_set] = 0
            j = j + 1
        i = i + 1
    return(combinations1, singles_count1)

######

# Gives stats on 1 pool set generation run (pool size of 5) 
if __name__ == "__main__":
    
    total_num_of_compounds_to_pool = 73

    import timeit
    start_time = timeit.default_timer()
    pool_set_size_out = pool_set_generator()
    elapsed = timeit.default_timer() - start_time

    perfect_pooling = (total_num_of_compounds_to_pool) * (total_num_of_compounds_to_pool - 1) / (2 * 10)

    for _ in [1]:
        p_waste = round(pool_set_size_out / (perfect_pooling) - 1, 3)
        print("Run Time: ", round(elapsed, 2), " s")
        print("Total number of unique combinations represented")
        print("Total pools: ", pool_set_size_out)
        print("pooling waste:", p_waste)
