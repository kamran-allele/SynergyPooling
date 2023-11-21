"""
This script generates pool sets in parallel with multicoreprocessors
"""

# Specificy Compound Library Size
compound_library_size = 73

# Adjust minutes_to_run (and/or the integer coefficient in loops_per_process) to change runtime
minutes_to_run = 1  
loops_per_process = 5 * minutes_to_run

# Set num of processors to run for parallelization of pool set generation
from multiprocessing import cpu_count
processors = cpu_count()


##############################################

import timeit
from multiprocessing import Pool, freeze_support
import pool_optimization_tools as pp


def set_looper(num_of_loops):
    set_scores = []

    for i in range(num_of_loops):
        pool_solution, performance_of_set = pp.pool_set_generator(compound_library_size)
        set_scores.append(performance_of_set)

    return(set_scores)



if __name__ == '__main__':
    freeze_support()

    loop_array = [loops_per_process] * processors

    start_time = timeit.default_timer()

    with Pool(processes=processors) as pooler:
        results = pooler.map(set_looper, loop_array)

        pooler.close()
        pooler.join()

        from itertools import chain
        flat_results = list(chain(*results))

    elapsed = timeit.default_timer() - start_time

    print("Number of Sets Generated: ", len(flat_results))
    print("Total Run Time: ", round(elapsed, 2), " s")

    print("Best Set:", min(flat_results))
    print("Worst Set:", max(flat_results))

    import numpy as np
    print("mean: ", round(np.mean(flat_results),2)) 