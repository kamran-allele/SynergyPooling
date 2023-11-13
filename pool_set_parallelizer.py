"""
This script generates pool sets in parallel with multicoreprocessors
"""

import timeit
from multiprocessing import Pool, cpu_count, freeze_support

import pool_optimization_tools as pp


def set_looper(num_of_loops):
    set_scores = []

    for i in range(num_of_loops):
        performance_of_set = pp.pool_set_generator()
        set_scores.append(performance_of_set)

    return(set_scores)


if __name__ == '__main__':
    freeze_support()

    # adjust minutes_to_run (and/or the integer coefficient in loops_per_process) to change runtime
    minutes_to_run = 1  
    loops_per_process = 5 * minutes_to_run

    # set # of processors to run for generating pool sets in parallel 
    processors = cpu_count()
    # processors = 8

    loop_array = [loops_per_process] * processors

    start_time = timeit.default_timer()

    with Pool(processes=processors) as pooler:
        results = pooler.map(set_looper, loop_array)

        pooler.close()
        pooler.join()

        from itertools import chain
        flat_results = list(chain(*results))

    elapsed = timeit.default_timer() - start_time
    max_pool_efficiency = ((min(flat_results) - 263) / 263)
    sets_per_second = round(len(flat_results) / elapsed, 2)
    sets_per_second_processor = round((len(flat_results) / elapsed) / processors, 2)

    print("Number of Sets Generated: ", len(flat_results))
    print("Total Run Time: ", round(elapsed, 2), " s")
    print("Total Processors: ", processors, "\n")

    print("Speed: ", sets_per_second, " sets/second")
    print("Speed:  ", sets_per_second_processor, " sets/second per processor\n")

    print("Best Set:", min(flat_results), "  - has inefficiency: ", round(max_pool_efficiency, 2))
    print("Worst Set:", max(flat_results))

    import numpy as np
    print("mean: ", round(np.mean(flat_results),2)," std: ",round(np.std(flat_results),2)) 