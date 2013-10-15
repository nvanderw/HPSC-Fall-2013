from __future__ import division
from matplotlib import pyplot as plt
import numpy as np

def lines(handle):
    line = handle.readline()
    while line:
        yield line
        line = handle.readline()

def average_results(results, ntrials):
    d = {}
    for (threads, trial, val) in results:
        if threads in d:
            d[threads] += val
        else:
            d[threads] = val

    for (threads, val) in d.iteritems():
        yield (threads, val / ntrials)


def parse_results(seq):
    return [int(seq[0]), int(seq[1]), float(seq[2])]

def get_speedup(results_matrix):
    single_threaded = results_matrix[0, 1]

    def _inner():
        for (nthreads, time) in results_matrix:
            yield (nthreads, single_threaded / time)
    return np.array(list(_inner()))

def get_efficiency(results_matrix):
    speedup = get_speedup(results_matrix)
    
    def _inner():
        for (nthreads, speedup) in results_matrix:
            yield (nthreads, speedup / nthreads)

    return np.array(list(_inner()))

def main():
    with open("jacobi_results.txt", "r") as jacobi_file:
        jacobi_results = (parse_results(line.rstrip().split(',')) for line in lines(jacobi_file))
        jacobi_averaged_matrix = np.array(list(average_results(jacobi_results, 5)))
        jacobi_speedup = get_speedup(jacobi_averaged_matrix)
        jacobi_efficiency = get_efficiency(jacobi_averaged_matrix)

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(jacobi_speedup[:,0], jacobi_speedup[:,1], marker="o")
        plt.show()

if __name__ == '__main__':
    main()
