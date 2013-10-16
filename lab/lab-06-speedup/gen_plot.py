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
    speedup_matrix = get_speedup(results_matrix)
    
    def _inner():
        for (nthreads, speedup) in speedup_matrix:
            yield (nthreads, speedup / nthreads)

    return np.array(list(_inner()))

def get_karp_flatt(results_matrix):
    speedup_matrix = get_speedup(results_matrix)

    def _inner():
        for (nthreads, speedup) in speedup_matrix:
            if nthreads > 1:
                yield (nthreads, (1/speedup - 1/nthreads) / (1 - 1 / nthreads))
    return np.array(list(_inner()))

def main():
    jacobi_averaged_matrix, multiply_averaged_matrix = None, None
    with open("jacobi_results.txt", "r") as jacobi_file:
        jacobi_results = (parse_results(line.rstrip().split(',')) for line in lines(jacobi_file))
        jacobi_averaged_matrix = np.array(list(average_results(jacobi_results, 5)))

    with open("matrix_results.txt", "r") as matrix_file:
        multiply_results = (parse_results(line.rstrip().split(',')) for line in lines(matrix_file))
        multiply_averaged_matrix = np.array(list(average_results(multiply_results, 5)))

    jacobi_speedup = get_speedup(jacobi_averaged_matrix)
    jacobi_efficiency = get_efficiency(jacobi_averaged_matrix)
    jacobi_kf = get_karp_flatt(jacobi_averaged_matrix)

    multiply_speedup = get_speedup(multiply_averaged_matrix)
    multiply_efficiency = get_efficiency(multiply_averaged_matrix)
    multiply_kf = get_karp_flatt(multiply_averaged_matrix)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title("Speedup vs. Number of threads")
    ax.set_xlabel("Number of threads")
    ax.set_ylabel("Speedup")

    ax.plot(jacobi_speedup[:,0], jacobi_speedup[:,1], marker="o", label="Jacobi iteration")
    ax.plot(multiply_speedup[:,0], multiply_speedup[:,1], marker="o", label="Matrix multiply")

    plt.legend(loc=2)
    plt.savefig("speedup.png", dpi=400)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title("Efficiency vs. Number of threads")
    ax.set_xlabel("Number of threads")
    ax.set_ylabel("Efficiency")

    ax.plot(jacobi_efficiency[:,0], jacobi_efficiency[:,1], marker="o", label="Jacobi iteration")
    ax.plot(multiply_efficiency[:,0], multiply_efficiency[:,1], marker="o", label="Matrix multiply")

    plt.legend()
    plt.savefig("efficiency.png", dpi=400)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title("Experimental serial fraction vs. Number of threads")
    ax.set_xlabel("Number of threads")
    ax.set_ylabel("Serial fraction")

    ax.plot(jacobi_kf[:,0], jacobi_kf[:,1], marker="o", label="Jacobi iteration")
    ax.plot(multiply_kf[:,0], multiply_kf[:,1], marker="o", label="Matrix multiply")

    plt.legend()
    plt.savefig("kf.png", dpi=400)

if __name__ == '__main__':
    main()
