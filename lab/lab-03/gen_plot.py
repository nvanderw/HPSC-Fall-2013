import sys
import numpy as np
import matplotlib.pyplot as plt

def lines(f):
    line = f.readline()
    while line:
        yield line.strip()
        line = f.readline()

def line_to_pair(line):
    def productOfInts(seq):
        acc = 1
        for e in seq:
            acc *= int(e)
        return acc

    parts = line.split(".")
    while parts[0] != "matrix":
        del parts[0]

    n = productOfInts(parts[1].split("x"))
    time = float(line.split(",")[1])
    
    return (n, time)

def io_speed(results):
    lns = list(lines(results))
    read =  [line for line in lns if line.startswith("matrix")]
    write = [line for line in lns if line.startswith("output")]

    read_txt = [line for line in read if ".txt" in line]
    read_txt_matrix = np.array(map(line_to_pair, read_txt))
    read_hd5 = [line for line in read if ".h5" in line]
    read_hd5_matrix = np.array(map(line_to_pair, read_hd5))

    write_txt = [line for line in write if ".txt" in line]
    write_txt_matrix = np.array(map(line_to_pair, write_txt))
    write_hd5 = [line for line in write if ".h5" in line]
    write_hd5_matrix = np.array(map(line_to_pair, write_hd5))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title('Time to read/write matrices of N elements.')
    ax.set_xlabel('N (number of elements in matrix)')
    ax.set_ylabel('Time (seconds)')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(which="both")

    ax.plot(read_txt_matrix[:,0], read_txt_matrix[:,1], marker="o", label="Text (reading)")
    ax.plot(write_txt_matrix[:,0], write_txt_matrix[:,1], marker="o", label="Text (writing)")
    ax.plot(read_hd5_matrix[:,0], read_hd5_matrix[:,1], marker="o", label="HDF5 (reading)")
    ax.plot(write_hd5_matrix[:,0], write_hd5_matrix[:,1], marker="o", label="HDF5 (writing)")

    plt.legend(loc=4)
    plt.savefig("io_speed.png", dpi=400)

def matrix_speed(results):
    lns = [line.split(',') for line in lines(results)]
    naive = np.array([line[1:] for line in lns if line[0] == 'naive'])
    opt = np.array([line[1:] for line in lns if line[0] == 'opt'])
    fig = plt.figure()

    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('Time to multiply NxN square matrices')
    ax.set_xlabel('Square matrix dimension')
    ax.set_ylabel('Time to multiply (seconds)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(which="both")

    ax.plot(naive[:,0], naive[:,1], marker="o", label="Naive, wall-clock")
    ax.plot(naive[:,0], naive[:,2], marker="o", label="Naive, CPU time")

    ax.plot(opt[:,0], opt[:,1], marker="o", label="Optimized code, wall-clock")
    ax.plot(opt[:,0], opt[:,2], marker="o", label="Optimized code, CPU time")

    plt.legend(loc=4)
    plt.savefig('matrix_speed.png', dpi=400)

def cache(results):
    lns = [line.split(',') for line in lines(results)]
    naive = np.array([line[1:] for line in lns if line[0] == 'naive'])
    opt = np.array([line[1:] for line in lns if line[0] == 'opt'])
    fig = plt.figure()

    ax = fig.add_subplot(1, 1, 1)
    ax.set_title('L3 Cache Misses vs. Matrix Dimension')
    ax.set_xlabel('Square matrix dimension')
    ax.set_ylabel('L3 Cache Misses')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(which="both")

    ax.plot(naive[:,0], naive[:,1], marker="o", label="Naive routine")
    ax.plot(opt[:,0], opt[:,1], marker="o", label="Optimized code")

    plt.legend(loc=4)
    plt.savefig('matrix_cache.png', dpi=400)

def main():
    with open(sys.argv[1]) as results:
        # io_speed(results)
        # matrix_speed(results)
        cache(results)

if __name__ == '__main__':
    main()
