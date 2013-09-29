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

def main():
    with open(sys.argv[1]) as results:
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

if __name__ == '__main__':
    main()
