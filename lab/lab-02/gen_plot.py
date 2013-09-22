import numpy as np
import matplotlib.pyplot as plt
import sys

def main():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    data = np.loadtxt("matrix_data.txt", skiprows=1)
    (ns,reals,users) = (data[:,0], data[:,1], data[:,2])

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_title('Time to multiply matrices vs. number of rows/cols')
    ax.set_xlabel('N (number of rows/columns)')
    ax.set_ylabel('Time (seconds)')

    ax.plot(ns, reals, color='blue', marker='o', linestyle='dashed', label='Wall clock time')
    ax.plot(ns, users, color='red', marker='o', linestyle='dashed', label='Usermode CPU time')

    plt.legend()
    plt.savefig('partb.png', dpi=400)

if __name__ == '__main__':
    sys.exit(main())
