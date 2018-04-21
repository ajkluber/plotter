import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('xdata', type=str, help='Filename of X data.')
    parser.add_argument('ydata', type=str, help='Filename of Y data.')
    parser.add_argument('xlabel', type=str, help='X axis label.')
    parser.add_argument('ylabel', type=str, help='Y axis label.')
    parser.add_argument('title', type=str, help='Title.')
    parser.add_argument('saveas', type=str, default=None, help='Name to save plot.')
    parser.add_argument('--nodisplay', action='store_true', help='Optional. Dont display to screen.')
    args = parser.parse_args()

    xdata = args.xdata
    ydata = args.ydata
    xlabel = args.xlabel
    ylabel = args.ylabel
    title = args.title
    saveas = args.saveas
    nodisplay = args.nodisplay

    if nodisplay:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    else:
        import matplotlib.pyplot as plt

    if xdata.endswith(".npy"):
        x = np.load(xdata)
    else:
        x = np.loadtxt(xdata)

    if ydata.endswith(".npy"):
        y = np.load(ydata)
    else:
        y = np.loadtxt(ydata)

    plt.plot(x, y)
    plt.xlabel(xlabel,fontsize=20)
    plt.ylabel(ylabel,fontsize=20)
    plt.title(title)

    plt.savefig(saveas + ".pdf")
    plt.savefig(saveas + ".png")
