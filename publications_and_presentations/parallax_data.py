import sys
import argparse

from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt


def main(infile):
    data = ascii.read(infile)
    x = data["X"]
    y = data["Y"]
    px = data["Parallax"]
    flux = data["Flux"]
    cflux = data["Cflux"]
    px = np.array([float(px[i]) for i in range(0, len(px)) if px[i] != "*"])
    x_ = np.array([float(x[i]) for i in range(0, len(x)) if x[i] != "*"])
    y_ = np.array([float(y[i]) for i in range(0, len(y)) if y[i] != "*"])
    flux = np.array([float(flux[i]) for i in range(0, len(flux)) if x[i] != "*"])
    cflux = np.array([float(cflux[i]) for i in range(0, len(cflux)) if y[i] != "*"])

    plt.figure(dpi=150)
    plt.scatter(x_, y_, s=flux)
    plt.xlabel("X [Kpc]")
    plt.ylabel("Y [Kpc]")
    plt.title("Flux")
    plt.show()

    plt.figure(dpi=150)
    plt.scatter(x_, y_, s=cflux)
    plt.xlabel("X [Kpc]")
    plt.ylabel("Y [Kpc]")
    plt.title("CFlux")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Visualize parallax data')
    parser.add_argument('inFILE', type=str, help='Specify the input csv file.')
    args = parser.parse_args()
    main(args.inFILE)
    sys.exit(0)
