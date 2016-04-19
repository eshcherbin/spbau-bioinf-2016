#!/usr/bin/python
from numpy import *
from matplotlib.pyplot import *
import sys


def main():
    if len(sys.argv) == 1:
        print('no file specified')
        sys.exit(0)

    with open(sys.argv[1]) as f:
        ref_name = f.readline().strip()
        max_cvg, avg, cvg = map(float, f.readline().split())
        max_cvg = int(max_cvg)
        x, y = zip(*map(lambda s: tuple(map(float, s.split())), f.readlines()))

    figure()

    title(ref_name)
    plot(x, y)
    grid(True)

    subplots_adjust(top=0.8, bottom=0.1)

    figtext(0.1, 0.95, 'maximal coverage: %d reads' % max_cvg)
    figtext(0.1, 0.9, 'average coverage: %.2f reads/position' % avg)
    figtext(0.1, 0.85, '%.2f%% positions covered' % (cvg * 100.))

    savefig('plot_' + ref_name + '.png')


if __name__ == '__main__':
    main()
