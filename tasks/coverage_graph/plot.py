#!/usr/bin/python

from numpy import *
from matplotlib.pyplot import *

with open('plot.dat') as f:
    max_cvg, avg, cvg = map(float, f.readline().split())
    max_cvg = int(max_cvg)
    x, y = zip(*map(lambda s: tuple(map(float, s.split())), f.readlines()))

figure()

plot(x, y)
grid(True)

subplots_adjust(top=0.8, bottom=0.1)

figtext(0.1, 0.95, 'maximal coverage: %d reads' % max_cvg)
figtext(0.1, 0.9, 'average coverage: %.2f reads/position' % avg)
figtext(0.1, 0.85, '%.2f%% positions covered' % (cvg * 100.))

savefig('plot.png')
