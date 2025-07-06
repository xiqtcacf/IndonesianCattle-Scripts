import math
from pathlib import Path
import argparse
from math import isinf

def add_ext(path, *args, sep=".", t=Path):
    p = str(path)
    for x in args:
        p = p + sep + x
    return t(p)

## DIRECT COPY FROM msmc-tools / plot_utils.py
class MSMCresult:
    """
    Simple class to read a MSMC result file. Constructor takes filename of MSMC result file
    """
    def __init__(self, filename):
        f = open(filename, "r")
        self.times_left = []
        self.times_right = []
        self.lambdas = []
        next(f) # skip header
        for line in f:
            fields = line.strip().split()
            nr_lambdas = len(fields) - 3
            if len(self.lambdas) == 0:
                self.lambdas = [[] for i in range(nr_lambdas)]
            time_left = float(fields[1])
            time_right = float(fields[2])
            self.times_left.append(time_left)
            self.times_right.append(time_right)
            for i in range(nr_lambdas):
                l = float(fields[3 + i])
                self.lambdas[i].append(l)
        self.T = len(self.times_left)
        # self.times_left[0] = self.times_left[1] / 4.0
        # self.times_right[-1] = self.times_right[-2] * 4.0

    def getInterval(self, t):
        for i, (tl, tr) in enumerate(zip(self.times_left, self.times_right)):
            if tl <= t and tr > t:
                return i
        else:
            raise ValueError("Should not happen!")

    def getLambdaAt(self, t, lambda_index=0):
        i = self.getInterval(t)
        return self.lambdas[lambda_index][i]

def popSizeStepPlot(filename, mu=1.25e-8, gen=30.0):
    """
    to be used with a step-plot function, e.g. matplotlib.pyplot.steps.
    returns (x, y), where x contains the left point of each step-segment in years, and y contains the effective population size. Note that there are two ways to make a step-plot. You should make sure that your step-plot routine moves right and then up/down instead of the other way around.
    If plotted on a logarithmic x-axis, you should adjust x[0] = x[1] / 4.0, otherwise the leftmost segment will start at 0 and won't be plotted on a log-axis.

    Options:
        mu: Mutation rate per generation per basepair (default=1.25e-8)
        gen: generation time in years (default=30)
    """
    M = MSMCresult(filename)
    x = [t * gen / mu for t in M.times_left]
    y = [(1.0 / l) / (2.0 * mu) for l in M.lambdas[0]]
    return (x, y)


def combinecoal(fh,withinCoalFile1, withinCoalFile2,crossCoalFile):
    print("time_index", "left_time_boundary", "right_time_boundary", "lambda_00", "lambda_01", "lambda_11", sep="\t", file=fh)
    cc = MSMCresult(crossCoalFile)
    within1 = MSMCresult(withinCoalFile1)
    within2 = MSMCresult(withinCoalFile2)
    res = 10
    for i, (tl, tr, l) in enumerate(zip(cc.times_left, cc.times_right, cc.lambdas[0])):
        lambda00 = 0.0
        lambda11 = 0.0
        if isinf(tr):
            lambda00 = within1.lambdas[0][-1]
            lambda11 = within2.lambdas[0][-1]
        else:
            for j in range(res):
                t = tl + j / float(res) * (tr - tl)
                lambda00 += within1.getLambdaAt(t) / float(res)
                lambda11 += within2.getLambdaAt(t) / float(res)
        print(i, tl, tr, lambda00, l, lambda11, sep="\t", file=fh)


def crossCoalPlotCombined(filename1, filename2, filename12, mu=1.25e-8, gen=30.0):
    """
    returns (x, y) where x is the time in years and y is the relative cross coalescence rate.
    This function should be used with msmc2, where the three files correspond to the three different
    cases within-population1, within-population2 and across populations.
    Check also the doc-string in popSizeStepPlot for Options and hints on how to plot.
    """
    M1 = MSMCresult(filename1)
    M2 = MSMCresult(filename2)
    M12 = MSMCresult(filename12)
    # I1 = M1.getInterp()
    # I2 = M2.getInterp()
    x = []
    y = []
    resolution = 10
    for i in range(M12.T - 1):
        tLeft = M12.times_left[i]
        tRight = M12.times_right[i]
        avgWithinRate = 0.0
        for j in range(resolution):
            t = tLeft + j / float(resolution) * (tRight - tLeft)
            lambda1 = M1.getLambdaAt(t)
            lambda2 = M2.getLambdaAt(t)
            avgWithinRate += 0.5 * (lambda1 + lambda2) / float(resolution)
        x.append(tLeft * gen / mu)
        y.append(M12.lambdas[0][i] / avgWithinRate)
    return (x, y)
