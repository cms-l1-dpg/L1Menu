#!/usr/bin/env python
# encoding: utf-8

# File        : CSV.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2016 Feb 12
#
# Description : 


import pandas as pd
import numpy as np
import glob
import math
import ROOT
import collections
import os

freq = 11245.6
nBunches = 2736
pubins = np.arange(0, 30, 0.2)
pumap = collections.defaultdict(list)

def DrawPU(f, l1seed):
    df = f[(f.L1Seed == l1seed )]
    # PileUp = pd.unique(df.PileUp)

    for i in range(0, len(pubins) -1):
        pumap[pubins[i]] = []
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i], df.PileUp <= pubins[i+1])].Fired.sum())
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i], df.PileUp <= pubins[i+1])].Total.sum())

    x = []
    y = []
    yerr = []
    for k, v in pumap.iteritems():
        x.append(k)
        if v[1] == 0:
            y.append(0)
            yerr.append(0)
        else:
            # print float(v[0])/v[1] * freq * nBunches
            y.append(float(v[0])/v[1] * freq * nBunches / 1000)
            yerr.append( math.sqrt(float(v[0]))/v[1] * freq * nBunches / 1000)

    glen = len([x_ for x_ in y if x_>0])
    ## Draw the plot
    graph = ROOT.TGraphErrors(glen)
    validx = [ i for i, e in zip(x, y) if e != 0]
    minx = min(validx)
    maxx = 31
    # maxx = max(validx)
    for i, (xx, yy, yee) in enumerate(zip(x, y, yerr)):
        # print i, xx, yy
        graph.SetPoint(i, xx, yy)
        graph.SetPointError(i, 0, yee)

    c1 = ROOT.TCanvas("fd","Fdf", 600, 500)
    # tdrstyle.setTDRStyle()
    graph.Draw("AP")
    graph.Fit("pol2", "QF", "", minx, maxx)
    # graph.Fit("pol2", minx, maxx)
    f = graph.GetFunction("pol2")
    f.SetLineColor(ROOT.kRed)
    graph.GetXaxis().SetTitle("PileUp")
    graph.GetYaxis().SetTitle("Rate (nBunches = %d) [kHz]" % nBunches)
    graph.SetTitle(l1seed)

    fun = "f = %.2f + %.2f*x + %.2f*x^2" % (f.GetParameter(0), f.GetParameter(1), f.GetParameter(2) )
    tex = ROOT.TLatex(0.15, 0.8, fun)
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.04)
    tex.SetLineWidth(2)
    tex.Draw()

    c1.SetGridy()
    c1.SetGridx()
    c1.Update()
    c1.SaveAs("plots/PU_%s.png" % l1seed)
    return


if __name__ == "__main__":
    allfiles = glob.glob("./*PU.csv")
    if not os.path.exists("plots"):
        os.mkdir("plots")

    df = pd.DataFrame()
    flist = [ ]
    for file_ in allfiles:
        df_ = pd.read_csv(file_, index_col=None, header=0)
        flist.append(df_)
    df = pd.concat(flist)

    # DrawPU(df, "L1APhysics")
    for seed in pd.unique(df.L1Seed):
        DrawPU(df, seed)
    exit()
