#!/usr/bin/env python
# encoding: utf-8

# File        : CompHLT.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2016 Aug 25
#
# Description :

import pandas as pd
import numpy as np
import glob
import math
import ROOT
import collections
import os
import re
import tdrstyle
import rootpy
from rootpy.interactive import wait
from rootpy.io import root_open
from matplotlib import pyplot as plt
from Config import DualMap, S1S2Map, S2S1Map

foldername = "Sep07fill_7118_and_7131_nanoDST_Prescale_2018_v2_1_0_Col_2.0_HATS"


fit_min = 22
fit_max = 55
plot_min = 0		#min of x axis
plot_max = 80		#max of x axis
minx = 20		#minx of points
maxx = 80		#maxx of points
maxy = 8
miny = 0
set_logy = False
if(set_logy): miny = 0.01


freq = 11245.6
nBunches = 2544
unit = "kHz"

#fitname = "pol4"
#fitname = "expo"
fitname = ROOT.TF1("fitname","[0]*x + [1]*x*x",0,80);

filedir = "/eos/uscms/store/user/huiwang/L1Menu2017/" + foldername + "/*Default_PU.csv"

pubins = np.arange(minx,maxx, 1)
pumap = collections.defaultdict(list)

PatMap = {  
#    "L1APhysics" : "L1APhysics",
   "L1_DoubleJet_110_35_DoubleJet35_Mass_Min620" : "L1_DoubleJet_110_35.*_DoubleJet35.*_Mass_Min620",
#    "L1_SingleJet_FWD3p0" : "L1_SingleJet\d+_FWD3p0",
#    "L1_ETM120" : "L1_ETM.*120",
#    "SingleMu22" : "L1_SingleMu22",
}

def DrawPU(canvas, f, l1seed, count, key=None):
    df = f[(f.L1Seed == l1seed )]
    RetVar = None

    for i in range(0, len(pubins) -1):
        pumap[pubins[i]] = []
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i], df.PileUp <= pubins[i+1])].Fired0.sum())
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i], df.PileUp <= pubins[i+1])].Total.sum())

    x = []
    y = []
    yerr = []
    for k, v in pumap.iteritems():
        if v[1] != 0:
            x.append(k)
            if unit == "Hz":
                y.append(float(v[0])/v[1] * freq * nBunches )
                yerr.append( math.sqrt(float(v[0]))/v[1] * freq * nBunches )
            if unit == "kHz":
                y.append(float(v[0])/v[1] * freq * nBunches / 1000)
                yerr.append( math.sqrt(float(v[0]))/v[1] * freq * nBunches / 1000)

    ## Draw the plot
    graph = ROOT.TGraphErrors(len(x))

    for i, (xx, yy, yee) in enumerate(zip(x, y, yerr)):
        graph.SetPoint(i, xx, yy)
        graph.SetPointError(i, 0, yee)

    graph.SetMarkerStyle(20+count)
    graph.SetMarkerSize(1.5)
    graph.SetMarkerColor(1+count)
    graph.SetLineColor(1+count)
    graph.SetLineWidth(2)
    tdrstyle.setTDRStyle()
    canvas.cd()
    canvas.Update()
    if count == 0:
        graph.Draw("AP")
        graph.GetXaxis().SetTitle("PileUp")
        graph.GetXaxis().SetLimits(plot_min, plot_max)
        graph.GetYaxis().SetRangeUser(miny, maxy)
        graph.GetYaxis().SetTitle("Rate (nBunches = %d) [%s]" % (nBunches, unit))
    else:
        graph.Draw("P")
    canvas.Update()
    leg.AddEntry(graph, l1seed, "p")

    result_ptr = graph.Fit(fitname, "SQ", "", fit_min, fit_max)
    f2 = graph.GetFunction("fitname").Clone()
    #f2 = graph.GetFunction(fitname).Clone()
    f2.SetLineColor(1+count)
    f2.SetLineWidth(2)
    f2.SetRange(plot_min, fit_min)
    f2.SetLineStyle(5)
    minChi = f2.GetChisquare() / f2.GetNDF()
    #fun = "Fit = %.2f + %.2f*x + %.3f*x^2" % (f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2) )
    fun = "Fit = %.3f*x + %.4f*x^2" % (f2.GetParameter(0), f2.GetParameter(1) )
    #fun = "Fit = exp(%.2f + %.3f*x)" % (f2.GetParameter(0), f2.GetParameter(1) )
    print "%s, Chi2/NDF = %.2f" % (fun, minChi)
    f2.Draw("same")
    #graph.Write()
    f2_2 = f2.Clone("dashline2")
    f2_2.SetRange(fit_max, plot_max)
    f2_2.Draw("same")

    if key is not None:
        tex = ROOT.TLatex(0.2, 0.9, key)
        tex.SetNDC()
        tex.SetTextFont(61)
        tex.SetTextSize(0.045)
        tex.SetTextColor(ROOT.kGreen+2)
        tex.SetLineWidth(2)
        tex.Draw()

    canvas.Update()


def DrawL1(key, pattern):
    c1.Clear()
    leg.Clear()

    inputlist = []
    pat = re.compile('^%s$' % pattern)

    for x in [x for x in pd.unique(df.L1Seed)]:
        if pat.match(x):
            inputlist.append(x)
    print key,

    for i, seed in enumerate(inputlist):
        DrawPU(c1, df, seed, i)
    leg.Draw()

    c1.SetGrid()

    box = ROOT.TBox(10, 8, 70, 12)
    box.SetFillColor(38)
    box.SetFillStyle(3002)

    if(set_logy): c1.SetLogy()
    c1.Update()
    #c1.SaveAs("plots/%s_%s.root" % (key, foldername))
    #c1.SaveAs("plots/%s_%s.C" % (key, foldername))
    c1.SaveAs("plots/%s_%s.png" % (key, foldername))
    #c1.SaveAs("plots/%s_%s.pdf" % (key, foldername))

if __name__ == "__main__":
    allfiles = glob.glob(filedir)
    if not os.path.exists("plots"):
        os.mkdir("plots")

    df = pd.DataFrame()
    flist = [ ]
    for file_ in allfiles:
        df_ = pd.read_csv(file_, index_col=None, header=0)
        flist.append(df_)
    df = pd.concat(flist)

    ## Redefine PatMap for each L1Seed in the dataframe
    #PatMap = {k:k for k in pd.unique(df.L1Seed)}

    ROOT.gStyle.SetOptStat(000000000)
    tdrstyle.setTDRStyle()
    c1 = ROOT.TCanvas("fd","Fdf", 800, 800)
    leg = ROOT.TLegend(0.15,0.6,0.45,0.85)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(62)
    leg.SetTextSize(0.05)
    for k, v in PatMap.items():
        DrawL1(k, v)
        # wait()
