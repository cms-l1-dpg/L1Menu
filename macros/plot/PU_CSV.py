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
import tdrstyle
from rootpy.interactive import wait
from Config import DualMap, S1S2Map, S2S1Map

freq = 11245.6
# nBunches = 1
nBunches = 2736
unit = "kHz"
# unit = "Hz"
pubins = np.arange(0, 30, 0.2)
pumap = collections.defaultdict(list)
filedir = "./PUresult/PUresult/*PU.csv"
s1csv = pd.read_csv("HLT_Fit_Run258425-260627_Tot10_fit.csv")

def GetStage1Fun(l1seed):
    s1seed = l1seed
    if l1seed in S2S1Map:
        s1seed = S2S1Map[l1seed]
    print l1seed, s1seed
    s1df = s1csv.loc[s1csv['path name'] == s1seed]
    if len(s1df) == 0:
        return None
    # if s1df['fit function'] != 'quad':
        # return None
    x0 = s1df.X0
    x1 = s1df.X1
    x2 = s1df.X2
    name = "S1 = %.2f + %.2f*x + %.2f*x^2" % (x0, x1, x2)
    fun = "%f + %f*x + %f*x*x" % (x0, x1, x2)
    # fun = s1df.X0 + s1
    s1fun = ROOT.TF1(name, fun, 0, 35 )
    return s1fun





def DrawPU(f, l1seed):
    df = f[(f.L1Seed == l1seed )]

    for i in range(0, len(pubins) -1):
        pumap[pubins[i]] = []
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i], df.PileUp <= pubins[i+1])].Fired.sum())
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i], df.PileUp <= pubins[i+1])].Total.sum())

    # # No merging
    # PileUp = pd.unique(df.PileUp)
    # for i in PileUp:
        # pumap[i] = []
        # pumap[i].append(df[df.PileUp == i].Fired.sum())
        # pumap[i].append(df[df.PileUp == i].Total.sum())

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
    minx = min(x)
    maxx = 31
    for i, (xx, yy, yee) in enumerate(zip(x, y, yerr)):
        # print i, xx, yy
        graph.SetPoint(i, xx, yy)
        graph.SetPointError(i, 0, yee)

    c1 = ROOT.TCanvas("fd","Fdf", 600, 500)
    ROOT.gStyle.SetOptStat(000000000)
    # tdrstyle.setTDRStyle()
    graph.Draw("AP")
    graph.SetTitle(l1seed)
    graph.GetXaxis().SetTitle("PileUp")
    graph.GetXaxis().SetLimits(0, 35)
    if unit == "Hz":
        graph.GetYaxis().SetTitle("Rate (nBunches = %d) [Hz]" % nBunches)
    if unit == "kHz":
        graph.GetYaxis().SetTitle("Rate (nBunches = %d) [kHz]" % nBunches)
 
    ## Get Stage1
    s1fun = GetStage1Fun(l1seed)
    if s1fun is not None:
        s1fun.SetLineColor(ROOT.kGreen+2)
        s1fun.SetLineWidth(2)
        s1fun.Draw("same")
        tex = ROOT.TLatex(0.19, 0.81, s1fun.GetName())
        tex.SetNDC()
        tex.SetTextAlign(13)
        tex.SetTextFont(61)
        tex.SetTextSize(0.04)
        tex.SetTextColor(ROOT.kGreen+2)
        tex.SetLineWidth(2)
        tex.Draw()


    tex = ROOT.TLatex(0.19, 0.9, l1seed)
    # f2.Draw()
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.05)
    tex.SetTextColor(ROOT.kBlack)
    tex.SetLineWidth(2)
    tex.Draw()
    ## Pol2
    fitname = "pol2"
    graph.Fit("pol2", "QF", "", minx, maxx)
    f2 = graph.GetFunction(fitname).Clone()
    f2.SetLineColor(ROOT.kBlue)
    f2.SetLineWidth(2)
    fun = "f2 = %.2f + %.2f*x + %.2f*x^2" % (f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2) )
    tex = ROOT.TLatex(0.19, 0.75, fun)
    f2.Draw("same")
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.04)
    tex.SetTextColor(ROOT.kBlue)
    tex.SetLineWidth(2)
    tex.Draw()

    ## Pol1
    fitname = "pol1"
    graph.Fit("pol1", "QF", "", minx, maxx)
    f1 = graph.GetFunction(fitname)
    f1.SetLineColor(ROOT.kRed)
    f1.SetLineWidth(2)
    # f1.Draw()
    fun = "f1 = %.2f + %.2f*x " % (f1.GetParameter(0), f1.GetParameter(1))
    tex = ROOT.TLatex(0.19, 0.69, fun)
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.04)
    tex.SetTextColor(ROOT.kRed)
    tex.SetLineWidth(2)
    tex.Draw()


    c1.SetGridy()
    c1.SetGridx()
    c1.Update()
    # c1.SaveAs("plots/PU_%s.root" % l1seed)
    c1.SaveAs("plots/PU_%s.png" % l1seed)
    # wait()
    return


if __name__ == "__main__":
    allfiles = glob.glob(filedir)
    if not os.path.exists("plots"):
        os.mkdir("plots")

    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetOptStat(000000000)
    df = pd.DataFrame()
    flist = [ ]
    for file_ in allfiles:
        df_ = pd.read_csv(file_, index_col=None, header=0)
        flist.append(df_)
    df = pd.concat(flist)

    DrawPU(df, "L1APhysics")
    # for seed in pd.unique(df.L1Seed):
        # DrawPU(df, seed)
    exit()
