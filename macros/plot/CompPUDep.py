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


label = "ZeroBias"
PU = 35
freq = 11245.6
config = 2017
filedir = "Nov192017EGTT28_v2/*CaloTower_PU.csv"
maxy = 40
if config == 2016:
    nBunches = 2208
    unit = "kHz"
if config == 2017:
    nBunches = 2592
    unit = "kHz"

if config == 1:
    nBunches = 1
    unit = "Hz"

pubins = np.arange(11,60, 1)
pumap = collections.defaultdict(list)

PatMap = {  
    # "EG" :  "L1_StrategyEG[34]\d+",
    "ETM" : "L1_ETM\d+",
    # "HTT" : "L1_HTT[23]\d+",
    # "DEGHTT" : "L1_DoubleEG6_HTT.*",
    # "dPhi80" : "L1_ETM\d+_Jet80_dPhi_Min0p4",
    # "dPhi60" : "L1_ETM\d+_Jet60_dPhi_Min0p4",
    # "dETM80" : "L1_ETM80_Jet\d+_dPhi_Min0p4",
    # "dETM100" : "L1_ETM100_Jet\d+_dPhi_Min0p4",
    # 'Mu6HTT': 'L1_Mu6_HTT2\d+',
    # 'Mu8HTT': 'L1_Mu8_HTT2\d+',
    # 'Mu5EG': 'L1_Mu5_EG2\d+',
    # 'Mu6EG': 'L1_Mu6_EG2\d+',
    # 'DMu': 'L1_DoubleMu0er.*_dEta_Max1p8_OS',
    # 'DEG12': 'L1_DoubleEG_\d+_12',
    # 'DEG25': 'L1_DoubleEG_25_\d+',
    # 'DTau': 'L1_DoubleIsoTau\d+er',
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
    minx = min(x)
    maxx = max(x)

    for i, (xx, yy, yee) in enumerate(zip(x, y, yerr)):
        # if yy != 0 and yee/yy >0.3:
            # continue
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
        graph.GetXaxis().SetLimits(10, 70)
        graph.GetYaxis().SetRangeUser(0, maxy)
        graph.GetYaxis().SetTitle("Rate (nBunches = %d) [%s]" % (nBunches, unit))
    else:
        graph.Draw("P")
    canvas.Update()
    leg.AddEntry(graph, l1seed, "p")

    ## Pol2
    fitname = "pol2"
    graph.Fit(fitname, "Q", "", minx, max(x))
    f2 = graph.GetFunction(fitname).Clone()
    f2.SetLineColor(1+count)
    f2.SetLineWidth(2)
    fun = "Fit = %.2f + %.2f*x + %.3f*x^2" % (f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2) )
    # f2.Draw("same")
    f2_2 = f2.Clone("dashline2")
    f2_2.SetRange(minx, 70)
    f2_2.SetLineStyle(5)
    # f2_2.Draw("same")
    key = "Rate(PU=%d): %.1f" % (PU, f2_2.Eval(PU))
    tex = ROOT.TLatex(0.15, 0.75, fun)
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.04)
    tex.SetTextColor(ROOT.kBlue)
    tex.SetLineWidth(2)
    # tex.Draw()


    if key is not None:
        tex = ROOT.TLatex(0.55, 0.85, key)
        tex.SetNDC()
        tex.SetTextFont(61)
        tex.SetTextSize(0.045)
        tex.SetTextColor(ROOT.kGreen+2)
        tex.SetLineWidth(2)
        # tex.Draw()

    # leg.Draw()
    canvas.Update()


def DrawL1(key, pattern):
    c1.Clear()
    leg.Clear()

    inputlist = []
    pat = re.compile('^%s$' % pattern)

    for x in [x for x in pd.unique(df.L1Seed)]:
        if pat.match(x):
            inputlist.append(x)
    print key, " : ",  inputlist

    for i, seed in enumerate(inputlist):
        DrawPU(c1, df, seed, i)
    leg.Draw()

    if config == 2016:
        l37 = ROOT.TLine(37, 0, 37, maxy)
        l37.SetLineColor(2)
        l37.SetLineWidth(2)
        l37.Draw()
        l40 = ROOT.TLine(40, 0, 40, maxy)
        l40.SetLineColor(2)
        l40.SetLineWidth(2)
        l40.Draw()
        l47 = ROOT.TLine(47, 0, 47, maxy)
        l47.SetLineColor(2)
        l47.SetLineWidth(2)
        l47.Draw()
        l52 = ROOT.TLine(52, 0, 52, maxy)
        l52.SetLineColor(2)
        l52.SetLineWidth(2)
        l52.Draw()

    if config == 2017:
        l47 = ROOT.TLine(47, 0, 47, maxy)
        l47.SetLineColor(2)
        l47.SetLineWidth(2)
        l47.Draw()
        l55 = ROOT.TLine(55, 0, 55, maxy)
        l55.SetLineColor(2)
        l55.SetLineWidth(2)
        l55.Draw()
        l60 = ROOT.TLine(60, 0, 60, maxy)
        l60.SetLineColor(2)
        l60.SetLineColor(2)
        l60.SetLineWidth(2)
        l60.Draw()

    tex = ROOT.TLatex(0.2, 0.3, "%d Thresholds" % config)
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.04)
    tex.SetTextColor(ROOT.kBlue)
    tex.SetLineWidth(2)
    # tex.Draw()
    c1.SetGrid()


    box = ROOT.TBox(10, 8, 70, 12)
    box.SetFillColor(38)
    box.SetFillStyle(3002)

    c1.Update()
    c1.SaveAs("plots/%s_%d.root" % (key, config))
    c1.SaveAs("plots/%s_%d.png" % (key, config))
    c1.SaveAs("plots/%s_%d.pdf" % (key, config))

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
    # PatMap = {k:k for k in pd.unique(df.L1Seed)}

    ROOT.gStyle.SetOptStat(000000000)
    tdrstyle.setTDRStyle()
    c1 = ROOT.TCanvas("fd","Fdf", 1200, 1000)
    leg = ROOT.TLegend(0.1711409,0.5412262,0.4714765,0.9238901)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(62)
    for k, v in PatMap.items():
        DrawL1(k, v)
        # wait()
