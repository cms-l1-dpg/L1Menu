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
import re
import tdrstyle
from rootpy.interactive import wait
from matplotlib import pyplot as plt
from Config import DualMap, S1S2Map, S2S1Map

freq = 11245.6
nBunches = 2736
unit = "kHz"
# unit = "Hz"
# nBunches = 1
pubins = np.arange(0, 50, 0.2)
pumap = collections.defaultdict(list)
filedir ="Output/MenuPU/*_tsgv4_PU.csv"
s1csv = pd.read_csv("HLT_Fit_Run258425-260627_Tot10_fit.csv")
DrawStage1= False
mcdf = None

def GetStage1Fun(l1seed):
    s1seed = l1seed
    if l1seed in S2S1Map:
        s1seed = S2S1Map[l1seed]
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
    s1fun = ROOT.TF1(name, fun, 0, 50 )
    return s1fun

def DrawPU(f, l1seed, key=None):
    df = f[(f.L1Seed == l1seed )]
    RetVar = None

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
        graph.SetPoint(i, xx, yy)
        graph.SetPointError(i, 0, yee)

    c1 = ROOT.TCanvas("fd","Fdf", 600, 500)
    ROOT.gStyle.SetOptStat(000000000)
    tdrstyle.setTDRStyle()
    graph.Draw("AP")
    graph.SetTitle(l1seed)
    graph.GetXaxis().SetTitle("PileUp")
    graph.GetXaxis().SetLimits(0, 50)
    if unit == "Hz":
        graph.GetYaxis().SetTitle("Rate (nBunches = %d) [Hz]" % nBunches)
    if unit == "kHz":
        graph.GetYaxis().SetTitle("Rate (nBunches = %d) [kHz]" % nBunches)

    leg = ROOT.TLegend(0.7432886,0.1733615,0.9949664,0.3530655)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(62)
    leg.AddEntry(graph, "Data", "p")

    ## Get Stage1
    s1fun = GetStage1Fun(l1seed)
    if DrawStage1 and s1fun is not None:
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
    fun = "f2 = %.2f + %.2f*x + %.3f*x^2" % (f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2) )
    RetVar = f2.GetParameter(2) 
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
    f1.Draw("same")
    # if l1seed == "L1APhysics":
    graph.GetYaxis().SetRangeUser(0, f1.Eval(50))
    fun = "f1 = %.2f + %.2f*x " % (f1.GetParameter(0), f1.GetParameter(1))
    tex = ROOT.TLatex(0.19, 0.69, fun)
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.04)
    tex.SetTextColor(ROOT.kRed)
    tex.SetLineWidth(2)
    tex.Draw()


    if key is not None:
        tex = ROOT.TLatex(0.79, 0.89, key)
        tex.SetNDC()
        tex.SetTextFont(61)
        tex.SetTextSize(0.05)
        tex.SetTextColor(ROOT.kGreen+2)
        tex.SetLineWidth(2)
        tex.Draw()

    if mcdf is not None:
        DrawMCPoint(mcdf[ mcdf.L1Seed == l1seed], leg)

    c1.SetGridy()
    c1.SetGridx()
    leg.Draw()
    c1.Update()
    # c1.SaveAs("plots/PU_%s.root" % l1seed)
    if key is not None:
        c1.SaveAs("plots/PU_%s_%s.png" % (l1seed, key))
        c1.SaveAs("plots/PU_%s_%s.root" % (l1seed, key))
    else:
        c1.SaveAs("plots/PU_%s.png" % l1seed)
        c1.SaveAs("plots/PU_%s.root" % l1seed)
    # wait()
    return RetVar

def GetMCDataFrame(allfiles):
    global mcdf
    mcdf = pd.DataFrame()
    mclist = [ ]
    for file_ in allfiles:
        if "MC" in file_:
            m = re.match(".*MC(\d*)PU.*", file_)
            if m is not None:
                df_ = pd.read_csv(file_, index_col=None, header=0)
                df_['PileUp'] = m.group(1)
                mclist.append(df_)
    if len(mclist) > 0:
        mcdf = pd.concat(mclist)
        return mcdf
    else:
        return None

def DrawMCPoint(mcdf, leg=None):
    x = []
    y = []
    yerr = []
    for pu in mcdf.PileUp:
        x.append(float(pu))
        fired = mcdf[mcdf.PileUp==pu].Fired.sum()
        total = mcdf[mcdf.PileUp==pu].Total.sum()
        rate = float(fired)/ total * freq * nBunches
        raterr=  math.sqrt(float(fired))/ total * freq * nBunches
        if unit == "Hz":
            y.append(rate)
            yerr.append(raterr)
        if unit == "kHz":
            y.append(rate/ 1000)
            yerr.append(raterr / 1000)

    ## Draw the plot
    graph = ROOT.TGraphErrors(len(x))
    for i, (xx, yy, yee) in enumerate(zip(x, y, yerr)):
        # print i, xx, yy
        graph.SetPoint(i, xx, yy)
        graph.SetPointError(i, 0, yee)
    graph.SetMarkerStyle(32)
    graph.SetMarkerSize(1.8)
    graph.SetMarkerColor(ROOT.kRed)
    graph.Draw("P")
    if leg is not None:
        leg.AddEntry(graph, "MC", "p")

    fitlowx = min(x)
    fithighx= max(x)
    fitname = "pol2"
    graph.Fit("pol2", "QF", "", fitlowx, fithighx)
    f2 = graph.GetFunction(fitname).Clone()
    f2.SetLineColor(ROOT.kOrange+1)
    f2.SetLineWidth(2)
    fun = "MC = %.2f + %.2f*x + %.2f*x^2" % (f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2) )
    tex = ROOT.TLatex(0.19, 0.62, fun)
    f2.Draw("same")
    tex.SetNDC()
    tex.SetTextAlign(13)
    tex.SetTextFont(61)
    tex.SetTextSize(0.04)
    tex.SetTextColor(ROOT.kOrange+1)
    tex.SetLineWidth(2)
    tex.Draw()

def DrawPUperFile(filedir, l1seed, key=None):
    global mcdf
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
    mcdf = GetMCDataFrame(allfiles)

    if l1seed == "All":
        for seed in pd.unique(df.L1Seed):
            DrawPU(df, seed)
    else:
        return DrawPU(df, l1seed, key)

def DrawMuonPU():
    ermap = {
        # 0.8 : "MuER0p8",
        # 1.25: "MuER1p2",
        # 1.5: "MuER1p5",
        # 1.8: "MuER1p8",
        # 2.0: "MuER2p0",
        # 2.1: "MuER2p1",
        # 2.2: "MuER2p2",
        # 2.3: "MuER2p3",
        # 2.4: "MuER2p4",
        2.5: "MuER2p5",
    }
    l1seeds = [
        # "L1_SingleMu14er",
        # "L1_SingleMu16er",
        # "L1_SingleMu18er",
        "L1_SingleMu20er",
        # "L1_SingleMu22er",
        # "L1_SingleMu25er",
        # "L1_SingleMu30er",
    ]


    l1plts=[]
    for l1 in l1seeds:
        quardmap = dict()
        for k,v in ermap.iteritems():
            # filedir ="Output/MuonBXAll/r259626_tsgv4_%s_PU.csv" % v
            filedir ="Output/MuonBXAll/r*_tsgv4_%s_PU.csv" % v
            # filedir ="Output/MuonBX2//r*_tsgv4_%s_PU.csv" % v
            # filedir ="Output/TSGv4Study/MuonER/Menu_MuonStudy_r*_tsgv4_%s_PU.csv" % v
            quardmap[k] = DrawPUperFile(filedir, l1, v)
            # quardmap[k] = abs(DrawPUperFile(filedir, l1))
        # print quardmap
        quardmap = collections.OrderedDict(sorted(quardmap.items()))
        l1plts.append(plt.plot(quardmap.keys(), quardmap.values(), '-o', label=l1))
        # plt.xticks(range(len(quardmap)), quardmap.keys())
    plt.legend(l1seeds, loc="lower left")
    plt.xlabel("MuonER")
    plt.ylabel("Fitted Quad Term")
    # plt.legend(l1seeds, loc="upper left")
    plt.show()

if __name__ == "__main__":
    DrawPUperFile(filedir, "All")
    # DrawPUperFile(filedir, "L1APhysics")
