#!/usr/bin/env python
# encoding: utf-8

# File        : CompDir.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2016 Dec 09
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
#import CMS_lumi
from rootpy.interactive import wait
from rootpy.io import root_open
from matplotlib import pyplot as plt
from Config import DualMap, S1S2Map, S2S1Map

#ROOT.gROOT.ForceStyle()

plot_min = 31
plot_max = 55
fit_min = 0.
fit_max = 60.
ymax = 180

def GetPUGraph(f, l1seed):
    df = pd.DataFrame()
    flist = [ ]
    for file_ in f:
        df_ = pd.read_csv(file_, index_col=None, header=0)
        flist.append(df_)
    df = pd.concat(flist)
    df = df[(df.L1Seed == l1seed )]


    for i in range(0, len(pubins) -1):
        pumap[pubins[i]] = []
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i],
                                                  df.PileUp <=
                                                  pubins[i+1])].Fired0.sum())
        pumap[pubins[i]].append(df[np.logical_and(df.PileUp > pubins[i], df.PileUp <= pubins[i+1])].Total.sum())

    x = []
    y = []
    tot = []
    yerr = []
    for k, v in pumap.iteritems():
        if v[1] != 0:
            x.append(k)
            tot.append(v[1])
	    if unit == "Hz":
                y.append(float(v[0])/v[1] * freq * nBunches )
                yerr.append( math.sqrt(float(v[0]))/v[1] * freq * nBunches )
	    if unit == "kHz":
                y.append(float(v[0])/v[1] * freq * nBunches /1000 )
                yerr.append( math.sqrt(float(v[0]))/v[1] * freq * nBunches /1000 )


    plt.plot(x, tot)

    plt.xlabel('PU')
    plt.ylabel('count')
    plt.savefig("MC.png")
    ## Draw the plot
    graph = ROOT.TGraphErrors(len(x))
    for i, (xx, yy, yee) in enumerate(zip(x, y, yerr)):
        # if yy != 0 and yee/yy >0.3:
            # continue
	if xx < plot_min or xx > plot_max:
		continue
        graph.SetPoint(i, xx, yy)
        graph.SetPointError(i, 0, yee)

    graph.SetLineWidth(2)
    graph.SetTitle(l1seed)
    return graph

def FittingGraphs(graph, fittype =""):
    ## Pol2
    pol2_c0 = ROOT.TF1("pol2_c0","[0]*x + [1]*x*x",0,80);
    minChi = 999
    mintf = None
    f2 = None
    formular = ""

    if fittype == "pol2_c0":
        fitname = pol2_c0
        graph.Fit(fitname, "QF0", "", fit_min, fit_max)
        f2 = graph.GetFunction(fittype).Clone()
        if f2.GetChisquare() / f2.GetNDF() < minChi:
            minChi = f2.GetChisquare() / f2.GetNDF() 
            mintf = f2
            formular = "y = %.2fx + %.3fx^{2}" % (f2.GetParameter(0), f2.GetParameter(1))
    
    if mintf is not None:
        mintf.SetLineColor(4+graph.GetLineColor())
        mintf.SetLineWidth(2)
        mintf.Draw("same")
        return formular, f2.Eval(56), minChi

'''
    fitnames = ["pol1", "pol2", "expo"]
    if fittype != "":
        fitnames = [fittype]

    for fitname in fitnames:
        try:
            # graph.Fit(fitname, "QF0", "")
            graph.Fit(fitname, "QF0", "", 0, 60)
            f2 = graph.GetFunction(fitname).Clone()
            if f2.GetChisquare() / f2.GetNDF() < minChi:
                minChi = f2.GetChisquare() / f2.GetNDF() 
                mintf = f2
                if fitname == "pol1":
                    formular = "R(#mu) = %.2f + %.2f#mu" % (f2.GetParameter(0), f2.GetParameter(1))
                if fitname == "pol2":
                    formular = "R(#mu) = %.2f + %.2f#mu + %.2f#mu^{2}" % (f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2) )
                if fitname == "expo":
                    formular = "R(#mu) = e^{%.2f + %.2f#mu}" % (f2.GetParameter(0), f2.GetParameter(1))

        except:
            continue
'''

def DrawGraphs(pad, gMap, fittype=""):
    pad.cd()
    leg = ROOT.TLegend(0.1,0.7,0.7,0.9)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    leg.SetTextFont(62)

    mg = ROOT.TMultiGraph()
    color = 1
    #ymax = 0
    seedname = ""
    for k,v  in gMap.items():
        if color == 3:
            color = 4
        v.SetMarkerColor(color)
        v.SetLineColor(color)
        v.SetMarkerStyle(20)
        v.SetMarkerSize(1)
        seedname = v.GetTitle()
        color +=1
        #ymax = max(ymax, v.Eval(60))
        mg.Add(v.Clone())

    mg.Draw("AP")
    #mg.GetXaxis().SetRangeUser(0, plot_max)
    #mg.GetYaxis().SetRangeUser(0, ymax * 1.5)
    mg.GetYaxis().SetRangeUser(0, ymax)
    mg.GetXaxis().SetTitle("PileUp")
    mg.GetXaxis().SetLimits(fit_min, fit_max)
    #print("object mg1 is called")
    #c.Update()
    #mg.GetYaxis().SetTitle("Rate (nBunches = %d) [%s]" % (nBunches, unit) )
    mg.GetYaxis().SetTitle("Rate [kHz]")
    mg.GetYaxis().SetTitleOffset(1.3)
    # mg.GetYaxis().SetRangeUser(0, 50)
    ROOT.SetOwnership(mg, False)
    pad.Update()

    for k,v  in gMap.items():
        formular, rate, minChi =FittingGraphs(v, fittype)
        print ",%s,%f" %(k, rate),
        #newleg = "%s : %s [%s]" % (k, formular, v.GetTitle())
        #newleg = "%s: #bf{%s,  #chi^{2}/NDF=%.2f}" % (k, formular, minChi)
        newleg = "%s: #bf{%s}" % (k, formular)
        leg.AddEntry(v, newleg, 'p')
    leg.Draw()

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextAlign(13)
    # tex.SetTextFont(61)
    #tex.SetTextColor(ROOT.kBlack)
    #tex.SetLineWidth(2)
    tex.SetTextSize(0.03)
    tex.DrawLatex(0.11,0.935,"#scale[1.5]{CMS}")
    tex.SetTextSize(0.035)
    tex.DrawLatex(0.63, 0.93, "#bf{2018 Data (13TeV)}")
    #CMS_lumi.CMS_lumi(pad, 13, 11, 0)


def DrawRatios(dpad, glist, denlabel="MC"):
    dpad.cd()
    dpad.Draw()
    # denlabel = "MC"
    dengraph = glist[denlabel]
    denx =  [dengraph.GetX()[i] for i in range(dengraph.GetN())]
    ratioMap = {}

    for k,v in glist.items():
        if k == denlabel:
            continue
        ratioMap[k] = ROOT.TGraphErrors(v.GetN())
        for i in range(1, v.GetN()):
            # print v.GetN(), v.GetX()[i], len(denx)
            if v.GetX()[i] <= len(denx):
                j = denx.index(v.GetX()[i])
            else: 
                j = 0
            # j = denx.index(v.GetX()[i])
            # j = np.where(denx == v.GetX()[i])
            if dengraph.GetY()[j] > 0:
                ratioMap[k].SetPoint(i, v.GetX()[i], v.GetY()[i]/dengraph.GetY()[j])
            else:
                ratioMap[k].SetPoint(i, v.GetX()[i], 0.0)
            ### Error
            sigA = v.GetEY()[i]/v.GetY()[i] if v.GetY()[i]> 0 else 0
            sigB = dengraph.GetEY()[j]/dengraph.GetY()[j] if dengraph.GetY()[j]> 0 else 0
            ratioMap[k].SetPointError(i, 0, ratioMap[k].GetY()[i] * math.sqrt(sigA**2 + sigB**2))

    mg2 = ROOT.TMultiGraph()
    leg = ROOT.TLegend(0.1556476,0.6198268,0.7018421,0.8051783)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(62)
    ymax = 0
    for k,v  in ratioMap.items():
        v.SetMarkerColor(glist[k].GetMarkerColor())
        v.SetLineColor(glist[k].GetLineColor())
        v.SetMarkerStyle(21)
        v.SetMarkerSize(1)
        #ymax = max(ymax, v.Eval(60))
        mg2.Add(v)
        leg.AddEntry(v, k, 'p')

    ROOT.SetOwnership(mg2, False)
    mg2.Draw("AP")
    mg2.GetXaxis().SetTitle("")
    mg2.GetXaxis().SetLabelSize(0.15)
    print("object mg2 is called")
    #mg2.GetXaxis().SetLimits(10, 60)
    mg2.GetXaxis().SetRangeUser(0, 60)
    mg2.GetYaxis().SetLabelSize(0.15)
    mg2.GetYaxis().SetTitleSize(0.1)
    mg2.GetYaxis().SetTitleOffset(0.3)
    mg2.GetYaxis().CenterTitle(True)
    mg2.GetYaxis().SetTitle("#frac{Data}{%s}" % denlabel)
    #mg2.GetYaxis().SetRangeUser(0, ymax * 1.5)
    mg2.GetYaxis().SetRangeUser(0, 250)
    tline = ROOT.TLine()
    tline.SetLineColor(4)
    tline.DrawLine(10, 1, 60, 1)
    leg.Draw()

freq = 11245.6
unit = "kHz"
nBunches = 2544
pubins = np.arange(0, 120, 1)
pumap = collections.defaultdict(list)

filedir = {

    "1.5E34":"/eos/uscms/store/user/huiwang/L1Menu2017/Dec11fill_7118_nanoDST_shifter_Prescale_2018_v2_1_0_Col_1.5/*_Default_PU.csv",
    "1.7E34":"/eos/uscms/store/user/huiwang/L1Menu2017/Dec11fill_7118_nanoDST_shifter_Prescale_2018_v2_1_0_Col_1.7/*_Default_PU.csv",
    "2.0E34":"/eos/uscms/store/user/huiwang/L1Menu2017/Dec11fill_7118_nanoDST_shifter_Prescale_2018_v2_1_0_Col_2.0/*_Default_PU.csv",
    #"simulated 8b":"/eos/uscms/store/user/huiwang/L1Menu2017/Mar13fill_7358_nanoDST_Prescale_2018_v2_1_0_Col_2.0_8b/*_Default_PU.csv",
    #"origin 12b":"/eos/uscms/store/user/huiwang/L1Menu2017/Mar13fill_7358_nanoDST_Prescale_2018_v2_1_0_Col_2.0/*_Default_PU.csv",
    #"bx 5 to 10":"/eos/uscms/store/user/huiwang/L1Menu2017/Mar13fill_7358_nanoDST_Prescale_2018_v2_1_0_Col_2.0_5to10/*_Default_PU.csv",
}


seedMap = {
    # "L1_SingleMu20"      : "pol1",
    # "L1_SingleMu20er"    : "pol1",
    # "L1_DoubleMu_12_5"   : "pol2",
    # "L1_DoubleMu18er2p1"   : "pol2",
    "L1APhysics"   : "pol2_c0",
    # "L1_SingleMu5"       : "pol1",
    # "L1_DoubleMu0"       : "pol2",
    # "L1_SingleEG32"      : "pol1",
    # "L1_SingleIsoEG30"   : "pol1",
    # "L1_SingleIsoEG28er" : "pol1",
    # "L1_DoubleEG_25_12"  : "pol2",
    # "L1_SingleJet180"    : "pol1",
    # "L1_SingleJet60"     : "pol1",
    # "L1_ETM100"          : "expo",
    # "L1_ETM105"          : "expo",
    # "L1_HTT220"          : "expo",
}





if __name__ == "__main__":
    c = ROOT.TCanvas("fd","Fdf", 600, 600 )
    ROOT.gStyle.SetOptStat(000000000)
    c.SetLeftMargin(0.11);
    #tdrstyle.setTDRStyle()
    doRatio = True
    doRatio = False



    # df = pd.read_csv("fROOT/Jul03v6PUS_v9/2017_v96p8_0_EmuPU_PU.csv", index_col=None, header=0)
    for seed, fit in seedMap.items():
    # for seed in pd.unique(df.L1Seed):
        # if seed == "L1_DoubleMu18er2p1" or seed == "L1_DoubleMu22er2p1"  :
            # continue
    # for seed in pd.unique(df[df.L1Seed.str.contains("ETMHF70")]['L1Seed']):
    # for seed in pd.unique(df[df.L1Seed.str.contains("ETM80")]['L1Seed']):
        fit = seedMap.get(seed, "")
        if fit =="" and "Single" in seed:
            fit = "pol1"
        if fit =="" and "Double" in seed or "Triple" in seed:
            fit = "pol2"
        if fit =="" and "ETM" in seed :
            fit = "expo"
        c.Clear()
        c.Update()
        upad = c
        if doRatio:
            upad = ROOT.TPad("pPlot","",0.05,0.26,0.99,0.99)
            upad.Draw()
            dpad = ROOT.TPad("pRatio","",0.05,0.01,0.99,0.25)
            dpad.SetTopMargin(0.0)
            dpad.Draw()
            dpad.SetGrid()
        glist = {}
        for k,v in filedir.items():
            if k == "Emulator":
                seed = "L1_ETM100"
            glist[k] = GetPUGraph(glob.glob(v), seed)
        #DrawGraphs(upad, glist, fit)
        print seed, DrawGraphs(upad, glist, fit)
        if doRatio:
            # DrawRatios(dpad, glist, denlabel="v2 PUS")
            DrawRatios(dpad, glist, denlabel="v96.8")
            # DrawRatios(dpad, glist, denlabel="2016 Unp Obj")
	c.SetGrid()
        c.Update()
        outdir = "plots"
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        c.SaveAs("%s/%s.root" % (outdir, seed))
        c.SaveAs("%s/%s.png" % (outdir, seed))
        c.SaveAs("%s/%s.pdf" % (outdir, seed))
        c.SaveAs("%s/%s.C" % (outdir, seed))

