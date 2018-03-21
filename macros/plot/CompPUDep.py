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

#foldername = "Mar12fill_6356_6360_Prescale_2018_v0_1_0_Col_1.0_408"
foldername = "Mar10fill_6356_6360_Prescale_2018_v0_1_0_Col_1.6"

fit_min = 20
fit_max = 75
plot_min = 0
plot_max = 90
maxy = 500
PU = 61
freq = 11245.6
config = 2017

fitname = "pol4"
#fitname = "expo"
#fitname = ROOT.TF1("fitname","[0]*x + [1]*x*x",0,80);
#fitname.SetParameters(0.1,0.001);

#filedir = "/eos/uscms/store/user/huiwang/L1Menu2017/Sep12v2.2_menu_2.2_v96p20_v8_run_301912_to_302029/*Default_PU.csv"
#filedir = "/eos/uscms/store/user/huiwang/L1Menu2017/Aug17v2.2_menu_2.2_v6_v96p20_IgnorePrescale/*Default_PU.csv"
#filedir = "/eos/uscms/store/user/huiwang/L1Menu2017/Sep19Prescale_Sets_RUN_302485_col_2.0/*Default_PU.csv"
#filedir = "/eos/uscms/store/user/huiwang/L1Menu2017/Jan21fill_6358_col_1.6_core/*Default_PU.csv"
filedir = "/eos/uscms/store/user/huiwang/L1Menu2017/" + foldername + "/*Default_PU.csv"

if config == 2016:
    nBunches = 2208
    unit = "kHz"
if config == 2017:
    #nBunches = 2544
    nBunches = 1866
    unit = "kHz"

if config == 1:
    nBunches = 1
    unit = "Hz"

pubins = np.arange(plot_min,plot_max, 1)
pumap = collections.defaultdict(list)

PatMap = {  
    "L1APhysics" : "L1APhysics",
    #"L1_SingleMu25" :  "L1_SingleMu25",
    #"L1_DoubleMu_15_7" : "L1_DoubleMu_15_7",
    #"L1_SingleEG40" : "L1_SingleEG40",
    #"L1_SingleIsoEG32" : "L1_SingleIsoEG32",
    #"L1_SingleIsoEG30er2p1" : "L1_SingleIsoEG30er2p1",
    #"L1_DoubleEG_25_14" : "L1_DoubleEG_25_14",
    #"L1_DoubleIsoTau34er2p1" : "L1_DoubleIsoTau34er2p1",
    #"L1_SingleJet180" : "L1_SingleJet180",
    # 'L1_DoubleJet150er2p7': 'L1_DoubleJet150er2p7',
    # 'L1_DoubleMu_15_5_SQ': 'L1_DoubleMu_15_5_SQ',
    # 'L1_HTT320er': 'L1_HTT320er',
    # 'L1_HTT340er': 'L1_HTT340er',
    # 'L1_HTT380er': 'L1_HTT380er',
    # 'L1_HTT400er': 'L1_HTT400er',
    # 'L1_ETMHF120': 'L1_ETMHF120',
    # 'L1_TripleJet_105_85_76_VBF': 'L1_TripleJet_105_85_76_VBF',
    # 'L1_TripleEG_14_10_8': 'L1_TripleEG_14_10_8',
    # 'L1_SingleIsoEG34er2p1': 'L1_SingleIsoEG34er2p1',
    # 'L1_SingleMu22_BMTF': 'L1_SingleMu22_BMTF',
    # 'L1_SingleMu22_EMTF': 'L1_SingleMu22_EMTF',
    # 'L1_SingleMu22_OMTF': 'L1_SingleMu22_OMTF',
    # 'L1_SingleEG40er2p5': 'L1_SingleEG40er2p5',
    # 'L1_SingleIsoEG32er2p5': 'L1_SingleIsoEG32er2p5',
    # 'L1_EG40er2p1_Tau20er2p1_dR_Min0p3': 'L1_EG40er2p1_Tau20er2p1_dR_Min0p3',
    # 'L1_TripleMu_0_0_0OQ_DoubleMu_0_0_OS_Mass_5to17': 'L1_TripleMu_0_0_0OQ_DoubleMu_0_0_OS_Mass_5to17',	#Hu Zhen seed 407
    # 'L1_TripleMu_0SQ_0SQ_0OQ_DoubleMu_0SQ_0SQ_OS_Mass_5to17': 'L1_TripleMu_0SQ_0SQ_0OQ_DoubleMu_0SQ_0SQ_OS_Mass_5to17',	#Hu Zhen seed 408
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
    #minx = min(x)
    #maxx = max(x)
    #maxyy = max(y)

    for i, (xx, yy, yee) in enumerate(zip(x, y, yerr)):
        # if yy != 0 and yee/yy >0.3:
            # continue
	#if i == 22 or i == 23 or i == 24:
	    # continue
        graph.SetPoint(i, xx, yy)
	#print (i,xx,yy,yee)
        #print "h1->SetBinContent( %d, %f);" %(xx,yy)
        #print "h1->SetBinError( %d, %f);" %(xx,yee)
        graph.SetPointError(i, 0, yee)

    #fit_p0_init = maxyy / fit_max
    #print "maxyy = ", maxyy, "fit_p0_init = ", fit_p0_init
    #fitname.SetParameters(fit_p0_init,0.001);

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
        graph.GetYaxis().SetRangeUser(0, maxy)
        graph.GetYaxis().SetTitle("Rate (nBunches = %d) [%s]" % (nBunches, unit))
    else:
        graph.Draw("P")
    canvas.Update()
    leg.AddEntry(graph, l1seed, "p")

    result_ptr = graph.Fit(fitname, "SQ", "", fit_min, fit_max)
    error_vec = result_ptr.GetConfidenceIntervals()
    #print ("error vec size = %d" % error_vec.size())
    #f2 = graph.GetFunction("fitname").Clone()
    f2 = graph.GetFunction(fitname).Clone()
    f2.SetLineColor(1+count)
    f2.SetLineWidth(2)
    #for i in range (fit_min, fit_max):
      #print("bin = %d, Obeserve = %.2f , Expect = %.2f\n" % (i, graph.Eval(i), f2.Eval(i)) )
      #print("bin = %d, Obeserve = %.2f , Expect = %.2f , error %.2f \n" % (i, graph.Eval(i), f2.Eval(i), error_vec.at(i-fit_min)) )
    if error_vec.size() > 60 - fit_min:
      #print f2.Eval(50), "+-", error_vec.at(50-fit_min), f2.Eval(60), "+-", error_vec.at(60-fit_min)
      print (",%.2f,+-,%.2f,%.2f,+-,%.2f") % (f2.Eval(50), error_vec.at(50-fit_min), f2.Eval(60), error_vec.at(60-fit_min))
    if error_vec.size() <= 60 - fit_min:
      #print f2.Eval(50), "+-", 0.00, f2.Eval(60), "+-", 0.00
      print (",%.2f,+-,%.2f,%.2f,+-,%.2f") % (f2.Eval(50), 0.00, f2.Eval(60), 0.00)
    minChi = f2.GetChisquare() / f2.GetNDF()
    #fun = "Fit = %.2f + %.2f*x + %.3f*x^2" % (f2.GetParameter(0), f2.GetParameter(1), f2.GetParameter(2) )
    fun = "Fit = %f*x + %f*x^2" % (f2.GetParameter(0), f2.GetParameter(1) )
    # f2.Draw("same")
    #print("chi2 = ", f2.GetChisquare())
    #print("NDF = ", f2.GetNDF())
    #print fun
    f2_2 = f2.Clone("dashline2")
    f2_2.SetRange(fit_max, plot_max)
    f2_2.SetLineStyle(5)
    f2_2.Draw("same")
    key = "Rate(PU=%d): %.2f +- %.2f, chi2/NDF=%.2f" %(PU, f2_2.Eval(PU), error_vec.at(60-fit_min), minChi)
    #tex = ROOT.TLatex(0.15, 0.75, fun)
    #tex.SetNDC()
    #tex.SetTextAlign(13)
    #tex.SetTextFont(61)
    #tex.SetTextSize(0.04)
    #tex.SetTextColor(ROOT.kBlue)
    #tex.SetLineWidth(2)
    # tex.Draw()

    if key is not None:
        tex = ROOT.TLatex(0.3, 0.85, key)
        tex.SetNDC()
        tex.SetTextFont(61)
        tex.SetTextSize(0.055)
        tex.SetTextColor(ROOT.kGreen+2)
        tex.SetLineWidth(2)
        tex.Draw()

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
    #print key, " : ",  inputlist
    print key,

    for i, seed in enumerate(inputlist):
        DrawPU(c1, df, seed, i)
    leg.Draw()

    if config == 2017:
        l_PU = ROOT.TLine(PU, 0, PU, maxy)
        l_PU.SetLineColor(2)
        l_PU.SetLineWidth(2)
        l_PU.Draw()
        #l56 = ROOT.TLine(50, 0, 50, maxy)
        #l56.SetLineColor(2)
        #l56.SetLineWidth(2)
        #l56.Draw()
        #l60 = ROOT.TLine(60, 0, 60, maxy)
        #l60.SetLineColor(2)
        #l60.SetLineColor(2)
        #l60.SetLineWidth(2)
        #l60.Draw()

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
    #c1.SaveAs("plots/%s_%d_%s_PU%d.root" % (key, config, foldername, PU))
    c1.SaveAs("plots/%s_%d_%s_PU%d.png" % (key, config, foldername, PU))
    #c1.SaveAs("plots/%s_%d_%s_PU%d.pdf" % (key, config, foldername, PU))

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
    c1 = ROOT.TCanvas("fd","Fdf", 1200, 1000)
    leg = ROOT.TLegend(0.15,0.7,0.45,0.75)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(62)
    leg.SetTextSize(0.05)
    for k, v in PatMap.items():
        DrawL1(k, v)
        # wait()
