#!/usr/bin/env python

# File        : GetRateStep.py
# Author      : Zhenbin Wu
# Contact     : zhenbin.wu@gmail.com
# Date        : 2016 Feb 29
#
# Description : 


from ROOT import *
import ROOT
from sys import argv
import os
import copy
from Config import DualMap
import re
import matplotlib.pyplot as plt

filename = "./r259721_tsgv3_rate.root"
folder = "Rate"
margin2D = 0.0001
fraction = [1, 0.55, 0.25, 0, -0.2, -0.35, -0.5, -0.6]
objectStart = {
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SingleEG ~~~~~
    "nEGVsPt"         : 40, # SingleEG
    # "nEGErVsPt"       : 100, # SingleEGer
    "nIsoEGVsPt"      : 27, # SingleIsoEG
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SingleMu ~~~~~
    "nMuVsPt"         : 20, # SingleMu
    "nMuErVsPt"       : 16, # SingleMu         |#eta|<2.1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SingleJet ~~~~~
    # "nJetCenVsPt"     : 100, # SingleJetCentral
    "nJetVsPt"        : 150, # SingleJet
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ SingleTau ~~~~~
    # "nIsoTauErVsPt"   : 100, # SingleIsoTauEr
    # "nIsoTauVsPt"     : 100, # SingleIsoTau
    # "nTauErVsPt"      : 100, # SingleTauer
    # "nTauVsPt"        : 100, # SingleTau
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sums ~~~~~
    # "nETMVsETM"       : 100, # ETM
    # "nETTVsETT"       : 100, # ETT
    # "nHTMVsHTM"       : 100, # HTM
    "nHTTVsHTT"       : 280, # HTT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DoubleTrigger on the low pt leg ~~~~~
    # "nDiCenJetVsPt"   : 100, # DiCenJet
    # "nDiEGVsPt"       : 100, # DiEG
    # "nDiIsoEGVsPt"    : 100, # DiIsoEG
    # "nDiJetVsPt"      : 100, # DiJet
    # "nDiTauVsPt"      : 100, # DiTau
    # "nQuadCenJetVsPt" : 100, # QuadCenJet
    # "nQuadJetVsPt"    : 100, # QuadJet
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ DoubleTrigger in 2D ~~~~~
    "nAsymDiCenJetVsPt" : (100, 100), # DiCenJet
    # "nAsymDiJetVsPt"    : [100, 100], # DiJet
    # "nEGPtVsPt"         : [100, 100], # DoubleEle
    # "nIsoEGPtVsPt"      : [100, 100], # DoubleIsolEle
    # "nMuPtVsPt"         : [100, 100], # DoubleMu
    # "nOniaMuPtVsPt"     : [100, 100], # DoubleMu_Er_HighQ_WdEta22 (Quarkonia)
}




def GetRateVariation2D(h, value):
    if len(value) != 2:
        print "Wrong value for 2D!"
        return False
    orgrate = h.GetBinContent(h.FindBin(value[0], value[1]))

    for frac in fraction:
        ibins = GetRate2D(h, orgrate * (1+frac) )
        # at %d, with rate %d+-%d Hz" % \
        if len(ibins) == 0:
            print "Varying %d%%(+-%d%%), thresholds and rates : None" % (frac *100, margin2D*100)
            continue
        print "Varying %d%%(+-%d%%), thresholds and rates :" % (frac *100, margin2D*100),
        for ibin in ibins:
            print "[%d, %d] with %d+-%d Hz" % \
            (h.GetXaxis().GetBinCenter(ibin[0]),h.GetYaxis().GetBinCenter(ibin[1]), \
             h.GetBinContent(ibin[0], ibin[1]), h.GetBinError(ibin[0], ibin[1]))





def GetRateVariation1D(h, value):
    bVali = {}
    for i in range(1, h.GetNbinsX()):
        bVali[h.GetBinCenter(i)] = i
        # print h.GetBinCenter(i), h.GetBinContent(i), h.GetBinError(i)
    orgrate =  h.GetBinContent(bVali[value])
    for frac in fraction:
        ibin = GetRate1D(h, orgrate * (1+frac) )
        print "Varying %d%%, threshold at %d, with rate %d+-%d Hz" % \
        (frac *100, h.GetBinCenter(ibin), h.GetBinContent(ibin), h.GetBinError(ibin))


def GetRate2D(h, rate):
    matchbin = []
    for i in range(1, h.GetNbinsX()):
        for j in range(1, h.GetNbinsY()):
            bincont = h.GetBinContent(i, j)
            if fabs(float(bincont) / rate -1 ) <= margin2D:
                matchbin.append([i, j])
    return matchbin

def GetRate1D(h, rate):
    valueMap = {}
    for i in range(1, h.GetNbinsX()):
        valueMap[h.GetBinContent(i)] = i
    minValue =  min(valueMap.keys(), key=lambda x:abs(x-rate))
    return valueMap[minValue]

def GetRateVariation(h, value):
    if isinstance(value, tuple):
        return GetRateVariation2D(h, value)
    else:
        pass
        return GetRateVariation1D(h, value)

if __name__ == "__main__":
    file = ROOT.TFile(filename, "OPEN")
    for k,v in objectStart.items():
        f = file.Get("Rate/%s" % k)
        print "------- for %s" % f.GetTitle()
        GetRateVariation(f, v)
