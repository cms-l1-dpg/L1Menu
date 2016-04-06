from ROOT import *
import ROOT
from sys import argv
import os
import copy
from Config import DualMap
import re
import matplotlib.pyplot as plt

NORM=False
LOG=True
# obj = "SingleJet"
# obj = "HTT"
# obj = "ETM22"
# obj = "ETM50"
# obj = "HTM"
obj = "SingleMu"
dump =[]

gStyle.SetOptStat(False)
fileNames=["../results/r259721_tsgv4plusLayer1_ETM.root"]
# gROOT.SetBatch(True)

class EFFFIT:
   def __call__( self, x, par ):
       value=par[0]/2.+par[0]/2.*ROOT.TMath.Erf((x[0]-par[1])/par[2]);
       return value

def getall(d, basepath=""):
    "Generator function to recurse into a ROOT file/dir and yield (path, obj) pairs"
    for key in d.GetListOfKeys():
        kname = key.GetName()
        if key.IsFolder():
            # TODO: -> "yield from" in Py3
            for i in getall(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)


def LineInterp(lowbin, highbin, y):
    x0, y0 = lowbin
    x1, y1 = highbin
    return x0 + (y - y0) * (x1 - x0) / (y1-y0)


def Get95Eff(g):

    ########## Beginning of LA  ############################    
    g.Draw();

    gx = [g.GetX()[i] for i in range(g.GetN())]
    gy = [g.GetY()[i] for i in range(g.GetN())]
    # xminf=0.
    # xminf=100.
    for x,y in zip(gx, gy):
        if y >= 0.1:
            xminf = x
            break
    # xminf=40
    xmaxf= gx[-1]
    width=20.
    mean=200.

    fitter = TF1("fitf",EFFFIT(),xminf,xmaxf,3);
    fitter.SetParameters(1,mean,width);
    fitter.FixParameter(0,1.0);

    fitter.SetLineColor(ROOT.kRed)
    fitter.SetLineStyle(1)
    g.Fit(fitter,"0RQ");
    dump.append(fitter)

    # print "Turnon is at 95% for pT = ",fitter.GetX(0.95)

    fitter.Draw("sames")

    c.Update();
    return fitter.GetX(0.95), fitter
    # raw_input('\npress return to continue...')

    ########## End of LA  ############################
    lowbin = None
    highbin = None
    per = 0.8
    gx = [g.GetX()[i] for i in range(g.GetN())]
    gy = [g.GetY()[i] for i in range(g.GetN())]
    # gy = list(g.GetY())
    for x, y in zip(gx, gy):
        if y < per:
            lowbin = (x, y)
        if y >= per and highbin is None:
            highbin = (x, y)
    return LineInterp(lowbin, highbin, per)

if __name__ == "__main__":
    try: fileNames=argv[1:]
    except:
        print "No files specified"
        exit()

    if not os.path.exists("plots"):
        os.mkdir("plots")

    files=[]
    for fileName in fileNames:
        files.append(TFile(fileName))


    c=TCanvas()

    l = 0
    x =[]
    y = []
    mg = TMultiGraph()

    for k, o in getall(files[0]):
        if "Eff" in k and obj in k and type(o)==type(TGraphAsymmErrors()):
            l = l +1
            o.SetLineColor(l)
            o.GetXaxis().SetTitle("offline %s [GeV]" % obj)
            o.GetYaxis().SetTitle("Efficiency")
            # mg.Add(o)
            # continue
            m = re.match("Eff/L1_%s(\d+)_Pt" % obj, k)
            f95, fit = Get95Eff(o)
            print l, k, f95
            fit.SetLineColor(l)
            if m is not None:
                x.append(m.group(1))
                y.append(f95)
            # fit.Draw()

    c.Update()
    # mg.Draw("AP")
    # mg.GetXaxis().SetTitle("offline %s [GeV]" % obj)
    # mg.GetYaxis().SetTitle("Efficiency")
    raw_input('\npress return to continue...')
    c.SaveAs("turnon.png")
    c.SaveAs("turnon.root")

    plt.scatter(x, y)
    plt.xlabel("Online %s" % obj)
    plt.ylabel("offline %s with 95%% eff" % "ETM")
    plt.show()
    plt.savefig("dia.png")
