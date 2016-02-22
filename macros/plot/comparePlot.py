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

gStyle.SetOptStat(False)
gROOT.SetBatch(True)

try: fileNames=argv[1:]
except:
    print "No files specified"
    exit()

if not os.path.exists("plots"): os.mkdir("plots")

files=[]
for fileName in fileNames:
    files.append(TFile(fileName))

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

c=TCanvas()

# allhist = getall(files[0])
mg = ROOT.TMultiGraph()

l = 0
x =[]
y = []
obj = "HTT"
# obj = "HTM"

for k, o in getall(files[0]):
    if "Eff" in k and obj in k and type(o)==type(TGraphAsymmErrors()):
    # if "Eff" in k and "HTT" in k and type(o)==type(TGraphAsymmErrors()):
        o.SetLineColor(l)
        m = re.match("Eff/L1_%s(\d+)_Pt" % obj, k)
        # m = re.match("Eff/L1_HTT(\d+)_Pt", k)
        print k ,  Get95Eff(o)
        if m is not None:
            x.append(m.group(1))
            y.append( Get95Eff(o))
        # o.SetTitle
        mg.Add(o)
        l = l +1

mg.Draw("APL")
mg.GetXaxis().SetTitle("offline %s [GeV]" % obj)
# mg.GetXaxis().SetTitle("offline HT [GeV]")
mg.GetYaxis().SetTitle("Efficiency")
# mg.Draw("a fb l3d")
c.SaveAs("ff.png")
plt.scatter(x, y, )
plt.xlabel("Online %s" % obj)
plt.ylabel("offline HT with 80\% eff")
plt.show()
