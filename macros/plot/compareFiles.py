from ROOT import *
from sys import argv
import os

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

keys=getall(files[0])

c=TCanvas()
for k, o in getall(files[0]):
    if type(files[0].Get(k))!=type(TH1F()) and \
       type(files[0].Get(k))!=type(TProfile()) and \
       type(files[0].Get(k))!=type(TH1D()) and \
       type(files[0].Get(k))!=type(TGraphAsymmErrors()):
        continue

    l=TLegend(.49,.69,.89,.89)
    l.SetFillStyle(0)
    hists=[]
    max=0

    for lp in range(len(files)):
        file=files[lp]
        fileName=fileNames[lp]

        file.cd()
        h=file.Get(k)
        #h.Rebin(4)
        hists.append(h)

        if NORM:h.Scale(1./h.Integral())

        if h.GetMaximum()>max: max=h.GetMaximum()

        h.SetLineColor(1+lp)
        h.SetLineWidth(3)
        legname = fileNames[lp].split('/')[-1]
        legname = legname.split('.')[0]
        legname = legname.replace("results_", "")
        legname = legname.replace("_Menu", "")
        legname = legname.replace("_RATE", "")
        l.AddEntry(h,legname.replace('_', ' '),"l")

    for lp in range(len(hists)):
        h=hists[lp]

        if lp==0:
            if LOG or k == "jet_eta":
              c.SetLogy(True)
              h.SetMaximum(100*max)
            else:
              c.SetLogy(False)
              # h.SetMaximum(2*max)
              h.SetMinimum(0)
            if type(h) == type(TGraphAsymmErrors()):
                c.SetLogy(False)
                h.Draw("AP")
            else:
                h.Draw("hist")
        else:
            if type(h) == type(TGraphAsymmErrors()):
                h.Draw("P")
            else:
                h.Draw("histSAME")
    l.Draw("SAME")
    kname = k.replace("/", "_")
    c.SaveAs("plots/"+kname+".pdf")
    c.SaveAs("plots/"+kname+".root")
    c.SaveAs("plots/"+kname+".png")

