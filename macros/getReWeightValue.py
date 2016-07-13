from ROOT import *

dataFile=TFile.Open("myDataDistributions_274199.root")
mcFile=TFile.Open("myMCDistributions.root")

dataHist=dataFile.Get("nVtex")
mcHist=mcFile.Get("nVtex")

dataHist.SetBinContent(1,0)
dataHist.Scale(1./dataHist.Integral())
mcHist.Scale(1./mcHist.Integral())

weights=[]

for b in range(1,52):
    if mcHist.GetBinContent(b)!=0:
        weight=dataHist.GetBinContent(b)/mcHist.GetBinContent(b)
    else:
        weight=0
    weights.append(weight)

#weights.append(0)

print weights
print len(weights)
