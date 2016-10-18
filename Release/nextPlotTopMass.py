#!/usr/bin/env python

from ROOT import *
import json, os, sys, math, getopt 
from CATTools.CatAnalyzer.histoHelper import *
gROOT.SetBatch(True)

gStyle.SetPalette(1)

xmin=0.0
xmax=2.2
nxbin=20

qmin=-100.5
qmax= 100.5
nqbin = 21

lqmin=-5
lqmax=6
nlqbin=11


xsec = [87.31, 18610, 6025.2]
lumi = 30 * 1000

def plotHist() :
  files = []
  files.append(TFile.Open("histOut_TTLL.root"))
  files.append(TFile.Open("histOut_DY10to50.root"))
  files.append(TFile.Open("histOut_DY50.root"))

  for idx, file in enumerate(files) :
    c1 = makeCanvas("top1Mass")
    h_top1_mass =[]
    h_top2_mass =[]
    h_top1_mass.append(file.Get("top1_mass_GenJet"))
    h_top1_mass.append(file.Get("top1_mass_Jet"))
    h_top1_mass.append(file.Get("top1_mass_charged_Jet"))

    h_top1_mass[0].SetLineColor(ROOT.kBlack)
    h_top1_mass[1].SetLineColor(ROOT.kRed)
    h_top1_mass[2].SetLineColor(ROOT.kBlue)

    h_top1_mass[0].SetMarkerColor(ROOT.kBlack)
    h_top1_mass[1].SetMarkerColor(ROOT.kRed)
    h_top1_mass[2].SetMarkerColor(ROOT.kBlue)

    h_top1_mass[0].Scale(xsec[idx]*lumi)
    h_top1_mass[1].Scale(xsec[idx]*lumi)
    h_top1_mass[2].Scale(xsec[idx]*lumi)

    h_top1_mass[1].Draw()
    h_top1_mass[2].Draw("same")
    h_top1_mass[0].Draw("same")
    c1.SaveAs("top1Mass_%d.png"%(idx))



if __name__ == "__main__" :
  if len(sys.argv) != 1 :
    print "Wrong argument!"
    sys.exit(-1)
  else :
      #plotOutFileName = "histOut_"+os.path.splitext(input)[0]+".root"  # [ histOut_TTLL.root, histOut_DY50.root, histOut_DY10to50.root ]
      plotHist()


