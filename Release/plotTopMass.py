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

def set3HistsColor( hist, hist2, hist3 ) :
  hist.SetMarkerColor(ROOT.kRed)
  hist.SetLineColor(ROOT.kRed)
  hist.SetMarkerStyle(21)
  hist2.SetMarkerColor(ROOT.kMagenta)
  hist2.SetLineColor(ROOT.kMagenta)
  hist2.SetMarkerStyle(22)
  hist3.SetMarkerColor(ROOT.kOrange)
  hist3.SetLineColor(ROOT.kOrange)
  hist3.SetMarkerStyle(23)
def set6HistsColor( hist, hist2, hist3,hist4,hist5,hist6 ) :
  setHistsColor(hist, hist2, hist3)
  
  hist4.SetMarkerColor(ROOT.kBlue)
  hist4.SetLineColor(ROOT.kBllue)
  hist5.SetMarkerColor(ROOT.kBlack)
  hist5.SetLineColor(ROOT.kBlack)
  hist6.SetMarkerColor(ROOT.kCyan)
  hist6.SetLineColor(ROOT.kCyan)



def anaTree( tree, tree2,nEvent ) :
  label="reco"

  pq_up = "(quality>0.0&& bjet_nCharged[0]>2 && bjet_nCharged[1]>2)"

  matched_jet1 = "&&(abs(bjet_partonPdgId[0])==5)"
  matched_jet2 = "&&(abs(bjet_partonPdgId[1])==5)"

  matched_or = "&&(abs(bjet_partonPdgId[0])==5 || abs(bjet_partonPdgId[1])==5)"
  matched_and = "&&(bjet_partonPdgId[0]*bjet_partonPdgId[1]==-25)"

  matched_first = "&&(bjet_partonPdgId[0]*lep_charge[0])"
  

  largeCharge = "&&((lep_charge[0]*bjet_charge[0]+lep_charge[1]*bjet_charge[1]) <0)"
  
  diffCharge = "&&( abs(bjet_charge[0]-bjet_charge[1])>=10)" 
  sameCharge = "&&( abs(bjet_charge[0]-bjet_charge[1])<10)" 

  btagged = "&&(bjet_btag[0]&&bjet_btag[1])"
  btagged_jet1 ="&&(bjet_btag[0])" 
  btagged_jet2 ="&&(bjet_btag[1])" 
  btagged_or   ="&&(bjet_btag[0] || bjet_btag[1])"

  lepJet1Pair = "&&(lep_charge[0]*bjet_partonPdgId[0]==5)"
  lepJet2Pair = "&&(lep_charge[1]*bjet_partonPdgId[1]==5)"
  lepJetPairs = lepJet1Pair+lepJet2Pair


  qualityBin = [nxbin, xmin, xmax]
  chargeBin = [nqbin, qmin, qmax]
  lqBin = [nlqbin, lqmin, lqmax]

 
  
  outFile= TFile("histOut.root","RECREATE")
  h1   = getTH1("top_mass ; M_{top}[GeV/c^2] ; Entries", [100,100,300],tree, "top_mass" ,pq_up)
  h2   = getTH1("top_mass_btagged ; M_{top}[GeV/c^2] ; Entries", [100,100,300],tree, "top_mass" ,pq_up+btagged_or)
  h3   = getTH1("top_mass_charged ; M_{top}[GeV/c^2] ; Entries", [100,100,300],tree2, "top_mass" ,pq_up)
  h4   = getTH1("top_mass_btagged_charged ; M_{top}[GeV/c^2] ; Entries", [100,100,300],tree2, "top_mass" ,pq_up+btagged_or)
  h1.SetName(h1.GetTitle())
  h2.SetName(h2.GetTitle())
  h3.SetName(h3.GetTitle())
  h4.SetName(h4.GetTitle())

  h1.Sumw2()
  h2.Sumw2()
  h3.Sumw2()
  h4.Sumw2()

  h1.Scale(1./nEvent)
  h2.Scale(1./nEvent)
  h3.Scale(1./nEvent)
  h4.Scale(1./nEvent)

  h1.Write()
  h2.Write()
  h3.Write()
  h4.Write()

  outFile.Write()
  outFile.Close()

  fitmodel = "gaus"
  c1 = makeCanvas("topMass")
  h1.Fit(fitmodel,"S")
  h1.SetLineColor(ROOT.kRed)
  h3.Fit(fitmodel,"S")
  h3.SetLineColor(ROOT.kBlue)
  h1.Draw()
  h3.Draw("same")

  c1.SaveAs("topMass.png")
  c1.SaveAs("plotCode/topMass.C")

  c2 = makeCanvas("topMass_btaggedOR")
  h2.Fit("gaus","S")
  h2.SetLineColor(ROOT.kRed)
  h4.Fit(fitmodel,"S")
  h4.SetLineColor(ROOT.kBlue)
  h2.Draw()
  h4.Draw()
  c2.SaveAs("topMass_btaggedOR.png")
  c2.SaveAs("plotCode/topMass_btaggedOR.C")


if __name__ == "__main__" :
  if len(sys.argv) != 2 :
    print "Wrong argument!"
    sys.exit(-1)
  else :
    file0 = TFile.Open(sys.argv[1])
    tree = file0.Get("JetTree")
    tree2 = file0.Get("JetTreeCharge")
    nEvent = file0.Get("nEvent").GetBinContent(1)
    anaTree(tree,tree2, nEvent)
