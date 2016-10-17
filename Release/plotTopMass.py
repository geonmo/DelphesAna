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



def anaTree( tree, tree2, nEvent, input, outFile, jetType ) :
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

  tags = { "": [tree,pq_up] ,"_btagged":[tree, pq_up+btagged],"_charged":[tree2, pq_up],"_btagged_charged":[tree2, pq_up+btagged] }


  qualityBin = [nxbin, xmin, xmax]
  chargeBin = [nqbin, qmin, qmax]
  lqBin = [nlqbin, lqmin, lqmax]

 
  
  massBin = [1000,100,200]
  largerTop = "top_mass[0]>top_mass[1]?top_mass[0]:top_mass[1]"
  smallTop = "top_mass[0]<top_mass[1]?top_mass[0]:top_mass[1]"

  hists =[]
  hists_norm=[]
  hists.append( getTH1("top1_mass_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin ,tree, largerTop ,pq_up))
  hists.append( getTH1("top1_mass_btagged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree, largerTop ,pq_up+btagged))
  hists.append( getTH1("top1_mass_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, largerTop ,pq_up))
  hists.append( getTH1("top1_mass_btagged_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, largerTop ,pq_up+btagged))

  hists.append( getTH1("top2_mass_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin ,tree, smallTop ,pq_up))
  hists.append( getTH1("top2_mass_btagged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree, smallTop ,pq_up+btagged))
  hists.append( getTH1("top2_mass_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, smallTop ,pq_up))
  hists.append( getTH1("top2_mass_btagged_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, smallTop ,pq_up+btagged))

  hists.append( getTH1("top_mass_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin ,tree, "top_mass" ,pq_up))
  hists.append( getTH1("top_mass_btagged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree, "top_mass" ,pq_up+btagged))
  hists.append( getTH1("top_mass_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, "top_mass" ,pq_up))
  hists.append( getTH1("top_mass_btagged_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, "top_mass" ,pq_up+btagged))

  hists.append( getTH1("ttbar_mass_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin ,tree, "ttbar_mass" ,pq_up))
  hists.append( getTH1("ttbar_mass_btagged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree, "ttbar_mass" ,pq_up+btagged))
  hists.append( getTH1("ttbar_mass_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, "ttbar_mass" ,pq_up))
  hists.append( getTH1("ttbar_mass_btagged_charged_%s ; M_{top}[GeV/c^2] ; Entries"%(jetType), massBin,tree2, "ttbar_mass" ,pq_up+btagged))

  for hist in hists :
   hist.SetName(hist.GetTitle())
   hist.Write()
   hist_norm = hist.Clone()
   hist_norm.Scale( 1./nEvent)
   hist_norm.SetName( hist.GetTitle().Strip()+"_Norm")
   hist_norm.Write() 


  histsForCal=[]
  for type in tags.keys() :
    temp1_hist = getTH1("correct_Pair_Lep_and_Jet1%s_%s; pdgId of Jets ; Entries"%(type, jetType),[11,-5.5,5.5] ,tags[type][0], "lep_charge[0]*bjet_partonPdgId[0]",tags[type][1])
    temp2_hist = getTH1("correct_Pair_Lep_and_Jet2%s_%s; pdgId of Jets ; Entries"%(type, jetType),[11,-5.5,5.5] ,tags[type][0], "lep_charge[1]*bjet_partonPdgId[1]",tags[type][1])
    sumHist = temp1_hist + temp2_hist
    sumHist.SetTitle("correct_Pair_Lep_and_Jet%s_%s"%(type, jetType)) 
    sumHist.SetName(sumHist.GetTitle())
    sumHist.Write()
 
  """
  fitmodel = "gaus"
  #fitmodel = "landau"
  c1 = makeCanvas("top1Mass_%s_%s"%(input,jetType))
  h1.Fit(fitmodel,"S")
  h1.SetLineColor(ROOT.kRed)
  h1.SetMarkerColor(ROOT.kRed)
  h3.Fit(fitmodel,"S")
  h3.SetLineColor(ROOT.kBlue)
  h3.SetMarkerColor(ROOT.kBlue)
  h1.Draw()
  h3.Draw("same")
  c1.BuildLegend(0.2393484,0.5348432,0.6190476,0.7456446)
  c1.SaveAs("topMass1_%s_%s.png"%(input,jetType))
  c1.SaveAs("plotCode/topMass1_%s_%s.C"%(input,jetType))

  c2 = makeCanvas("top1Mass_btaggedOR_%s_%s"%(input,jetType))
  h2.Fit(fitmodel,"S")
  h2.SetLineColor(ROOT.kRed)
  h2.SetMarkerColor(ROOT.kRed)
  h4.Fit(fitmodel,"S")
  h4.SetLineColor(ROOT.kBlue)
  h4.SetMarkerColor(ROOT.kBlue)
  h2.Draw()
  h4.Draw("same")
  c2.BuildLegend(0.2393484,0.5348432,0.6190476,0.7456446)
  c2.SaveAs("topMass1_btaggedOR_%s_%s.png"%(input,jetType))
  c2.SaveAs("plotCode/topMass1_btaggedOR_%s_%s.C"%(input,jetType))

  c3 = makeCanvas("top2Mass_%s_%s"%(input,jetType))
  h5.Fit(fitmodel,"S")
  h5.SetLineColor(ROOT.kRed)
  h5.SetMarkerColor(ROOT.kRed)
  h7.Fit(fitmodel,"S")
  h7.SetLineColor(ROOT.kBlue)
  h7.SetMarkerColor(ROOT.kBlue)
  h5.Draw()
  h7.Draw("same")
  c3.BuildLegend(0.2393484,0.5348432,0.6190476,0.7456446)
  c3.SaveAs("topMass2_%s_%s_Norm.png"%(input,jetType))
  c3.SaveAs("plotCode/topMass2_%s_%s_Norm.C"%(input,jetType))

  c4 = makeCanvas("top2Mass_btaggedOR_%s_%s"%(input,jetType))
  h6.Fit(fitmodel,"S")
  h6.SetLineColor(ROOT.kRed)
  h6.SetMarkerColor(ROOT.kRed)
  h8.Fit(fitmodel,"S")
  h8.SetLineColor(ROOT.kBlue)
  h8.SetMarkerColor(ROOT.kBlue)
  h6.Draw()
  h8.Draw("same")
  c4.BuildLegend(0.2393484,0.5348432,0.6190476,0.7456446)
  c4.SaveAs("topMass2_btaggedOR_%s_%s.png"%(input,jetType))
  c4.SaveAs("plotCode/topMass2_btaggedOR_%s_%s.C"%(input,jetType))
  """

if __name__ == "__main__" :
  if len(sys.argv) != 1 :
    print "Wrong argument!"
    sys.exit(-1)
  else :
    inputFileList = ["TTLL.root","DY50.root","DY10to50.root"]
    for input in inputFileList :
      gROOT.Reset()
      file0 = TFile.Open(input)
      outFileName = "histOut_"+os.path.splitext(input)[0]+".root"  # [ histOut_TTLL.root, histOut_DY50.root, histOut_DY10to50.root ]
      outFile= TFile.Open(outFileName,"RECREATE")

      for jetType in ["Jet","GenJet"] :
        treeName = jetType+"Tree"
        treeNameCharge = jetType+"TreeCharge"
        nEventName = "nEvent_"+jetType
        tree = file0.Get(treeName)
        tree2 = file0.Get(treeNameCharge)
        nEvent = file0.Get(nEventName).GetBinContent(1)
        anaTree(tree,tree2, nEvent, os.path.splitext(input)[0], outFile, jetType)
      outFile.Close()
      file0.Close()
