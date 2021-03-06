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



def anaTree( tree ) :
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

  h0_p1 = getTH1("Jet1Charge_from_b_quark",chargeBin,tree,"bjet_charge[0]","bjet_partonPdgId[0]==5")
  h0_p2 = getTH1("Jet2Charge_from_b_quark",chargeBin,tree,"bjet_charge[1]","bjet_partonPdgId[1]==5")
  h0 = h0_p1+h0_p2;

  h0_1_p1 = getTH1("Jet1Charge_from_bbar_quark",chargeBin,tree,"bjet_charge[0]","bjet_partonPdgId[0]==-5")
  h0_1_p2 = getTH1("Jet2Charge_from_bbar_quark",chargeBin,tree,"bjet_charge[1]","bjet_partonPdgId[1]==-5")
  h0_1 = h0_1_p1+h0_1_p2;

  h0_2_p1 = getTH1("Jet1Charge_from_non_b_quark",chargeBin,tree,"bjet_charge[0]","abs(bjet_partonPdgId[0])!=5")
  h0_2_p2 = getTH1("Jet2Charge_from_non_b_quark",chargeBin,tree,"bjet_charge[1]","abs(bjet_partonPdgId[1])!=5")
  h0_2 = h0_2_p1+h0_2_p2;


  print pq_up+diffCharge+matched_and
  print pq_up+sameCharge+matched_and
  

  h1 = getTH1("Quality_for_all_jets; Quality weight ; Entries", qualityBin, tree,"quality*1",pq_up)
  h1_1 = getTH1("Quality_with_largeCharge; Quality Weight ; Entries", qualityBin,tree,"quality*1",pq_up+largeCharge)
  h1_2 = getTH1("Quality_with_diffCharge", qualityBin,tree,"quality*1",pq_up+diffCharge)
  h1_3 = getTH1("Quality_with_sameCharge", qualityBin,tree,"quality*1",pq_up+sameCharge)

  h1_4 = getTH1("Quality_btagged", qualityBin,tree,"quality*1",pq_up+btagged)
  h1_5 = getTH1("Quality_btagged_with_largeCharge", qualityBin,tree,"quality*1",pq_up+largeCharge+btagged)
  h1_6 = getTH1("Quality_btagged_with_diffCharge", qualityBin,tree,"quality*1",pq_up+diffCharge+btagged)
  h1_7 = getTH1("Quality_btagged_with_sameCharge", qualityBin,tree,"quality*1",pq_up+sameCharge+btagged)

  h2 = getTH1("Efficiency_OR", qualityBin,tree,"quality*1",pq_up+matched_or)
  h2_1 = getTH1("Efficiency_OR_with_largeCharge", qualityBin,tree,"quality*1",pq_up+largeCharge+matched_or)
  h2_2 = getTH1("Efficiency_OR_with_diffCharge", qualityBin,tree,"quality*1",pq_up+diffCharge+matched_or)

  h3 = getTH1("Efficiency_AND", qualityBin,tree,"quality*1",pq_up+matched_and)
  h3_1 = getTH1("Efficiency_AND_with_largeCharge", qualityBin,tree,"quality*1",pq_up+largeCharge+matched_and)
  h3_2 = getTH1("Efficiency_AND_with_diffCharge", qualityBin,tree,"quality*1",pq_up+diffCharge+matched_and)
  h3_3 = getTH1("Efficiency_AND__with_sameCharge", qualityBin,tree,"quality*1",pq_up+sameCharge+matched_and)

  h3_4 = getTH1("Efficiency_AND_btagged", qualityBin,tree,"quality*1",pq_up+matched_and+btagged)
  h3_5 = getTH1("Efficiency_AND_btagged_with_largeCharge", qualityBin,tree,"quality*1",pq_up+largeCharge+matched_and+btagged)
  h3_6 = getTH1("Efficiency_AND_btagged_with_diffCharge", qualityBin,tree,"quality*1",pq_up+diffCharge+matched_and+btagged)
  h3_7 = getTH1("Efficiency_AND_btagged_with_sameCharge", qualityBin,tree,"quality*1",pq_up+sameCharge+matched_and+btagged)

  h4   = getTH1("j1_pdgId_lq ; pdgId * LeptonCharge ; Entries", lqBin , tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up)
  h4_1 = getTH1("j1_pdgId_lq_only_largeCharge",lqBin ,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up + largeCharge)
  h4_2 = getTH1("j1_pdgId_lq_only_diffCharge", lqBin ,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up + diffCharge)

  h5   = getTH1("j2_pdgId_lq ; pdgId * LeptonCharge ; Entries", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up)
  h5_1 = getTH1("j2_pdgId_lq_only_largeCharge", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up + largeCharge)
  h5_2 = getTH1("j2_pdgId_lq_only_diffCharge", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up + diffCharge)

  h6   = getTH1("j1_pdgId_lq_btagged ; pdgId * LeptonCharge ; Entries", lqBin,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up+btagged_jet1)
  h6_1 = getTH1("j1_pdgId_lq_btagged_only_large_charge", lqBin,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up + largeCharge+btagged_jet1)
  h6_2 = getTH1("j1_pdgId_lq_btagged_only_diffcharge", lqBin,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up + diffCharge+btagged_jet1)

  h7   = getTH1("j2_pdgId_lq_btagged ; pdgId * LeptonCharge ; Entries", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up+btagged_jet2)
  h7_1 = getTH1("j2_pdgId_lq_btagged_only_largeCharge", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up + largeCharge+btagged_jet2)
  h7_2 = getTH1("j2_pdgId_lq_btagged_only_diffCharge", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up + diffCharge+btagged_jet2)


  h8   = getTH1("diffCharge",chargeBin,tree,"(bjet_charge[0]-bjet_charge[1])",pq_up)
  h8_1 = getTH1("diffCharge_matched",chargeBin,tree,"(bjet_charge[0]-bjet_charge[1])",pq_up+matched_and)
  h8_2 = getTH1("diffCharge_btagged",chargeBin,tree,"(bjet_charge[0]-bjet_charge[1])",pq_up+btagged)

  h9   = getTH2("diffCharge2D",[chargeBin,chargeBin],tree,"bjet_charge[1]:bjet_charge[0]",pq_up)
  h9_1 = getTH2("diffCharge2D_matched",[chargeBin,chargeBin],tree,"bjet_charge[1]:bjet_charge[0]",pq_up+matched_and)

  h10   = getTH1("largeCharge",chargeBin,tree,"lep_charge[0]*bjet_charge[0]+lep_charge[1]*bjet_charge[1]",pq_up)
  h10_1 = getTH1("largeCharge_matched",chargeBin,tree,"lep_charge[0]*bjet_charge[0]+lep_charge[1]*bjet_charge[1]",pq_up+matched_and)
  h10_2 = getTH1("largeCharge_btagged",chargeBin,tree,"lep_charge[0]*bjet_charge[0]+lep_charge[1]*bjet_charge[1]",pq_up+btagged)

  h11   = getTH2("largeCharge2D",[chargeBin,chargeBin],tree,"lep_charge[0]*bjet_charge[0]:lep_charge[1]*bjet_charge[1]",pq_up)
  h11_1 = getTH2("largeCharge2D_matched",[chargeBin,chargeBin],tree,"lep_charge[0]*bjet_charge[0]:lep_charge[1]*bjet_charge[1]",pq_up+matched_and)

  h12   = getTH1("nTracks",[30,0,30],tree,"bjet_nCharged",pq_up)
  h12_1   = getTH1("nTracks_from_bquark",[30,0,30],tree,"bjet_nCharged[0]",pq_up+"&&(bjet_partonPdgId[0]==5)")
  h12_2   = getTH1("nTracks_from_b2",[30,0,30],tree,"bjet_nCharged[1]",pq_up+"&&(bjet_partonPdgId[1]==5)")
  h12_3   = getTH1("nTracks_from_bbar1",[30,0,30],tree,"bjet_nCharged[0]",pq_up+"&&(bjet_partonPdgId[0]==-5)")
  h12_4   = getTH1("nTracks_from_bbar2",[30,0,30],tree,"bjet_nCharged[1]",pq_up+"&&(bjet_partonPdgId[1]==-5)")
  h12_5 = h12_1+h12_2+h12_3+h12_4

  h13   = getTH1("j1_pdgId_lq_bquark ; pdgId * LeptonCharge ; Entries", lqBin,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up+matched_jet1)
  h13_1 = getTH1("j1_pdgId_lq_bquark_only_large_charge", lqBin,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up + largeCharge+matched_jet1)
  h13_2 = getTH1("j1_pdgId_lq_bquark_only_diffcharge", lqBin,tree, "bjet_partonPdgId[0]*lep_charge[0]" ,pq_up + diffCharge+matched_jet1)

  h14   = getTH1("j2_pdgId_lq_bquark ; pdgId * LeptonCharge ; Entries", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up+matched_jet2)
  h14_1 = getTH1("j2_pdgId_lq_bquark_only_largeCharge", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up + largeCharge+matched_jet2)
  h14_2 = getTH1("j2_pdgId_lq_bquark_only_diffCharge", lqBin,tree, "bjet_partonPdgId[1]*lep_charge[1]" ,pq_up + diffCharge+matched_jet2)

  h15   = getTH1("top_mass ; M_{top}[GeV/c^2] ; Entries", [100,160,180],tree, "top_mass" ,pq_up+btagged_or)

  c0 = makeCanvas("Jet_Charge_of_all_jets")
  set3HistsColor(h0,h0_1,h0_2)
  h0.SetFillColor(ROOT.kRed)
  h0_1.SetFillColor(ROOT.kBlue)
  h0_2.SetFillColor(ROOT.kBlack)
  h0.SetFillStyle(3004)
  h0_1.SetFillStyle(3005)
  h0_2.SetFillStyle(3006)
  #h0_2.Draw("hist")
  h0.SetTitle("Jet Charge Distribution from various parton ; Jet Charge x 100 ; Entries") 
  h0.Draw("hist")
  h0_1.Draw("histsame")
  c0.BuildLegend(0.6779449,0.7700348,0.9511278,0.9111498)
  c0.SaveAs("jetCharge.png")
  c0.SaveAs("plotCode/jetCharge.C")

  c0_1 = makeCanvas("JetCharge_of_all_jets_Norm")
  h0.Scale(1./h0.GetEntries())
  h0_1.Scale(1./h0_1.GetEntries())
  h0_2.Scale(1./h0_2.GetEntries())
  h0.SetTitle("Normalized Jet Charge Distribution from various parton ; Jet Charge x 100 ; Normalized Entries") 
  h0.Draw("hist")
  h0_1.Draw("histsame")
  h0_2.Draw("histsame")
  c0_1.BuildLegend(0.6779449,0.7700348,0.9511278,0.9111498)
  c0_1.SaveAs("jetChargeNorm.png")
  c0_1.SaveAs("plotCode/jetChargeNorm.C")
  

  c1 = makeCanvas("Quality_for_matched_OR_for_all")
  set3HistsColor( h2,h2_1,h2_2)

  h2.Draw() 
  h2_1.Draw("same") 
  h2_2.Draw("same") 
  c1.BuildLegend()
  c1.SaveAs("quality_dist_or_all.png")
  c1.SaveAs("plotCode/quality_dist_or_all.C")
  print "c1"

  c2 = makeCanvas("Eff_for_matched_OR_for_all")
  h2.Divide(h1)
  h2_1.Divide(h1_1)
  h2_2.Divide(h1_2)

  h2.Draw() 
  h2_1.Draw("same") 
  h2_2.Draw("same") 
  c2.BuildLegend(0.4072682,0.434669,0.7869674,0.6454704)
  c2.SaveAs("eff_dist_or_all.png")
  c2.SaveAs("plotCode/eff_dist_or_all.C")
  print "c2"

  c4 = makeCanvas("Eff_for_matched_AND_for_all")
  set3HistsColor( h3,h3_1,h3_2)
  h3.Divide(h1)
  h3_1.Divide(h1_1)
  h3_2.Divide(h1_2)
  h3_3.Divide(h1_3)

  h3_2.Draw() 
  h3.Draw("same" )
  h3_1.Draw("same") 
  h3_3.Draw("same") 
  c4.BuildLegend(0.3790727,0.6289199,0.7587719,0.8397213)
  c4.SaveAs("eff_dist_and_all.png")
  c4.SaveAs("plotCode/eff_dist_and_all.C")

  c4_1 = makeCanvas("Eff_for_matched_AND_btagged_for_all")
  set3HistsColor( h3_4,h3_5,h3_6)
  h3_4.Divide(h1_4)
  h3_5.Divide(h1_5)
  h3_6.Divide(h1_6)
  h3_7.Divide(h1_7)

  h3_4.Draw() 
  h3_5.Draw("same" )
  h3_6.Draw("same") 
  h3_7.Draw("same") 
  c4_1.BuildLegend(0.3790727,0.2926829,0.7587719,0.5034843)
  c4_1.SaveAs("eff_dist_btagged_AND_all.png")
  c4_1.SaveAs("plotCode/eff_dist_btagged_AND_all.C")



  c5 = makeCanvas("leptonCharge_for_jet1")
  set3HistsColor( h4,h4_1,h4_2)
  h4.Scale(1./h4.GetEntries())
  h4_1.Scale(1./h4_1.GetEntries())
  h4_2.Scale(1./h4_2.GetEntries())

  h4.SetMaximum(1.0)
  h4.Draw() 
  h4_1.Draw("same") 
  h4_2.Draw("same") 
  c5.BuildLegend()
  c5.SaveAs("jet1_leptonCharge.png")
  c5.SaveAs("plotCode/jet1_leptonCharge.C")

  c6 = makeCanvas("leptonCharge_for_jet2")
  set3HistsColor( h5,h5_1,h5_2)
  h5.Scale(1./h5.GetEntries())
  h5_1.Scale(1./h5_1.GetEntries())
  h5_2.Scale(1./h5_2.GetEntries())

  h5.SetMaximum(1.0)
  h5.Draw() 
  h5_1.Draw("same") 
  h5_2.Draw("same") 
  #c6.BuildLegend(0.3,0.3,0.8,0.5)
  c6.BuildLegend()
  c6.SaveAs("jet2_leptonCharge.png")
  c6.SaveAs("plotCode/jet2_leptonCharge.C")

  c7 = makeCanvas("leptonCharge_for_bjet1")
  set3HistsColor( h6,h6_1,h6_2)
  h6.Scale(1./h6.GetEntries())
  h6_1.Scale(1./h6_1.GetEntries())
  h6_2.Scale(1./h6_2.GetEntries())

  h6.SetMaximum(1)
  h6.Draw() 
  h6_1.Draw("same") 
  h6_2.Draw("same") 

  c7.BuildLegend(0.3934837,0.4703833,0.773183,0.6811847)
  c7.SaveAs("bjet1_leptonCharge_btagged.png")
  c7.SaveAs("plotCode/bjet1_leptonCharge_btagged.C")




  c8 = makeCanvas("leptonCharge_for_bjet2")
  set3HistsColor( h7,h7_1,h7_2)
  h7.Scale(1./h7.GetEntries())
  h7_1.Scale(1./h7_1.GetEntries())
  h7_2.Scale(1./h7_2.GetEntries())

  h7.SetMaximum(1)
  h7.Draw() 
  h7_1.Draw("same") 
  h7_2.Draw("same") 
  c8.BuildLegend(0.3934837,0.4703833,0.773183,0.6811847)
  c8.SaveAs("bjet2_leptonCharge_btagged.png")
  c8.SaveAs("plotCode/bjet2_leptonCharge_btagged.C")

  c7_1 = makeCanvas("leptonCharge_for_bquark1")
  set3HistsColor( h13, h13_1, h13_2)
  h13.Scale(1./h13.GetEntries())
  h13_1.Scale(1./h13_1.GetEntries())
  h13_2.Scale(1./h13_2.GetEntries())
  h13.SetMaximum(1)
  h13.Draw() 
  h13_1.Draw("same") 
  h13_2.Draw("same") 

  c7_1.BuildLegend(0.3934837,0.4703833,0.773183,0.6811847)
  c7_1.SaveAs("bjet1_leptonCharge_bquark.png")
  c7_1.SaveAs("plotCode/bjet1_leptonCharge_bquark.C")

  c8_1 = makeCanvas("leptonCharge_for_bquark2")
  set3HistsColor( h14, h14_1, h14_2)
  h14.Scale(1./h14.GetEntries())
  h14_1.Scale(1./h14_1.GetEntries())
  h14_2.Scale(1./h14_2.GetEntries())
  h14.SetMaximum(1)
  h14.Draw() 
  h14_1.Draw("same") 
  h14_2.Draw("same") 

  c8_1.BuildLegend(0.3934837,0.4703833,0.773183,0.6811847)
  c8_1.SaveAs("bjet2_leptonCharge_bquark.png")
  c8_1.SaveAs("plotCode/bjet2_leptonCharge_bquark.C")

  c9_0 = makeCanvas("diffCharge")
  h8_1.SetLineColor(ROOT.kRed)
  h8.Draw()
  h8_1.Draw("same")
  c9_0.SaveAs("diffCharge.png")
  c9_0.SaveAs("plotCode/diffCharge.C")

  c9 = makeCanvas("diffChargeNorm")
  h8_1.Divide(h8)
  h8_1.Draw()
  h8_1.SetTitle("Truth Matching Rate vs JetCharge difference ; #Delta Charge ; Rate")
  c9.SaveAs("diffChargeNom_truthMatching.png")
  c9.SaveAs("plotCode/diffChargeNom_truthMatching.C")

  c10_1 = makeCanvas("diffCharge2D_truthMatching")
  #gStyle.SetPaintTextFormat("2g")
  h9_1.Draw("colztext")
  c10_1.SaveAs("diffCharge2D_truthMatching.png")
  c10_1.SaveAs("plotCode/diffCharge2D_truthMatching.C")

  c10_2 = makeCanvas("diffCharge2DEntries")
  #gStyle.SetPaintTextFormat("2g")
  h9.Draw("colztext")
  c10_2.SaveAs("diffCharge2DEntries.png")
  c10_2.SaveAs("plotCode/diffCharge2DEntries.C")



  c11 = makeCanvas("largeCharge")
  h10_1.SetLineColor(ROOT.kRed)
  h10.Draw()
  h10_1.Draw("same") 
  c11.SaveAs("largeCharge.png")
  c11.SaveAs("plotCode/largeCharge.C")

  c10 = makeCanvas("diffCharge2DNorm")
  gStyle.SetPaintTextFormat(".2g")
  h9_1.Divide(h9)
  h9_1.Draw("colztext")
  c10.SaveAs("diffCharge2DNorm_truthMatching.png")
  c10.SaveAs("plotCode/diffCharge2DNorm_truthMatching.C")
 
  c12 = makeCanvas("largeChargeNorm")
  h10_1.SetLineColor(ROOT.kRed)
  h10_1.Divide(h10)
  h10_1.Draw()
  h10_1.SetTitle("Truth Matching Rate vs leptonCharge*JetCharge ; lep_q*jet_q ; Rate")
  c12.SaveAs("largeChargeNorm.png")
  c12.SaveAs("plotCode/largeChargeNorm.C")

  c13 = makeCanvas("largeCharge2DEntries")
  h11.Draw("colztext")
  c13.SaveAs("largeCharge2DEntries.png")
  c13.SaveAs("plotCode/largeCharge2DEntries.C")
  
  c14 = makeCanvas("largeCharge2D_truthMatching")
  h11_1.Draw("colztext")
  c14.SaveAs("largeCharge2D_truthMatching.png")
  c14.SaveAs("plotCode/largeCharge2D_truthMatching.C")

  c15 = makeCanvas("largeCharge2DNorm_truthMatching")
  h11_1.Divide(h11)
  h11_1.SetTitle("Fraction of correct b-pair per 2D lepton * Jet charge; Lep1_q*jet1_q ; Lep2_q*Jet2_q")
  h11_1.Draw("colztext")
  c15.SaveAs("largeCharge2DNorm_truthMatching.png")
  c15.SaveAs("plotCode/largeCharge2DNorm_truthMatching.C")

  c16 = makeCanvas("number_of_tracks")
  c16.SetLogy()
  h12.SetLineColor(ROOT.kRed)
  h12.Draw()
  h12_5.SetLineColor(ROOT.kBlue)
  h12_5.Draw("same")
  c16.BuildLegend(0.2042607,0.1968641,0.5839599,0.4076655)
  c16.SaveAs("nTracks.png")
  c16.SaveAs("plotCode/nTracks.C")


if __name__ == "__main__" :
  if len(sys.argv) != 2 :
    print "Wrong argument!"
    sys.exit(-1)
  else :
    file0 = TFile.Open(sys.argv[1])
    tree = file0.Get("JetTree")
    anaTree(tree)
