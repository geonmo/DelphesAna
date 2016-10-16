#!/usr/bin/env python

from ROOT import *
import json, os, sys, math, getopt 
from CATTools.CatAnalyzer.histoHelper import *
from array import array
gROOT.SetBatch(True)

gStyle.SetPalette(1)


def anaTree(inFile, num) :
  tags = [ "" ,"_btagged","_charged","_btagged_charged"]
  inFileName = [ "TTLL","DY50","DY10to50"]
  #for jetType in ["Jet","GenJet"] :
  c1 = TCanvas("pieChart","pieChart",600,600)
  c1.Divide(2,2)
  idx=0
  for jetType in ["Jet"] :
    for type in tags :
      idx+=1
      c1.cd(idx)
      coHist = inFile.Get("correct_Pair_Lep_and_Jet%s_%s"%(type,jetType))
      nXbin = coHist.GetNbinsX()
      coHist.SetTitle( coHist.GetTitle()+"_%s"%(inFileName[num]) )
      print coHist.GetTitle()
      if ( coHist.GetEntries() ==0 ) : continue
      coHist.Scale(1./ coHist.Integral()*100)
      """
      binCenter=[]
      binValue= []
      binColor = []
      for x in range(nXbin) :
        binCenter.append( int(coHist.GetBinCenter(x+1)))
        binValue.append( coHist.GetBinContent(x+1) )
        binColor.append( x+1 )
        
      piePlot = TPie("pie",coHist.GetTitle(), len(binCenter), binValue, binColor )
      """
      piePlot = TPie(coHist)
      piePlot.SetEntryRadiusOffset(nXbin+1 , 1)
      piePlot.MakeLegend()
      piePlot.SetLabelFormat("%txt(%perc)")
      for x in range(nXbin) :
        piePlot.SetEntryFillColor(x, x%8+2)
        piePlot.SetEntryLabel(x,"lf") 
      piePlot.SetEntryLabel( coHist.FindBin(0)-1, "Fake")
      piePlot.SetEntryLabel( coHist.FindBin(-5)-1, "Wrong")
      piePlot.SetEntryLabel( coHist.FindBin(5)-1, "Correct")
      piePlot.SetLabelsOffset(-0.2)
      piePlot.Draw("nol <")
  c1.SaveAs("ChartPie_%s.png"%(inFileName[num]))

if __name__ == "__main__" :
  if len(sys.argv) != 1 :
    print "Wrong argument!"
    sys.exit(-1)
  else :
    inputFileList = ["histOut_TTLL.root","histOut_DY50.root","histOut_DY10to50.root"]
    for num, input in enumerate(inputFileList) :
      gROOT.Reset()
      file0 = TFile.Open(input)
      anaTree( file0,num )

