#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys, copy
from CATTools.CatAnalyzer.histoHelper import *
#import DYestimation
ROOT.gROOT.SetBatch(True)


datalumi = 100
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

#topMassList = ['TT_powheg_mtop1665','TT_powheg_mtop1695','TT_powheg_mtop1755','TT_powheg_mtop1785','TTLL_powheg']
#topMassList = ['TT_powheg_mtop1695','TT_powheg_mtop1755','TTLL_powheg']

mcfilelist = ['TTLL_powheg','SingleTbar_tW', 'SingleTop_tW','DYJets', 'DYJets_10to50','WJets','WW','WZ','ZZ']
#mcfilelist = ['TTLL_powheg','SingleTbar_tW', 'SingleTop_tW','DYJets', 'DYJets_10to50']
rootfileDir = "./tupleOut_"

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))


#saving mc histos
mchistList = []

print "DataSet  &  CrossSection & num of expected events & nLep2 & zVeto & nJet2 & MET40 & nbjet1"
for i, mcname in enumerate(mcfilelist):
  data = findDataSet(mcname, datasets)
  scale = datalumi*data["xsec"]
  colour = data["colour"]
  title = data["title"]

  rfname = rootfileDir + mcname +".root"
  tfile = ROOT.TFile(rfname)
  nEventHist = tfile.Get("nEvent")
  wentries = nEventHist.GetBinContent(1)
  #print wentries
  scale = scale/wentries
  
  h1 = nEventHist.Clone() 
  h1.Scale(scale)
  
  row =  "%s & %s & "%(mcname, data["xsec"], )
 
  for idx, xbin in enumerate(range(h1.GetNbinsX())) :
    if ( idx ==1 ) : continue
    entry =  h1.GetBinContent( xbin+1 )
    row += "& %8.3f"%(entry)
  print row
    

  

