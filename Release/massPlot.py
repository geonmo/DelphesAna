#!/usr/bin/env python
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys, copy
from CATTools.CatAnalyzer.histoHelper import *
#import DYestimation
ROOT.gROOT.SetBatch(True)


datalumi = 100
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb

topMassList = ['TT_powheg_mtop1665','TT_powheg_mtop1695','TT_powheg_mtop1755','TT_powheg_mtop1785','TTLL_powheg']
#topMassList = ['TT_powheg_mtop1695','TT_powheg_mtop1755','TTLL_powheg']
mcfilelist = ['SingleTbar_tW', 'SingleTop_tW','DYJets', 'DYJets_10to50','WJets','WW','WZ','ZZ']
rootfileDir = "./tupleOut_"
#channel_name = ['Combined', 'MuEl', 'ElEl', 'MuMu']

datasets = json.load(open("%s/src/CATTools/CatAnalyzer/data/dataset/dataset.json" % os.environ['CMSSW_BASE']))

#defalts
step = 1
channel = 1
cut = '1'
weight = '1'
binning = [60, 20, 320]
plotvar = 'll_m'
x_name = 'mass [GeV]'
y_name = 'Events'
dolog = False
overflow = False
binNormalize = False
suffix = ''
#get input
try:
    opts, args = getopt.getopt(sys.argv[1:],"hdnoc:w:b:p:x:y:a:s:f:",["binNormalize","overflow","cut","weight","binning","plotvar","x_name","y_name","dolog","channel","step","suffix"])
except getopt.GetoptError:          
    print 'Usage : ./.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print 'Usage : ./topDraw.py -c <cut> -w <weight> -b <binning> -p <plotvar> -x <x_name> -y <y_name> -d <dolog> -f <suffix>'
        sys.exit()
    elif opt in ("-c", "--cut"):
        cut = arg
    elif opt in ("-a", "--channel"):
        channel = int(arg)
    elif opt in ("-s", "--step"):
        step = int(arg)
    elif opt in ("-w", "--weight"):
        weight = arg
    elif opt in ("-b", "--binning"):
        binning = eval(arg)
    elif opt in ("-p", "--plotvar"):
        plotvar = arg
    elif opt in ("-x", "--x_name"):
        x_name = arg
    elif opt in ("-y", "--y_name"):
        y_name = arg
    elif opt in ("-d", "--dolog"):
        dolog = True
    elif opt in ("-o", "--overflow"):
        overflow = True
    elif opt in ("-n", "--binNormalize"):
        binNormalize = True
    elif opt in ("-f", "--suffix"):
        suffix = "_"+arg

tname = "tree"

#cut define
tcut="(%s)"%(cut)

#namming
#x_name = "Dilepton channel "+x_name
if len(binning) <= 3:
    num = (binning[2]-binning[1])/float(binning[0])
    if num != 1:
        if x_name.endswith(']'):
            unit = "["+x_name.split('[')[1]
        else: unit = ""
        y_name = y_name + "/%g%s"%(num,unit)


#saving mc histos
mchistList = []
for i, mcname in enumerate(mcfilelist):
  data = findDataSet(mcname, datasets)
  scale = datalumi*data["xsec"]
  colour = data["colour"]
  title = data["title"]

  rfname = rootfileDir + mcname +".root"
  tfile = ROOT.TFile(rfname)
  wentries = tfile.Get("nEvent").GetBinContent(1)
  scale = scale/wentries
    
  mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
  mchist.SetLineColor(colour)
  mchist.SetFillColor(colour)
  mchistList.append(mchist)
  

#overflow
if overflow:
  if len(binning) == 3 : nbin = binning[0]
  else : nbin = len(binnin)-1
  for hist in mchistList:
    hist.SetBinContent(nbin, hist.GetBinContent(nbin+1))

#bin normalize
if binNormalize and len(binning)!=3:
  for hist in mchistList:
    for i in range(len(binning)):
      hist.SetBinContent(i, hist.GetBinContent(i)/hist.GetBinWidth(i))
      hist.SetBinError(i, hist.GetBinError(i)/hist.GetBinWidth(i))

hs_bkg = ROOT.THStack("bkg_hs","bkg_hs")
for hist in mchistList :
  hs_bkg.Add( hist)
#hs_bkg.Draw()

#bkgs = hs_bkg.GetStack().Last()
#bkgs.SetName("bkg")
#bkgs.Draw()
#bkgs.Write()
#print "bkg entries: ",bkgs.GetEntries()
matchingvar=''
matchingtype=''
matchingcut='&&(1)'
if '443' in cut :
  #matchingcut = "jpsi_delPtTrue<0.1&&jpsi_dRTrue<0.1&&abs(jpsi_isFromTop)==6"
  #matchingcut = "jpsi_delPtTrue<0.1&&jpsi_dRTrue<0.1&&abs(jpsi_isFromTop)==6"
  matchingvar='J/#psi'
  matchingtype='jpsi'
 
elif '421' in cut :
  #matchingcut = "&&(sv_delPtTrue<0.2&&sv_delPtTrue>0&&sv_dRTrue<0.2&&sv_dRTrue>0)"
  #matchingcut = "&&(sv_delPtTrue>0&&sv_dRTrue>0)"
  matchingvar='D0'
  matchingtype='d0'
elif '413' in cut :
  #matchingcut = "&&(sv_delPtTrue<0.2&&sv_delPtTrue>0&&sv_dRTrue<0.2&&sv_dRTrue>0)"
  #matchingcut = "&&(sv_delPtTrue>0&&sv_dRTrue>0)"
  matchingvar='D*'
  matchingtype='dstar'
print "matchingvar : ",matchingvar
outMassHist = ROOT.TFile.Open("invMass_%s_%s.root"%(matchingtype,suffix),"RECREATE")
#outMassHist.cd()

lsvHists=[]
for topMass in topMassList :
  if ( topMass.find("mtop") != -1 ) :  massValue = topMass.split("mtop")[-1]
  else : massValue = "nominal" 
  sum_hs =  hs_bkg.Clone()
  sum_hs.SetTitle("Signal + Bkg.")
  data = findDataSet(topMass, datasets)
  scale = datalumi*data["xsec"]
  colour = data["colour"]
  title = data["title"]

  rfname = rootfileDir + topMass +".root"
  tfile = ROOT.TFile(rfname)
  #wentries = data['nevt']
  nEvent = tfile.Get("nEvent")
  wentries = nEvent.GetBinContent(1)
  scale = scale/wentries
  print topMass, scale, wentries, colour, title
   
    
  mchist = makeTH1(rfname, tname, title, binning, plotvar, tcut, scale)
  mchist.SetLineColor(colour)
  mchist.SetFillColor(colour)
  mcTruthHist = makeTH1(rfname, tname, title+' %s Signal'%(matchingvar), binning, plotvar, tcut+matchingcut, scale)
  mcTruthHist.SetLineColor(906)
  mcTruthHist.SetFillColor(906)
  if ( matchingcut != '&&(1)') : 
    mchist.Add(mcTruthHist,-1)


  print "topmass hsit : ",mchist.Integral()
  #overflow
  if overflow:
    if len(binning) == 3 : nbin = binning[0]
    else : nbin = len(binnin)-1
    mchist.SetBinContent(nbin, mchist.GetBinContent(nbin+1))

  #bin normalize
  if binNormalize and len(binning)!=3:
    for i in range(len(binning)):
      mchist.SetBinContent(i, mchist.GetBinContent(i)/mchist.GetBinWidth(i))
      mchist.SetBinError(i, mchist.GetBinError(i)/mchist.GetBinWidth(i))
      mcTruthHist.SetBinContent(i, mcTruthHist.GetBinContent(i)/mcTruthHist.GetBinWidth(i))
      mcTruthHist.SetBinError(i, mcTruthHist.GetBinError(i)/mcTruthHist.GetBinWidth(i))

  sum_hs.Add( mchist )
  if ( matchingcut != '&&(1)') : 
    sum_hs.Add( mcTruthHist )

  masshist = sum_hs.GetStack().Last()
  #masshist.Sumw2()
  masshist.Draw()
  print masshist.GetEntries()
  masshist.SetName("invMass_%s"%(massValue))
  masshist.SetTitle("Invariant mass mtop_%s; M_{%s}  ;Entries/%f"%( massValue, matchingvar, masshist.GetBinWidth(1) ))
  if ( plotvar in ["lsv_mass"] ) :
    lsvHists.append(masshist.Clone())
    outMassHist.cd()
    masshist.Write()
    mchist.SetName("ttbar_mtop%s"%(massValue))
    mchist.Write()

  c1= ROOT.TCanvas("c1","c1")
  """
  if plotvar in ["sv_mass"] :
    sum_hs.SetTitle("Invariant mass; M_{%s}[GeV/c^2];Entries/%f"%( matchingvar, masshist.GetBinWidth(1) ))
  elif plotvar in ["lsv_mass"] :
    sum_hs.SetTitle("Invariant mass; M_{l+%s}[GeV/c^2];Entries/%f"%( matchingvar, masshist.GetBinWidth(1) ))
  elif plotvar in ["sv_diffmass"] :
    sum_hs.SetTitle("Invariant mass; diff. M;Entries/%f"%( masshist.GetBinWidth(1) ))
  """
  sum_hs.SetTitle("Invariant mass distribution;%s;Entries/%f"%(x_name,masshist.GetBinWidth(1)))

  sum_hs.Draw("hist")
  sum_hs.SetName("hstack_%s"%(massValue))
  sum_hs.Write()
  leg = ROOT.TLegend(0.7,0.75,0.95,0.9)
  leg.AddEntry(mchistList[0],mchistList[0].GetTitle())
  leg.AddEntry(mchistList[2],mchistList[2].GetTitle())
  leg.AddEntry(mchistList[4],mchistList[4].GetTitle())
  leg.AddEntry(mchistList[5],mchistList[5].GetTitle())
  leg.AddEntry(mchist,mchist.GetTitle())
  leg.Draw()
  #leg = c1.BuildLegend(0.7,0.4,0.9,0.7)
  c1.SaveAs("plot_%s_%s_%s.png"%(matchingtype,plotvar, topMass))

if ( len(lsvHists) != 0 ) : 
  c1 = ROOT.TCanvas("c1","c1")
  lsvHists[0].SetLineColor(ROOT.kBlue)
  lsvHists[0].SetMarkerColor(ROOT.kBlue)
  lsvHists[0].SetFillColor(0)
  lsvHists[1].SetLineColor(ROOT.kBlack)
  lsvHists[1].SetMarkerColor(ROOT.kBlack)
  lsvHists[1].SetFillColor(0)
  lsvHists[2].SetFillColor(0)
  lsvHists[2].SetMarkerColor(ROOT.kRed)
  #lsvHists[1].SetLineColor(ROOT.kYellow)
  for idx, lsvHist in enumerate(lsvHists) :
    if ( idx==0 ) : lsvHist.Draw("e")
    else : lsvHist.Draw("esame")
  c1.BuildLegend(0.4,0.2,0.8,0.4)
  c1.SaveAs("plot_lsv_%s_%s.png"%(matchingtype, plotvar))

outMassHist.Write()
outMassHist.Close()
