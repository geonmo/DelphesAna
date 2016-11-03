#!/usr/bin/env python
import glob
from ROOT import *
import array as ar
import sys
from CATTools.CatAnalyzer.histoHelper import *
import ROOT, CATTools.CatAnalyzer.CMS_lumi, json, os, getopt, sys
import ROOT.TRandom3 as ran3
import glob
random = ran3.TRandom3()
random.SetSeed(1)


hist_for_PE_mean = TH1F("dist_mean_pe","Distribution of PE's mean value",120,0,120)
hist_for_PE_err  = TH1F("dist_err_pe","Distribution of PE's mean err value",100,0,50)


if ( len(sys.argv) != 1 ) : 
  sys.exit(-1)

#infiles = glob.glob("invMass_*.root")
infiles =["invMass_jpsi_.root","invMass_d0_.root","invMass_dstar_.root"] 
"""
types=[]
for filetype in infiles :
  types.append(filetype.split('invMass_')[-1].split('__')[0])
print types
"""
datalumi = 100
CMS_lumi.lumi_sqrtS = "%.1f fb^{-1}, #sqrt{s} = 13 TeV"%(datalumi)
datalumi = datalumi*1000 # due to fb
CMS_lumi.extraText   = "Private work"

CMS_lumi.relPosX = 0.06
CMS_lumi.cmsTextSize = 0.6

sv_mgs=[]
sv_vars=[]
for infile in infiles: 
  type=infile.split('invMass_')[-1].split('_')[0]
  sv = ""
  print type
  if ( type == "jpsi") : sv="J/#psi" 
  elif ( type == "d0") : sv="D0" 
  elif ( type == "dstar") : sv="D*(2010)" 
  sv_vars.append(sv)
  massList = [  
    "invMass_1665",
    "invMass_1695",
    #"invMass_1715",
    #"invMass_1735",
    "invMass_1755",
    "invMass_1785",
    "invMass_nominal"
  ]

  mean_gauss=[]
  error_gauss = []
  x_error_gauss=[]

  file0 = TFile(infile)

  for mass in massList :
    hist1= file0.Get(mass)
    #hist1.SetStats(False)
    clone_hist = hist1.Clone()
    tf1 = TF1("f1","landau",5,245)
    #tf1 = TF1("f1","gaus",40,80)
           
    clone_hist.Fit(tf1)
    c1 = makeCanvas(mass,False)
    c1 = setMargins(c1, False)
    c1.SetTopMargin(30)
    setDefAxis( clone_hist.GetXaxis(), "M_{l+%s} [GeV/c^2]"%sv,1)
    setDefAxis( clone_hist.GetYaxis(), "Entries",1)
    clone_hist.SetTitle("Invariant Mass of Lepton + Secondary Vertex")
    gStyle.SetTitleFontSize(0.2);
    clone_hist.Draw()
    leg = TLegend(0.4,0.15, 0.7,0.25)
    if ( mass != "invMass_nominal") : leg.AddEntry( clone_hist, "M_{top}=%.1fGeV/c^2"%(float(mass.split("invMass_")[-1])/10.)  )
    else : leg.AddEntry(clone_hist, "M_{top} nominal")
    leg.SetTextSize(0.03)
    leg.Draw()
    c1.Update()
    st = clone_hist.FindObject("stats")
    st.SetY1NDC(0.8)
    st.SetY2NDC(0.9)
    CMS_lumi.CMS_lumi( c1, 0, 11 )
    c1.Modified()
    c1.Update()

    fitresult = TVirtualFitter.GetFitter()
    sig = fitresult.GetParameter(1)
    error = fitresult.GetParError(1)

    print "sig-Highest point ", sig , "-",  tf1.GetMaximumX()
    mean_gauss.append(sig)
    #mean_gauss.append(tf1.GetMaximumX())
    error_gauss.append(error)

    x_error_gauss.append(0)
    c1.SaveAs(mass+"_%s_landau.png"%(type))

  #mass_various = [ 166.5, 169.5, 171.5, 173.5, 175.5, 178.5]
  mass_various = [ 166.5, 169.5, 175.5, 178.5, 172.5]
  x_array = ar.array("f", mass_various)
  x_error_array = ar.array("f",x_error_gauss)
  error_array = ar.array("f",error_gauss)


  y_array = ar.array("f",mean_gauss)
  print "nb. x_array : ",len(x_array)
  h1 = TGraphErrors( len(x_array), x_array, y_array, x_error_array, error_array)
  h1.SetName("M_{t} MC samples")
  h1.SetMarkerStyle( 23)
  h1.SetMarkerColor ( kRed )
  h1.SetMarkerSize( 2 )
  #h1.Draw("ALP")
  sv_mgs.append(h1)

  """
  x_array = ar.array("f", mass_various[:-1])
  x_error_array = ar.array("f",x_error_gauss[:-1])
  error_array = ar.array("f",error_gauss[:-1])
  """

c1 = makeCanvas("M_{t} vs M_{l+sv}",False)
leg = TLegend(0.2,0.55, 0.6,0.75)
leg.SetTextSize(0.05)
leg.SetBorderSize(0)
mg = TMultiGraph();
mg.SetTitle("M_{t} vs Invariant mass of (Lepton + Secondary Vertex); M_{t} [GeV/c^{2}] ; Invariant M_{l+SV}[GeV/c^{2}]")

colors =[ ROOT.kRed, ROOT.kGreen+3, ROOT.kBlue  ]
markers = [ 21, 22, 23]
y0s=[]
slopes=[]
slopeErrors=[]
for idx, sv_graph in enumerate(sv_mgs) :
  sv_graph.SetLineColor(colors[idx])
  sv_graph.SetMarkerColor(colors[idx])
  sv_graph.SetMarkerStyle(markers[idx])

  tf1 = TF1("linear","pol1")
  tf2 = TF1("linear2","pol1")
  tf3 = TF1("linear3","pol1")
  tf1.FixParameter(0,0)
  tf2.FixParameter(0,0)
  tf3.FixParameter(0,0)
  tf1.SetLineWidth(2)
  sv_graph.Fit(tf1)
  fitresult2 = TVirtualFitter.GetFitter()
  y_0  = fitresult2.GetParameter(0)
  slope    = fitresult2.GetParameter(1)
  slopeError    = fitresult2.GetParError(1)
  
  print slope,slopeError, slope+slopeError
  tf2.FixParameter(1, slope+slopeError)
  tf2.SetLineColor(ROOT.kBlue)
  tf2.SetLineStyle(2)
  sv_graph.Fit(tf2,"OB+")
  tf3.FixParameter(1, slope-slopeError)
  tf3.SetLineColor(ROOT.kBlue)
  tf3.SetLineStyle(2)
  sv_graph.Fit(tf3,"OB+")

  y0s.append(y_0)
  slopes.append(slope)
  slopeErrors.append(slopeError)
  mg.Add(sv_graph)
  leg.AddEntry( sv_graph, "M_{top} using %s"%(sv_vars[idx]),"p")

gStyle.SetOptStat(0)
gStyle.SetOptFit(0)
mg.Draw("AP")
mg.GetXaxis().SetRangeUser(160, 180)


mg.GetYaxis().SetRangeUser( 40,80)

setDefAxis( mg.GetXaxis(), "M_{top} [GeV/c^{2}]",1)
setDefAxis( mg.GetYaxis(), "Invariant M_{l+sv}[GeV/c^{2}]",1)


leg.Draw()


pavet = TPaveText()
pavet.SetBorderSize(0)
pavet.SetFillColor(0)
for x in range(len(sv_mgs)) :
  pavet.AddText("M_{l+%s}= %.2f M_{t} + %.1f"%(sv_vars[x], slopes[x],y0s[x]))
#pavet.AddText("Nominal M_{top} : %3.2f #pm %.2f GeV/c^2 (Stat.)"%(nominal_mass,nominal_mass_error))
pavet.SetTextSize(0.05)
pavet.SetX1NDC(0.55)
pavet.SetX2NDC(0.9)
pavet.SetY1NDC(0.75)
pavet.SetY2NDC(0.9)

pavet.Draw()

st = h1.FindObject("stats")
#st.SetY1NDC(0.3)
#st.SetY2NDC(0.4)


c1.Modified()
c1.Update()
c1.SaveAs("top_mass_landau.png")

"""
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.cmsTextSize = 0.9
CMS_lumi.CMS_lumi( c1,0,11)
"""
