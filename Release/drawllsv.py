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
norm_graphs=[]
for infile in infiles: 
  type=infile.split('invMass_')[-1].split('_')[0]
  sv = ""
  print type
  if ( type == "jpsi") : sv="J/#psi" 
  elif ( type == "d0") : sv="(D^{0}+l)" 
  elif ( type == "dstar") : sv="(D^{*}+l)" 
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
    #tf1 = TF1("f1","gaus",40,80)
           
    mean = clone_hist.GetMean()#fitresult.GetParameter(1)
    
    tf0 = None
    tf1 = None
    tf2 = None
    tf3 = None
    if   ( type=="jpsi") :  
      tf0 = TF1("f0","gaus",20,100)
    elif ( type=="d0") :    
      tf0 = TF1("f0","gaus",10,70)
    elif ( type=="dstar") : 
      tf0 = TF1("f0","gaus",20,80)


    total = clone_hist.GetEntries()
    #isComposite = False
    isComposite = ( type=="jpsi" or type=="dstar" or type=="d0")
  
    if isComposite :
      if ( type == "jpsi") : 
        tf2 = TF1("f2","gaus",20,100)
      if ( type == "d0") :    
        tf2 = TF1("f2","gaus",10,70)
      if ( type == "dstar") : 
        tf2 = TF1("f2","gaus",20,80)

      #clone_hist.Fit("f0","R")
      #fitresult = TVirtualFitter.GetFitter()

      nGaus = total
      """
      nLandau = total
      if ( type == "jpsi") : 
        nGaus = total*0.5
        nLandau = total*0.5
      if ( type == "d0") : 
        nGaus = total*0.8
        nLandau = total*0.2
      if ( type == "dstar") : 
        nGaus = total*0.01
        nLandau = total*0.99

      if ( type=="d0" and "1665" in mass) :
        nLandau = total*10
      if ( type=="d0" and "nominal" in mass) :
        nLandau = total
      """


      #GausMean = fitresult.GetParameter(1)
      #tf2.SetParameter(0,nGaus)
      #tf2.SetParameter(1,GausMean )
      #tf2.SetParameter(2,fitresult.GetParameter(2) )



      #tf2.SetParLimits(0,0,clone_hist.GetEntries() )
      #tf2.SetParLimits(3,0,clone_hist.GetEntries() )
      #tf2.SetParLimits(1,GausMean-7, GausMean+7)

    tf2.SetLineWidth(2)
    clone_hist.Fit("f2","R")
    fitresult = TVirtualFitter.GetFitter()


    c1 = makeCanvas(mass,False)
    c1 = setMargins(c1, False)
    c1.SetTopMargin(30)
    setDefAxis( clone_hist.GetXaxis(), "M_{l+%s} [GeV/c^2]"%sv,1)
    setDefAxis( clone_hist.GetYaxis(), "Entries",1)
    clone_hist.SetTitle("Invariant Mass of Lepton + Secondary Vertex")
    gStyle.SetTitleFontSize(0.2);
    clone_hist.Draw()
    if isComposite : 
      tf0.SetLineColor(ROOT.kBlue)
      #tf2.SetLineColor(ROOT.kRed)
      tf0.Draw("same")
      #tf2.Draw("same")
    c1.Update()
    st = clone_hist.FindObject("stats")
    st.SetY1NDC(0.6)
    st.SetY2NDC(0.9)
    #CMS_lumi.CMS_lumi( c1, 0, 11 )
    c1.Modified()
    c1.Update()


    fitresult = TVirtualFitter.GetFitter()
    sig = fitresult.GetParameter(1)
    error = fitresult.GetParError(1)

    print "sig ", sig 
    mean_gauss.append(sig)
    error_gauss.append(error)

    x_error_gauss.append(0)
    c1.SaveAs(mass+"_%s_landau.png"%(type))
    c1.SaveAs(mass+"_%s_landau.eps"%(type))

  #mass_various = [ 166.5, 169.5, 171.5, 173.5, 175.5, 178.5]
  mass_various = [ 166.5, 169.5, 175.5, 178.5, 172.5]
  #mass_various = [ 166.5, 169.5, 175.5, 178.5]
  x_array = ar.array("f", mass_various)
  x_error_array = ar.array("f",x_error_gauss)
  error_array = ar.array("f",error_gauss)


  y_array = ar.array("f",mean_gauss)
  print "nb. x_array : ",len(x_array)
  h1 = ROOT.TGraphErrors( len(x_array), x_array, y_array, x_error_array, error_array)
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

  normalized_hist = file0.Get("invMass_nominal").Clone()
  for x in range( normalized_hist.GetNbinsX()) : 
    normalized_hist.SetBinError( x+1, TMath.Sqrt( normalized_hist.GetBinContent(x+1)))
  if ( type == "jpsi") :
    tf3 = TF1("f3","gaus",20,100)
  if ( type == "d0") :
    tf3 = TF1("f3","gaus",10,70)
  if ( type == "dstar") :
    tf3 = TF1("f3","gaus",20,80)
  normalized_hist.Fit(tf3,"RS")

  fitresult =  TVirtualFitter.GetFitter()
  sig = fitresult.GetParameter(1)
  error = fitresult.GetParError(1)
  norm_graph = ROOT.TGraphErrors(1,ar.array("f",[172.5]),ar.array("f",[sig]),ar.array("f",[0]),ar.array("f",[error]))
  norm_graph.SetMarkerSize(3)
  norm_graph.SetMarkerStyle(25)
  norm_graph.SetMarkerColor(ROOT.kBlack)

  norm_graphs.append(norm_graph)
  c2 = TCanvas("c2","c2")
  normalized_hist.Draw()
  c2.SaveAs("normHist%s.png"%(type))
  c2.SaveAs("normHist%s.eps"%(type))
  print "norm (172.5) y value :  ",sig,"  y error : ", error




c1 = makeCanvas("M_{t} vs M_{l+sv}",False)
leg = TLegend(0.25,0.67, 0.75,0.92)
leg.SetTextSize(0.05)
leg.SetBorderSize(0)
mg = TMultiGraph();
mg.SetTitle("M_{t} vs Invariant mass of (Lepton + Secondary Vertex); M_{t} (GeV/c^{2}) ; Invariant M_{l+SV}(GeV/c^{2})")

colors =[ ROOT.kRed, ROOT.kGreen+3, ROOT.kBlue  ]
markers = [ 21, 22, 23]
styles  = [ 3335,3353, 3244]
y0s=[]
y0Errors=[]
slopes=[]
slopeErrors=[]
for idx, sv_graph in enumerate(sv_mgs) :
  sv_graph.SetLineColor(colors[idx])
  sv_graph.SetMarkerColor(colors[idx])
  sv_graph.SetMarkerStyle(markers[idx])

  tf1 = TF1("linear","pol1")
  tf2 = TF1("linear2","pol1")
  tf3 = TF1("linear3","pol1")
  #tf1.FixParameter(0,0)
  #tf2.FixParameter(0,0)
  #tf3.FixParameter(0,0)
  tf1.SetLineWidth(2)
  sv_graph.Fit(tf1)
  fitresult2 = TVirtualFitter.GetFitter()
  y_0  = fitresult2.GetParameter(0)
  y_0Error  = fitresult2.GetParError(0)
  slope    = fitresult2.GetParameter(1)
  slopeError    = fitresult2.GetParError(1)
  
  print slope,slopeError, slope+slopeError
  tf2.FixParameter(0, y_0-y_0Error)
  tf2.FixParameter(1, slope)
  tf2.SetLineColor(colors[idx])
  #tf2.SetFillStyle(styles[idx])
  #tf2.SetFillColor(colors[idx])
  tf2.SetLineStyle(2)
  #sv_graph.Fit(tf2,"OB+")
  tf3.FixParameter(0, y_0+y_0Error)
  tf3.FixParameter(1, slope)
  tf3.SetLineColor(colors[idx])
  #tf3.SetFillStyle(styles[idx])
  #tf3.SetFillColor(colors[idx])
  tf3.SetLineStyle(2)
  #sv_graph.Fit(tf3,"OB+")

  y0s.append(y_0)
  y0Errors.append(y_0Error)
  slopes.append(slope)
  slopeErrors.append(slopeError)
 
  sv_graph.SetFillColor(0) 
  mg.Add(sv_graph)
  leg.AddEntry( sv_graph, "M_{l+%s}= (%.2f#pm%.2f) M_{t} + (%.1f #pm %.1f)"%(sv_vars[idx], slope,slopeError, y_0, y_0Error))
  leg.SetTextFont(12)
  

gStyle.SetOptStat(0)
gStyle.SetOptFit(0)


#mg.Add(norm_graphs[0])
#mg.Add(norm_graphs[1])
#mg.Add(norm_graphs[2])
mg.Draw("AP")

hh1 = mg.GetHistogram()
hh1.GetXaxis().SetLimits(160,185)
#hh1.GetXaxis().SetRangeUser(160,180)


mg.GetYaxis().SetRangeUser( 30,90)

setDefAxis( mg.GetXaxis(), "M_{top} (GeV/c^{2})",1)
setDefAxis( mg.GetYaxis(), "M^{inv.}_{l+sv}(GeV/c^{2})",1)

mg.GetXaxis().SetTitleFont(12)
mg.GetYaxis().SetTitleFont(12)

leg.Draw()
"""
norm_graphs[0].Draw("sameAP")
norm_graphs[1].Draw("sameAP")
norm_graphs[2].Draw("sameAP")
"""

pavet = TPaveText()
pavet.SetBorderSize(0)
pavet.SetFillColor(0)
"""
for x in range(len(sv_mgs)) :
  pavet.AddText("M_{l+%s}= %.2f#pm%.1f M_{t} + %.1f #pm %.1f"%(sv_vars[x], slopes[x],slopeErrors[x], y0s[x], y0Errors[x]))
#pavet.AddText("Nominal M_{top} : %3.2f #pm %.2f GeV/c^2 (Stat.)"%(nominal_mass,nominal_mass_error))
pavet.SetTextSize(0.05)
pavet.SetX1NDC(0.55)
pavet.SetX2NDC(0.9)
pavet.SetY1NDC(0.75)
pavet.SetY2NDC(0.9)
pavet.Draw()
"""

st = h1.FindObject("stats")
#st.SetY1NDC(0.3)
#st.SetY2NDC(0.4)


c1.Modified()
c1.Update()
c1.SaveAs("top_mass_landau.png")
c1.SaveAs("top_mass_landau.eps")

"""
CMS_lumi.lumiTextSize = 0.7
CMS_lumi.cmsTextSize = 0.9
CMS_lumi.CMS_lumi( c1,0,11)
"""


