void fitTopMassNorm()
{

  gROOT->SetBatch(true);
  gSystem->Load("libRooFit");

  std::vector<const char*> inFiles;
  inFiles.push_back("histOut_TTLL.root"); 
  inFiles.push_back("histOut_DY10to50.root");
  inFiles.push_back("histOut_DY50.root");

  const char* inFileType[3] = {"TTLL","DY10to50","DY50"};
  int colors[4] = {kRed, kBlue, kBlue, kBlack};

  const char* histType[4] = { "","_btagged","_charged","_btagged_charged"};
  double xsec[3] = {87.31, 18610, 6025.2};  
  double lumi=50 * 1000;

  std::vector<std::vector<TH1F*> > GenJetHists;
  std::vector<std::vector<TH1F*> > JetHists;


  TH1F* blankHist = nullptr;
  for( int i= 0 ; i <3 ; i++) {
    TFile* inFile = TFile::Open(inFiles[i]);
    std::vector<TH1F*> hists_GenJet;
    std::vector<TH1F*> hists_Jet;
    for( int j= 0 ; j <4 ; j++) {
      TString histName = TString::Format("top_mass%s_Jet_Norm",histType[j]);
      TH1F* tempJet = (TH1F*)(inFile->Get( histName.Data() ))->Clone();
      if ( tempJet != nullptr) hists_Jet.push_back(tempJet); 
      TH1F* tempGenJet = (TH1F*)(inFile->Get( TString::Format("top_mass%s_GenJet_Norm",histType[j]).Data()))->Clone();
      if ( tempGenJet != nullptr) hists_GenJet.push_back(tempGenJet); 
    }
    GenJetHists.push_back(hists_GenJet);
    JetHists.push_back(hists_Jet);
  }
  for(int j=0 ; j< 4 ; j++) {
    THStack* ts = new THStack("ts","ts");
    TH1F* ttll     = (TH1F*)JetHists[0][j]->Clone();
    ttll->Scale(xsec[0]*lumi);
    ttll->SetFillColor( colors[0]);
    ttll->SetMarkerStyle(21);
    ttll->SetMarkerColor( colors[0]);
    //ttll->Rebin(20);
    ttll->SetTitle("t#bar{t} #rightarrow l^{+} l^{-}");
    ts->Add(ttll);
 
    TH1F* dy10to50 = (TH1F*)JetHists[1][j]->Clone();
    dy10to50->Scale(xsec[1]*lumi);
    dy10to50->SetFillColor( colors[1]);
    dy10to50->SetMarkerStyle(21);
    dy10to50->SetMarkerColor( colors[1]);
    //dy10to50->Rebin(20);
    dy10to50->SetTitle("DY");

    TH1F* dy50     = (TH1F*)JetHists[2][j]->Clone();
    dy50->Scale(xsec[2]*lumi);
    dy50->SetFillColor( colors[2]);
    dy50->SetMarkerStyle(21);
    dy50->SetMarkerColor( colors[2]);
    //dy50->Rebin(20);

    dy10to50->Add(dy50);
    ts->Add(dy10to50);

    TH1F* sumHist = (TH1F*)ts->GetStack()->Last();
    int nTotal = sumHist->Integral();
    
    RooRealVar x("mass","M_{top}",100,200);
    std::pair<float, float> range(make_pair(150,180)); 
    RooDataHist* dh = new RooDataHist("dh","data histogram",x,sumHist);

    auto BW_mean  = RooRealVar("BW_mean","BW_mean",172.5,100,200);
    auto BW_sigma = RooRealVar("BW_sigma","BW_sigma",10,0,100);
    auto nBW = RooRealVar("nBW","nBW",nTotal*0.7, 0, nTotal);
    auto BW_pdf = RooBreitWigner("sig","signal p.d.f",x,BW_mean, BW_sigma );

    RooFormulaVar mirrorX("neg_x","-mass", RooArgSet(x)); 
    RooRealVar landau_mean("lan_mean","Landau_mean",-160,-200,-100);
    RooRealVar landau_sigma("lan_sigma","Landau_sigma",2.5, 0,100);
    RooRealVar nLandau("nLan","nlan",nTotal*0.3, 0, nTotal);
    RooLandau  landau_pdf("lx","lx",mirrorX,landau_mean, landau_sigma);

    RooAddPdf* model;
    RooFitResult * fitResult;
    TPaveText* t1 = new TPaveText(0.15,0.5,0.45,0.8,"brNDC");
    model = new RooAddPdf("model","model",RooArgList(BW_pdf, landau_pdf),RooArgList(nBW, nLandau));
    fitResult = model->fitTo(*dh, RooFit::Extended(true), RooFit::Save(), RooFit::Range( range.first, range.second ));
    t1->AddText(TString::Format("Breit-Wigner ) Mean : %3.3f(+-%3.3f)", BW_mean.getVal(), BW_mean.getError()));
    t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", BW_sigma.getVal(), BW_sigma.getError()));
    t1->AddText(TString::Format("Landau ) Mean : %3.3f(+-%3.3f)", landau_mean.getVal(), landau_mean.getError()));
    t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", landau_sigma.getVal(), landau_sigma.getError()));
   
    t1->Draw();
    //ts->Draw("hist");
    RooPlot* top_mass_frame = x.frame();
    dh->plotOn(top_mass_frame);
    model->plotOn(top_mass_frame);
    model->plotOn(top_mass_frame, RooFit::Components(BW_pdf),RooFit::LineStyle(kDashed));
    model->plotOn(top_mass_frame, RooFit::Components(landau_pdf),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

    TCanvas* c1 = new TCanvas("c1","c1");
    top_mass_frame->addObject(t1);
    //top_mass_frame->addObject(sumHist);
    top_mass_frame->Draw();
    top_mass_frame->SetTitle(histType[j]);

    //c1->BuildLegend(0.2,0.5,0.5,0.65);
    c1->SaveAs( TString::Format("topMass_TTLL_DY_%s.png",histType[j]));
    c1->SaveAs( TString::Format("plotCode/topMass_TTLL_DY_%s.C",histType[j]));
  }    
}
  

