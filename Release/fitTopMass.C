void fitTopMass(const char * inFileName)
{

  gROOT->SetBatch(true);
  gSystem->Load("libRooFit");
  TFile* inFile = TFile::Open(inFileName);

  RooRealVar x("mass","M_{top}",100,200);
  using namespace std;
  vector<TH1F*> hists;
  vector<int> nTotal;
  std::vector<std::pair<float, float> >  range;

  const char* types[] = { "top1","top2","top"};
  for( auto type : types ) {
    hists.push_back (  (TH1F*)inFile->Get(TString::Format("%s_mass_GenJet",type)) ) ;
    hists.push_back (  (TH1F*)inFile->Get(TString::Format("%s_mass_Jet",type)) ) ;
    hists.push_back (  (TH1F*)inFile->Get(TString::Format("%s_mass_btagged_Jet",type)) ) ;
    hists.push_back (  (TH1F*)inFile->Get(TString::Format("%s_mass_charged_Jet",type)) ) ;
    hists.push_back (  (TH1F*)inFile->Get(TString::Format("%s_mass_btagged_charged_Jet",type)) ) ;
  }
  for( int i= 0 ; i<hists.size() ; i++) {
    if ( hists[i] == nullptr) std::cout<<"Error! Hist is nullptr."<<std::endl;
  }


  for( int i = 0 ; i < hists.size() ; i++) {
    nTotal.push_back(hists[i]->Integral());
  } 

  // For top1,
  for(int i= 0 ; i < 5 ; i++) {
    range.push_back( make_pair(150,180) ); 
  }
  // For top2,
  for(int i= 0 ; i < 5 ; i++) {
    range.push_back( make_pair(150,180) ); 
  }
  // For all top,
  for(int i= 0 ; i < 5 ; i++) {
    range.push_back( make_pair(150,180) ); 
  }
 
  for( int i= 0 ; i< hists.size() ; i++) {
    std::cout<<"hist["<<i<<"]"<<std::endl;
    RooDataHist* dh = new RooDataHist("dh","data histogram",x,hists[i]);

    auto gaus_mean = RooRealVar("gaus_mean","gaus_mean",172.5, 100,200);
    auto gaus_sigma = RooRealVar("gaus_sigma","gaus_sigma",1, 0,100);
    auto nGaus = RooRealVar("nGaus","nGaus",nTotal[i]*0.2, 0., nTotal[i]);
    auto gaus_pdf = RooGaussian("gaus_pdf","Gaussian p.d.f",x,gaus_mean, gaus_sigma);

    
    auto expo_tau = RooRealVar("exp_const","exp_const",0,-10000,10000);
    auto nExp = RooRealVar("nExp","nExp",1,0,nTotal[i]);
    auto expo_pdf = RooExponential("bkg","bkg p.d.f",x,expo_tau);

    auto CB_mean  = RooRealVar("CB_mean","CB_mean",172.5,100,200);
    auto CB_sigma = RooRealVar("CB_sigma","CB_sigma",5,0,100);
    auto CB_alpha = RooRealVar("CB_alpha","CB_alpha",1,0,100);
    auto CB_n     = RooRealVar("CB_n","CB_n",1,0,100);
    auto nCB = RooRealVar("nCB","nCB",nTotal[i]*0.8,0.,nTotal[i]);
    auto CB_pdf = RooCBShape("sig","signal p.d.f",x,CB_mean, CB_sigma, CB_alpha, CB_n );

    auto BW_mean  = RooRealVar("BW_mean","BW_mean",172.5,100,200);
    auto BW_sigma = RooRealVar("BW_sigma","BW_sigma",10,0,100);
    auto nBW = RooRealVar("nBW","nBW",nTotal[i]*0.7, 0, nTotal[i]);
    auto BW_pdf = RooBreitWigner("sig","signal p.d.f",x,BW_mean, BW_sigma );
  
    RooFormulaVar mirrorX("neg_x","-mass", RooArgSet(x)); 
    //RooRealVar landau_mean("lan_mean","Landau_mean",150.5,100,200);
    RooRealVar landau_mean("lan_mean","Landau_mean",-160,-200,-100);
    RooRealVar landau_sigma("lan_sigma","Landau_sigma",2.5, 0,100);
    RooRealVar nLandau("nLan","nlan",nTotal[i]*0.3, 0, nTotal[i]);
    RooLandau  landau_pdf("lx","lx",mirrorX,landau_mean, landau_sigma);
  
    //TF1* f1 = new TF1("landau","nLan*ROOT::Math::landau_pdf(-x,lan_sigma,-lan_mean)",100,300);
    //RooAbsPdf* invlandau_pdf = RooFit::bindPdf(f1,x,nLandau, landau_mean,landau_sigma);
    //RooGenericPdf invlandau_pdf("invLandau","ROOT::Math::landau_pdf(-x,lan_sigma,-lan_mean)",RooArgSet(x,landau_sigma,landau_mean)); 
 
    //sig_pdf = cb_pdf;
    //bkg_pdf = bw_pdf;
    
    //nsig = nCB;
    //nbkg = nBW;

    //RooFFTConvPdf model("lxg","CB (x) gauss",x,CB_pdf,gaus_pdf);
    //RooFFTConvPdf model("lxg","CB (x) BW",x,CB_pdf,BW_pdf);
    //auto model = RooAddPdf("model","model",RooArgList(landau),RooArgList(nlandau));

    //auto model = RooAddPdf("model","model",RooArgList(sig_pdf,  bkg_pdf),RooArgList(nsig, nbkg));
    //auto model = RooAddPdf("model","model",RooArgList(CB_pdf,  gaus_pdf),RooArgList(nCB, nGaus));
    RooAddPdf* model;
    RooFitResult * fitResult;
    TPaveText* t1 = new TPaveText(0.15,0.5,0.45,0.8,"brNDC");
    //auto model = RooAddPdf("model","model",RooArgList(CB_pdf,  expo_pdf),RooArgList(nCB, nExp));
    //auto model = RooAddPdf("model","model",RooArgList(CB_pdf,  BW_pdf),RooArgList(nCB, nBW));
    //auto model = RooAddPdf("model","model",RooArgList(CB_pdf),RooArgList(nCB));
    //auto model = RooAddPdf("model","model",RooArgList(cb_pdf,  gaus_pdf),RooArgList(nCB, nGaus));
    //auto model = sig_pdf;
    std::cout<<"Hello"<<std::endl;

    if ( i< 5 ) {
      /*
      model = new RooAddPdf("model","model",RooArgList(BW_pdf),RooArgList(nBW));
      fitResult = model->fitTo(*dh, RooFit::Extended(true), RooFit::Save(), RooFit::Range( range[i].first, range[i].second ));
      t1->AddText(TString::Format("Breit-Wigner ) Mean : %3.3f(+-%3.3f)", BW_mean.getVal(), BW_mean.getError()));
      t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", BW_sigma.getVal(), BW_sigma.getError()));
      */
      model = new RooAddPdf("model","model",RooArgList(BW_pdf, landau_pdf),RooArgList(nBW, nLandau));
      fitResult = model->fitTo(*dh, RooFit::Extended(true), RooFit::Save(), RooFit::Range( range[i].first, range[i].second ));
      t1->AddText(TString::Format("Breit-Wigner ) Mean : %3.3f(+-%3.3f)", BW_mean.getVal(), BW_mean.getError()));
      t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", BW_sigma.getVal(), BW_sigma.getError()));
      t1->AddText(TString::Format("Landau ) Mean : %3.3f(+-%3.3f)", landau_mean.getVal(), landau_mean.getError()));
      t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", landau_sigma.getVal(), landau_sigma.getError()));

    }
    if ( i>=5 && i<10 ) {
      model = new RooAddPdf("model","model",RooArgList(BW_pdf, landau_pdf),RooArgList(nBW, nLandau));
      fitResult = model->fitTo(*dh, RooFit::Extended(true), RooFit::Save(), RooFit::Range( range[i].first, range[i].second ));
      t1->AddText(TString::Format("Breit-Wigner ) Mean : %3.3f(+-%3.3f)", BW_mean.getVal(), BW_mean.getError()));
      t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", BW_sigma.getVal(), BW_sigma.getError()));
      t1->AddText(TString::Format("Landau ) Mean : %3.3f(+-%3.3f)", landau_mean.getVal(), landau_mean.getError()));
      t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", landau_sigma.getVal(), landau_sigma.getError()));

    }
    if ( i>=10 ) {
      model = new RooAddPdf("model","model",RooArgList(BW_pdf, landau_pdf),RooArgList(nBW, nLandau));
      fitResult = model->fitTo(*dh, RooFit::Extended(true), RooFit::Save(), RooFit::Range( range[i].first, range[i].second ));
      t1->AddText(TString::Format("Breit-Wigner ) Mean : %3.3f(+-%3.3f)", BW_mean.getVal(), BW_mean.getError()));
      t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", BW_sigma.getVal(), BW_sigma.getError()));
      t1->AddText(TString::Format("Landau ) Mean : %3.3f(+-%3.3f)", landau_mean.getVal(), landau_mean.getError()));
      t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", landau_sigma.getVal(), landau_sigma.getError()));
    }

    //t1->AddText(TString::Format("Breit-Wigner ) Mean : %3.3f(+-%3.3f)", BW_mean.getVal(), BW_mean.getError()));
    //t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)", BW_sigma.getVal(), BW_sigma.getError()));
    //t1->AddText(TString::Format("CrystalBall ) Mean : %3.3f(+-%3.3f)", CB_mean.getVal(), CB_mean.getError()));
    //t1->AddText(TString::Format("#sigma : %0.3f(+-%0.3f)  #alpha : %0.2f ", CB_sigma.getVal(), CB_sigma.getError(), CB_alpha.getVal()));
    t1->Draw();

    RooPlot* top_mass_frame = x.frame();
    dh->plotOn(top_mass_frame);

    model->plotOn(top_mass_frame);
    model->plotOn(top_mass_frame, RooFit::Components(BW_pdf),RooFit::LineStyle(kDashed));
    //model.plotOn(top_mass_frame, RooFit::Components(CB_pdf),RooFit::LineStyle(kDashed));
    //model.plotOn(top_mass_frame, RooFit::Components(gaus_pdf),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
    model->plotOn(top_mass_frame, RooFit::Components(landau_pdf),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

    TString canvas_name = TString::Format("FitResult_%s",hists[i]->GetName());
    TCanvas* c1 = new TCanvas(canvas_name.Data(),canvas_name.Data(),800,600);
    top_mass_frame->addObject(t1);
    top_mass_frame->Draw();
    top_mass_frame->SetTitle(canvas_name.Data());


    c1->SaveAs( (TString(c1->GetName())+".png").Data());
    c1->SaveAs( (TString("plotCode/")+TString(c1->GetName()).Strip()+".C").Data());
    delete c1;
    delete dh;
  }
}


