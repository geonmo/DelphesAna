void fitTopMass(const char * inFileName)
{
  gROOT->SetBatch(true);
  gSystem->Load("libRooFit");
  TFile* inFile = TFile::Open(inFileName);

  RooRealVar x("mass","M_{top}",140,200);

  TH1F* h1 = (TH1F*)inFile->Get("top_mass_btagged_Jet");
  RooDataHist* dh = new RooDataHist("dh","data histogram",x,h1);

  int nTotal = h1->Integral();

  /*
  auto gaus_mean = RooRealVar("mean","mean",172.5, 100,200);
  auto gaus_sigma = RooRealVar("sigma","sigma",0,100);
  auto nsig = RooRealVar("nTrue","nTrue",0,nTotal);
  auto gaus_pdf = RooGaussian("sig","signal p.d.f",x,gaus_mean, gaus_sigma);
  auto sig_pdf = gaus_pdf;

  
  auto expo_tau = RooRealVar("exp_const","exp_const",0,-100,100);
  auto nbkg = RooRealVar("nFake","nFake",0.1,0.,1.0);
  auto expo_pdf = RooExponential("bkg","bkg p.d.f",x,expo_tau);

  auto bkg_pdf = expo_pdf; 
  */

  
  auto CB_mean  = RooRealVar("mean","mean",172.5,140,200);
  auto CB_sigma = RooRealVar("sigma","sigma",10);
  auto CB_alpha = RooRealVar("alpha","alpha",10);
  auto CB_n     = RooRealVar("n","n",5);
  auto nsig = RooRealVar("nTrue","nTrue",0.1,0.,1.0);
  
  auto sig_pdf = RooCBShape("sig","signal p.d.f",x,CB_mean, CB_sigma, CB_alpha, CB_n );
  /*
  auto BW_mean  = RooRealVar("mean","mean",172.5,160,180);
  auto BW_sigma = RooRealVar("sigma","sigma",10);
  auto nsig = RooRealVar("nTrue","nTrue",0.1,0.,1.0);
  auto sig_pdf = RooBreitWigner("sig","signal p.d.f",x,BW_mean, BW_sigma );
  */
   

  //auto model = RooAddPdf("model","model",RooArgList(sig_pdf));

  //auto model = RooAddPdf("model","model",RooArgList(sig_pdf,  bkg_pdf),RooArgList(nsig, nbkg));
  auto model = sig_pdf;


  auto* fitResult = model.fitTo(*dh, RooFit::Extended(true), RooFit::Save());
  if ( fitResult != nullptr) {
    //float fom = nsig.getVal() / TMath::Sqrt( nsig.getVal()+ nbkg.getVal());
    //std::cout<<"nSig : "<<nsig.getVal()<<"  nBkg : "<<nbkg.getVal()<<" Figure of merit : "<<fom<<std::endl; 
  }
  RooPlot* top_mass_frame = x.frame();
  dh->plotOn(top_mass_frame);

  model.plotOn(top_mass_frame);
  //model.plotOn(top_mass_frame, RooFit::Components(bkg_pdf),RooFit::LineStyle(kDashed));
  model.plotOn(top_mass_frame, RooFit::Components(sig_pdf),RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

  TString canvas_name = TString::Format("topMass_FitResult");
  TCanvas* c1 = new TCanvas(canvas_name.Data(),canvas_name.Data(),600,600);
  top_mass_frame->Draw();
  top_mass_frame->SetTitle(canvas_name.Data());
  c1->SaveAs( (TString(c1->GetName())+".png").Data());
  delete c1;
  delete dh;
}


