Double_t IntegrandBGPi(const double *x, const double *p)
{

  double x0   = x[0];
  double mass = 0.140;
  double pT   = p[0];
  double beta = p[1];
  double T    = p[2];
  double mT   = TMath::Sqrt(mass*mass + pT*pT);

  double n    = 1.;
  double rho0 = TMath::ATanH(beta*pow(x0,n));
  double a0   = pT*TMath::SinH(rho0)/T;
  double a1   = mT*TMath::CosH(rho0)/T;

  return x0*mT*TMath::BesselI0(a0)*TMath::BesselK1(a1);
}

Double_t StaticBGdNdPtPi(const double *x, const double *p)
{
  double pT   = x[0];
  double mass = 0.140;
  double beta = p[0];
  double T    = p[1];

  TF1 *fIntBG = 0;
  if(!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGPi, 0., 1., 3);

  fIntBG->SetParameters(pT, beta, T);
  //fIntBG->SetNpx(100000);
  return fIntBG->Integral(0.,1.) * p[2]; //p[2] norm
}
Double_t IntegrandBGP(const double *x, const double *p)
{

  double x0   = x[0];
  double mass = 0.938;
  double pT   = p[0];
  double beta = p[1];
  double T    = p[2];
  double mT   = TMath::Sqrt(mass*mass + pT*pT);

  double n    = 1.;
  double rho0 = TMath::ATanH(beta*pow(x0,n));
  double a0   = pT*TMath::SinH(rho0)/T;
  double a1   = mT*TMath::CosH(rho0)/T;

  return x0*mT*TMath::BesselI0(a0)*TMath::BesselK1(a1);
}

Double_t StaticBGdNdPtP(const double *x, const double *p)
{
  double pT   = x[0];
  double mass = 0.938;
  double beta = p[0];
  double T    = p[1];

  TF1 *fIntBG = 0;
  if(!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGP, 0., 1., 3);

  fIntBG->SetParameters(pT, beta, T);
  //fIntBG->SetNpx(100000);
  return fIntBG->Integral(0.,1.) * p[2]; //p[2] norm
}
void star_p_pi_bw(const int scan = 0)
{
  const Double_t R_PPi = 0.1; //rescale to P/Pi = 
  const Double_t T = 0.09;
  const Double_t betaMax = 0.9;

  TF1 *funBW = new TF1("funBW", StaticBGdNdPtPi, 0.2, 10., 3);
  funBW->SetParameters(betaMax, T, 1.);
  funBW->SetLineStyle(7);
  funBW->SetLineColor(4);
  funBW->SetLineWidth(4);

  TF1 *funBWP = new TF1("funBWP", StaticBGdNdPtP, 0.2, 10., 3);
  funBWP->SetParameters(betaMax, T, 1.);

  double norm_Pi = funBW->Integral(0,10.);
  double norm_P = funBWP->Integral(0,10.);
  funBWP->SetParameter(2, R_PPi*norm_Pi/norm_P);

  cout << " norm_Pi = " << norm_Pi << "\t norm_P = " << norm_P << "\t scale = " << R_PPi*norm_Pi/norm_P << endl;
  
  funBWP->SetLineStyle(1);
  funBWP->SetLineColor(1);
  funBWP->SetLineWidth(4);

  cout << " Check yields:\t Pi = " << funBW->Integral(0,10.) << "\t P = " << funBWP->Integral(0,10.) << endl;

  double pt[100],r[100];
  for(int i=0;i<100;i++) {
    pt[i] = (i+0.5)*0.1+0.2;
    r[i] = funBWP->Eval(pt[i])/funBW->Eval(pt[i]);
  }
  TGraph *gr = new TGraph(100,pt,r);
  gr->SetLineWidth(2);
  gr->SetLineColor(2);

  

  TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(0);
  c1->SetFrameBorderMode(0);
  c1->SetLeftMargin(0.18);
  c1->SetBottomMargin(0.15);
  c1->SetTopMargin(0.025);
  c1->SetRightMargin(0.025);
  c1->Draw();
  
  const float xmin   = 0.;
  const float xmax   = 5.;
  float ymin = 0;
  float ymax = 5.;

  TH1D *htmp = new TH1D("htmp","",1,xmin, xmax);
  htmp->GetYaxis()->SetTitle("Baryon/Meson");
  htmp->GetYaxis()->SetTitleSize(0.07);
  htmp->GetYaxis()->SetTitleOffset(1.3);
  htmp->GetYaxis()->SetLabelSize(0.055);
  htmp->SetMinimum(ymin);
  htmp->SetMaximum(ymax);
  htmp->GetXaxis()->SetNdivisions(208);
  //  htmp->GetXaxis()->CenterTitle();
  htmp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  htmp->GetXaxis()->SetTitleOffset(1.1);
  htmp->GetXaxis()->SetTitleSize(0.07);
  htmp->GetXaxis()->SetLabelOffset(0.01);
  htmp->GetXaxis()->SetLabelSize(0.055);
  htmp->GetXaxis()->SetLabelFont(42);
  htmp->GetXaxis()->SetTitleFont(42);
  htmp->GetYaxis()->SetNdivisions(208);
  htmp->Draw();
   
  gr->Draw("c");  

}
