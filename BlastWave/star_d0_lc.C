Double_t IntegrandBGD0(const double *x, const double *p)
{

  double x0   = x[0];
  double mass = p[3];
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

Double_t StaticBGdNdPtD0(const double *x, const double *p)
{
  double pT   = x[0];
  double mass = p[2];
  double beta = p[0];
  double T    = p[1];

  TF1 *fIntBG = 0;
  if(!fIntBG) fIntBG = new TF1("fIntBG", IntegrandBGD0, 0., 1., 4);

  fIntBG->SetParameters(pT, beta, T, mass);
  //fIntBG->SetNpx(100000);
  return fIntBG->Integral(0.,1.) * p[3]; //p[3] norm
}
void star_d0_lc(const int scan = 0)
{
  const Double_t massD = 1.865;
  const Double_t massLc = 2.286;
  TF1 *funBW = new TF1("funBW", StaticBGdNdPtD0, 0.2, 5., 3);
  funBW->FixParameter(0,0.488197);
  funBW->FixParameter(1,0.215270);
  funBW->FixParameter(2,massD);
  funBW->FixParameter(3,3606.473492);
  funBW->SetLineStyle(7);
  funBW->SetLineColor(4);
  funBW->SetLineWidth(4);

  TCanvas *cc = new TCanvas("cc", "cc",0,0,800,600);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  cc->SetFillColor(10);
  cc->SetBorderMode(0);
  cc->SetBorderSize(2);
  cc->SetFrameFillColor(0);
  cc->SetFrameBorderMode(0);
  cc->SetLogy();
  //  cc->SetGridx();
  //  cc->SetGridy();
  cc->SetLeftMargin(0.18);
  cc->SetBottomMargin(0.15);
  cc->SetTopMargin(0.025);
  cc->SetRightMargin(0.025);
  
  const float xmin   = 0.;
  const float xmax   = 5.1;
  float ymin = 1e-6;
  float ymax = 1.6;

  TH1D *htmp = new TH1D("htmp","",1,xmin, xmax);
  htmp->GetYaxis()->SetTitle("d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
  //htmp->GetYaxis()->CenterTitle();
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
   
   TGraphErrors *gr0 = new TGraphErrors(6);
   gr0->SetPoint(0,0.314591,0.1221096);
   gr0->SetPointError(0,0,0.0264988);
   gr0->SetPoint(1,0.967794,0.1085676);
   gr0->SetPointError(1,0,0.0206796);
   gr0->SetPoint(2,1.35482,0.0482164);
   gr0->SetPointError(2,0,0.00833648);
   gr0->SetPoint(3,1.88638,0.0203288);
   gr0->SetPointError(3,0,0.00340368);
   gr0->SetPoint(4,2.56187,0.0041562);
   gr0->SetPointError(4,0,0.000910948);
   gr0->SetPoint(5,3.74846,0.000330338);
   gr0->SetPointError(5,0,0.0000869312);
   gr0->SetMarkerColor(1);
   gr0->SetMarkerSize(1.0);
   gr0->SetLineColor(1);
   gr0->SetMarkerStyle(20);
   gr0->Draw("P");

  
  funBW->Draw("same");
  
  const float txx = 0.5;


  TLatex tx;
  tx.SetNDC();
  tx.DrawLatex(txx,0.9,"(D^{0}+#bar{D^{0}})/2, |y| < 1");
  tx.DrawLatex(txx,0.84,"10-40% Au+Au 200 GeV");
  char txbeta[100];
  sprintf(txbeta,"<#beta^{1}_{T}> = %4.2f #pm %4.2f",funBW->GetParameter(0)*2./3.,0.11);
  tx.DrawLatex(txx,0.8,txbeta);
  char txtemp[100];
  sprintf(txtemp,"T^{1} = %d #pm %d MeV",funBW->GetParameter(1)*1000,89.);
  tx.DrawLatex(txx,0.75,txtemp);
  
  TLegend *lg = new TLegend(0.23,0.28,0.62,0.48);
  lg->SetFillStyle(0);
  lg->SetFillColor(10);
  lg->SetBorderSize(0);
  //lg->SetTextFont(22);
  lg->SetTextSize(0.055);
  lg->AddEntry(funBW,"Blast-wave fit","l");
  lg->Draw();
  
}
