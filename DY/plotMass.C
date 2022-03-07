void plotMass(const char *ver = "RHIC")
{
  gROOT->Reset();

  const Float_t Xsec[2][3] = { 3.735e-4, 1.355e-5, 2.876e-6,
			       4.319e-4, 1.322e-6, 3.586e-8};   // pb for qqbar->gamma* and to e+e-, mu+mu- only.  sqrt(s) = 0.2-2 GeV, 2-4, 4-10, RHIC, SPS

  const Int_t Nbins = 1000;
  const Float_t MassMax = 10.;
  const Float_t Nevts = 1.e6;

  const Int_t nRebin = 10;
  const Float_t ppXsec[2] = {42.0, 32.6};

  Int_t nBins = Nbins/nRebin;
  Float_t dm = MassMax/(Nbins/nRebin);

  int index = -1;
  if(strcmp(ver,"RHIC")==0) {
    index = 0;
  } else if(strcmp(ver,"SPS")==0) {
    index = 1;
  } else {
    return;
  }

  char inname[100];
  sprintf(inname,"root/test_%s_0_2.root",ver);
  TFile *f1 = new TFile(inname);
  TH1F *invMassE1 = (TH1F *)f1->Get("invMassEMid");
  TH1F *invMassMu1 = (TH1F *)f1->Get("invMassMuMid");  
  invMassE1->Rebin(nRebin);
  invMassMu1->Rebin(nRebin);
  invMassE1->Sumw2();
  invMassMu1->Sumw2();
  invMassE1->Scale(Xsec[index][0]/dm/Nevts/ppXsec[index]);
  invMassMu1->Scale(Xsec[index][0]/dm/Nevts/ppXsec[index]);

  sprintf(inname,"root/test_%s_2_4.root",ver);
  TFile *f2 = new TFile(inname);
  TH1F *invMassE2 = (TH1F *)f2->Get("invMassEMid");
  TH1F *invMassMu2 = (TH1F *)f2->Get("invMassMuMid");  
  invMassE2->Rebin(nRebin);
  invMassMu2->Rebin(nRebin);
  invMassE2->Sumw2();
  invMassMu2->Sumw2();
  invMassE2->Scale(Xsec[index][1]/dm/Nevts/ppXsec[index]);
  invMassMu2->Scale(Xsec[index][1]/dm/Nevts/ppXsec[index]);

  sprintf(inname,"root/test_%s_4_10.root",ver);
  TFile *f3 = new TFile(inname);
  TH1F *invMassE3 = (TH1F *)f3->Get("invMassEMid");
  TH1F *invMassMu3 = (TH1F *)f3->Get("invMassMuMid");  
  invMassE3->Rebin(nRebin);
  invMassMu3->Rebin(nRebin);
  invMassE3->Sumw2();
  invMassMu3->Sumw2();
  invMassE3->Scale(Xsec[index][2]/dm/Nevts/ppXsec[index]);
  invMassMu3->Scale(Xsec[index][2]/dm/Nevts/ppXsec[index]);


  TH1F invMassEAll = (*invMassE1) + (*invMassE2) + (*invMassE3);
  TH1F invMassMuAll = (*invMassMu1) + (*invMassMu2) + (*invMassMu3);


   TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,600);
   gStyle->SetOptFit(0);
   gStyle->SetOptStat(0);
   gStyle->SetEndErrorSize(0.01);
   c1->SetFillColor(10);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameFillColor(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SetLogy();
   c1->SetGridx();
   c1->SetGridy();
   c1->SetLeftMargin(0.14);
   c1->SetBottomMargin(0.15);
   c1->SetTopMargin(0.02);
   c1->SetRightMargin(0.02);

   double x1 = 0.;
   double x2 = 10.;
   double y1 = 5.e-12;
   double y2 = 20.e-8;

   TH1D *h0 = new TH1D("h0","",1,x1,x2);
   h0->SetMaximum(y2);
   h0->SetMinimum(y1);
   h0->GetXaxis()->SetNdivisions(208);
   h0->GetXaxis()->CenterTitle();
   h0->GetXaxis()->SetTitle("M_{l^{+}l^{-}} (GeV/c^{2})");
   h0->GetXaxis()->SetTitleOffset(0.9);
   h0->GetXaxis()->SetTitleSize(0.07);
   h0->GetXaxis()->SetLabelOffset(0.01);
   h0->GetXaxis()->SetLabelSize(0.05);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetYaxis()->SetNdivisions(204);
   h0->GetYaxis()->SetTitle("dN/dM (c^{2}/GeV)");
   h0->GetYaxis()->SetTitleOffset(1.0);
   h0->GetYaxis()->SetTitleSize(0.07);
   h0->GetYaxis()->SetLabelOffset(0.01);
   h0->GetYaxis()->SetLabelSize(0.05);
   h0->GetYaxis()->SetLabelFont(42);
   h0->GetYaxis()->SetTitleFont(42);
   h0->Draw();


   TLine *l1 = new TLine(x1,y1,x2,y1);
   l1->SetLineWidth(3);
   l1->Draw("same");
   TLine *l2 = new TLine(x1,y2,x2,y2);
   l2->SetLineWidth(3);
   l2->Draw("same");
   TLine *l3 = new TLine(x1,y1,x1,y2);
   l3->SetLineWidth(3);
   l3->Draw("same");
   TLine *l4 = new TLine(x2,y1,x2,y2);
   l4->SetLineWidth(3);
   l4->Draw("same");

   invMassEAll.SetLineWidth(2);
   invMassEAll.SetLineColor(2);
   invMassEAll.Draw("same");

   invMassE1->SetLineColor(4);
   invMassE1->SetLineWidth(2);
   invMassE1->Draw("same");
   invMassE2->SetLineColor(4);
   invMassE2->SetLineWidth(2);
   invMassE2->Draw("same");
   invMassE3->SetLineColor(4);
   invMassE3->SetLineWidth(2);
   invMassE3->Draw("same");


//    invMassMuAll.SetLineWidth(2);
//    invMassMuAll.SetLineColor(4);
//    invMassMuAll.Draw("same");

   TLegend *leg = new TLegend(0.70, 0.64, 0.96, 0.94);
   leg->SetFillColor(10);
   leg->SetFillStyle(0);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetTextFont(42);
   leg->SetTextSize(0.05);
//    leg->AddEntry(&invMassE," e^{+}e^{-}","l");
//    leg->AddEntry(&invMassMu," #mu^{+}#mu^{-}","l");
//    leg->Draw();

   c1->Update();
   char outname[100];
   sprintf(outname,"fig/DY_Mass_%s.eps",ver);
   c1->SaveAs(outname);
   sprintf(outname,"fig/DY_Mass_%s.gif",ver);
   c1->SaveAs(outname);


}
