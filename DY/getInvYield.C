void getInvYield(const char* ver="RHIC")
{
  gROOT->Reset();

  const Float_t Xsec[2][3] = { 3.735e-4, 1.355e-5, 2.876e-6,
			       4.319e-4, 1.322e-6, 3.586e-8};   // pb for qqbar->gamma* and to e+e-, mu+mu- only.  sqrt(s) = 0.2-2 GeV, 2-4, 4-10, RHIC, SPS

  const Int_t Nbins = 1000;
  const Float_t MassMax = 10.;
  const Float_t Nevts = 1.e6;
  const Int_t nRebin = 10;
  const Int_t CountCut = 1.5;  // bin count below cut is ignored to avoid large fluctuations.

  const Float_t ppXsec[2] = {42.0, 32.6};

  int index = -1;
  if(strcmp(ver,"RHIC")==0) {
    index = 0;
  } else if(strcmp(ver,"SPS")==0) {
    index = 1;
  } else {
    return;
  }


  Float_t Nevts_e;
  Float_t Nevts_mu;

  TH1F *invmTE1_slice[50];
  TH1F *invmTE2_slice[50];
  TH1F *invmTE3_slice[50];
  TH1F *invmTMu1_slice[50];
  TH1F *invmTMu2_slice[50];
  TH1F *invmTMu3_slice[50];

  TH1F *invmTEAll_slice[50];
  TH1F *invmTMuAll_slice[50];
  TH1F *invmTEAll = new TH1F("invmTEAll","",Nbins/nRebin,0.,MassMax);
  TH1F *invmTMuAll = new TH1F("invmTMuAll","",Nbins/nRebin,0.,MassMax);


  char inname[100];
  sprintf(inname,"root/test_%s_0_2.root",ver);
  TFile *f1 = new TFile(inname);
  TH2F *invmTE2D1 = (TH2F *)f1->Get("dNdmTEMid");
  TH2F *invmTMu2D1 = (TH2F *)f1->Get("dNdmTMuMid");
  TH1F *invmTE1 = (TH1F *)invmTE2D1->ProjectionX();
  TH1F *invmTMu1 = (TH1F *)invmTMu2D1->ProjectionX();

  Int_t nBins = Nbins/nRebin;
  Float_t dm = MassMax/(Nbins/nRebin);

  invmTE1->Rebin(nRebin);
  invmTMu1->Rebin(nRebin);

  sprintf(inname,"root/test_%s_2_4.root",ver);
  TFile *f2 = new TFile(inname);
  TH2F *invmTE2D2 = (TH2F *)f2->Get("dNdmTEMid");
  TH2F *invmTMu2D2 = (TH2F *)f2->Get("dNdmTMuMid");
  TH1F *invmTE2 = (TH1F *)invmTE2D2->ProjectionX();
  TH1F *invmTMu2 = (TH1F *)invmTMu2D2->ProjectionX();

  invmTE2->Rebin(nRebin);
  invmTMu2->Rebin(nRebin);

  sprintf(inname,"root/test_%s_4_10.root",ver);
  TFile *f3 = new TFile(inname);
  TH2F *invmTE2D3 = (TH2F *)f3->Get("dNdmTEMid");
  TH2F *invmTMu2D3 = (TH2F *)f3->Get("dNdmTMuMid");
  TH1F *invmTE3 = (TH1F *)invmTE2D3->ProjectionX();
  TH1F *invmTMu3 = (TH1F *)invmTMu2D3->ProjectionX();

  invmTE3->Rebin(nRebin);
  invmTMu3->Rebin(nRebin);

  for(int i=0;i<Nbins/nRebin;i++) {
    double a1 = invmTE1->GetBinContent(i+1);
    double a2 = invmTE2->GetBinContent(i+1);
    double a3 = invmTE3->GetBinContent(i+1);
    if(a1<CountCut) a1 = 0;
    if(a2<CountCut) a2 = 0;
    if(a3<CountCut) a3 = 0;
    invmTEAll->SetBinContent(i+1,(a1*Xsec[index][0]+a2*Xsec[index][1]+a3*Xsec[index][2])/dm/Nevts/ppXsec[index]);

    a1 = invmTMu1->GetBinContent(i+1);
    a2 = invmTMu2->GetBinContent(i+1);
    a3 = invmTMu3->GetBinContent(i+1);
    if(a1<CountCut) a1 = 0;
    if(a2<CountCut) a2 = 0;
    if(a3<CountCut) a3 = 0;
    invmTMuAll->SetBinContent(i+1,(a1*Xsec[index][0]+a2*Xsec[index][1]+a3*Xsec[index][2])/dm/Nevts/ppXsec[index]);
  }

  for(int i=0;i<50;i++) {
    double invM = (i+0.5)*0.2;

    char name[100];
    sprintf(name,"mTInvE1_%d",i);
    invmTE1_slice[i] = new TH1F(name,"",Nbins,0.,MassMax);
    invmTE1_slice[i] = (TH1F *)invmTE2D1->ProjectionX(name,i*2+1,i*2+2);
    invmTE1_slice[i]->SetName(name);
    invmTE1_slice[i]->Rebin(nRebin);

    sprintf(name,"mTInvE2_%d",i);
    invmTE2_slice[i] = new TH1F(name,"",Nbins,0.,MassMax);
    invmTE2_slice[i] = (TH1F *)invmTE2D2->ProjectionX(name,i*2+1,i*2+2);
    invmTE2_slice[i]->SetName(name);
    invmTE2_slice[i]->Rebin(nRebin);

    sprintf(name,"mTInvE3_%d",i);
    invmTE3_slice[i] = new TH1F(name,"",Nbins,0.,MassMax);
    invmTE3_slice[i] = (TH1F *)invmTE2D3->ProjectionX(name,i*2+1,i*2+2);
    invmTE3_slice[i]->SetName(name);
    invmTE3_slice[i]->Rebin(nRebin);

    sprintf(name,"mTInvMu1_%d",i);
    invmTMu1_slice[i] = new TH1F(name,"",Nbins,0.,MassMax);
    invmTMu1_slice[i] = (TH1F *)invmTMu2D1->ProjectionX(name,i*2+1,i*2+2);
    invmTMu1_slice[i]->SetName(name);
    invmTMu1_slice[i]->Rebin(nRebin);

    sprintf(name,"mTInvMu2_%d",i);
    invmTMu2_slice[i] = new TH1F(name,"",Nbins,0.,MassMax);
    invmTMu2_slice[i] = (TH1F *)invmTMu2D2->ProjectionX(name,i*2+1,i*2+2);
    invmTMu2_slice[i]->SetName(name);
    invmTMu2_slice[i]->Rebin(nRebin);
    
    sprintf(name,"mTInvMu3_%d",i);
    invmTMu3_slice[i] = new TH1F(name,"",Nbins,0.,MassMax);
    invmTMu3_slice[i] = (TH1F *)invmTMu2D3->ProjectionX(name,i*2+1,i*2+2);
    invmTMu3_slice[i]->SetName(name);
    invmTMu3_slice[i]->Rebin(nRebin);

    sprintf(name,"mTInvEAll_%d",i);
    invmTEAll_slice[i]= new TH1F(name,"",Nbins/nRebin,0.,MassMax);
    for(int j=0;j<Nbins/nRebin;j++) {
      double mTm0 = invmTEAll_slice[i]->GetBinCenter(j+1);
      double mT = mTm0 + invM;

      double ra1 = invmTE1_slice[i]->GetBinContent(j+1)*mT;
      double ra2 = invmTE2_slice[i]->GetBinContent(j+1)*mT;
      double ra3 = invmTE3_slice[i]->GetBinContent(j+1)*mT;
      if(ra1<CountCut||i>10) ra1 = 0.;
      if(ra2<CountCut||i<10||i>20) ra2 = 0.;
      if(ra3<CountCut||i<20) ra3 = 0.;
     
      double ra = ra1*Xsec[index][0] + ra2*Xsec[index][1] + ra3*Xsec[index][2];
      double re = sqrt(ra1*pow(Xsec[index][0],2.)+ra2*pow(Xsec[index][1],2.)+ra3*pow(Xsec[index][2],2.));
      ra *= 1./mT/dm/Nevts/ppXsec[index];
      re *= 1./mT/dm/Nevts/ppXsec[index];
      invmTEAll_slice[i]->SetBinContent(j+1,ra);
      invmTEAll_slice[i]->SetBinError(j+1,re);
    }

    sprintf(name,"mTInvMuAll_%d",i);
    invmTMuAll_slice[i]= new TH1F(name,"",Nbins/nRebin,0.,MassMax);
    for(int j=0;j<Nbins/nRebin;j++) {
      double mTm0 = invmTMuAll_slice[i]->GetBinCenter(j+1);
      double mT = mTm0 + invM;

      double ra1 = invmTMu1_slice[i]->GetBinContent(j+1)*mT;
      double ra2 = invmTMu2_slice[i]->GetBinContent(j+1)*mT;
      double ra3 = invmTMu3_slice[i]->GetBinContent(j+1)*mT;
      if(ra1<CountCut||i>10) ra1 = 0.;
      if(ra2<CountCut||i<10||i>20) ra2 = 0.;
      if(ra3<CountCut||i<20) ra3 = 0.;
     
      double ra = ra1*Xsec[index][0] + ra2*Xsec[index][1] + ra3*Xsec[index][2];
      double re = sqrt(ra1*pow(Xsec[index][0],2.)+ra2*pow(Xsec[index][1],2.)+ra3*pow(Xsec[index][2],2.));
      ra *= 1./mT/dm/Nevts/ppXsec[index];
      re *= 1./mT/dm/Nevts/ppXsec[index];
      invmTMuAll_slice[i]->SetBinContent(j+1,ra);
      invmTMuAll_slice[i]->SetBinError(j+1,re);
    }


  }


  invmTE1->Scale(Xsec[index][0]/dm/Nevts/ppXsec[index]);
  invmTMu1->Scale(Xsec[index][0]/dm/Nevts/ppXsec[index]);
  invmTE2->Scale(Xsec[index][1]/dm/Nevts/ppXsec[index]);
  invmTMu2->Scale(Xsec[index][1]/dm/Nevts/ppXsec[index]);
  invmTE3->Scale(Xsec[index][2]/dm/Nevts/ppXsec[index]);
  invmTMu3->Scale(Xsec[index][2]/dm/Nevts/ppXsec[index]);

  for(int i=0;i<50;i++) {
    invmTE1_slice[i]->Scale(Xsec[index][0]/dm/Nevts/ppXsec[index]);
    invmTMu1_slice[i]->Scale(Xsec[index][0]/dm/Nevts/ppXsec[index]);
    invmTE2_slice[i]->Scale(Xsec[index][1]/dm/Nevts/ppXsec[index]);
    invmTMu2_slice[i]->Scale(Xsec[index][1]/dm/Nevts/ppXsec[index]);
    invmTE3_slice[i]->Scale(Xsec[index][2]/dm/Nevts/ppXsec[index]);
    invmTMu3_slice[i]->Scale(Xsec[index][2]/dm/Nevts/ppXsec[index]);
  }


  char outroot[100];
  sprintf(outroot,"mTInv_slice_%s.root",ver);
  TFile *fout = new TFile(outroot,"RECREATE");
  invmTE1->Write();
  invmTE2->Write();
  invmTE3->Write();
  invmTMu1->Write();
  invmTMu2->Write();
  invmTMu3->Write();
  invmTEAll->Write();
  invmTMuAll->Write();
  for(int i=0;i<50;i++) {
    invmTE1_slice[i]->Write();
    invmTE2_slice[i]->Write();
    invmTE3_slice[i]->Write();
    invmTEAll_slice[i]->Write();
    invmTMu1_slice[i]->Write();
    invmTMu2_slice[i]->Write();
    invmTMu3_slice[i]->Write();
    invmTMuAll_slice[i]->Write();
  }
  fout->Close();


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
   double y1 = 5.e-18;
   double y2 = 1.e-5;

   TH1D *h0 = new TH1D("h0","",1,x1,x2);
   h0->SetMaximum(y2);
   h0->SetMinimum(y1);
   h0->GetXaxis()->SetNdivisions(208);
   h0->GetXaxis()->CenterTitle();
   h0->GetXaxis()->SetTitle("m_{T} - m_{inv} (GeV/c^{2})");
   h0->GetXaxis()->SetTitleOffset(0.9);
   h0->GetXaxis()->SetTitleSize(0.07);
   h0->GetXaxis()->SetLabelOffset(0.01);
   h0->GetXaxis()->SetLabelSize(0.05);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetYaxis()->SetNdivisions(204);
   h0->GetYaxis()->SetTitle("d^{2}N/m_{T}dm_{T}dy ((GeV/c^{2})^{-2})");
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

   invmTE1->SetLineWidth(2);
   invmTE1->SetLineColor(2);
   invmTE1->Draw("same h");

   invmTMu1->SetLineWidth(2);
   invmTMu1->SetLineStyle(2);
   invmTMu1->SetLineColor(2);
   invmTMu1->Draw("same h");

   invmTE2->SetLineWidth(2);
   invmTE2->SetLineColor(4);
   invmTE2->Draw("same h");

   invmTMu2->SetLineWidth(2);
   invmTMu2->SetLineStyle(2);
   invmTMu2->SetLineColor(4);
   invmTMu2->Draw("same h");

   invmTE3->SetLineWidth(2);
   invmTE3->SetLineColor(6);
   invmTE3->Draw("same h");

   invmTMu3->SetLineWidth(2);
   invmTMu3->SetLineStyle(2);
   invmTMu3->SetLineColor(6);
   invmTMu3->Draw("same h");

   invmTEAll->SetLineWidth(2);
   invmTEAll->SetLineColor(1);
   invmTEAll->Draw("same h");

   invmTMuAll->SetLineWidth(2);
   invmTMuAll->SetLineStyle(2);
   invmTMuAll->SetLineColor(1);
   invmTMuAll->Draw("same h");


   TLegend *leg = new TLegend(0.70, 0.64, 0.96, 0.94);
   leg->SetFillColor(10);
   leg->SetFillStyle(0);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetTextFont(42);
   leg->SetTextSize(0.05);
//    leg->AddEntry(invMassE," e^{+}e^{-}","l");
//    leg->AddEntry(invMassMu," #mu^{+}#mu^{-}","l");
//    leg->Draw();

   c1->Update();
   char outname[100];
   sprintf(outname,"fig/mTInv_%s.eps",ver);
   c1->SaveAs(outname);
   sprintf(outname,"fig/mTInv_%s.gif",ver);
   c1->SaveAs(outname);


}
