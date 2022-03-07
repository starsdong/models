void plotMeanPt(const char *ver = "RHIC")
{
  gROOT->Reset();

  const Float_t Xsec[2][3] = { 3.735e-4, 1.355e-5, 2.876e-6,
			       4.319e-4, 1.322e-6, 3.586e-8};   // pb for qqbar->gamma* and to e+e-, mu+mu- only.  sqrt(s) = 0.2-2 GeV, 2-4, 4-10, RHIC, SPS

  char inname[100];
  sprintf(inname,"root/test_%s_0_2.root",ver);
  TFile *f1 = new TFile(inname);
  TH2F *pTMassE1 = (TH2F *)f1->Get("dNdpTEMid");
  TH2F *pTMassMu1 = (TH2F *)f1->Get("dNdpTMuMid");

  sprintf(inname,"root/test_%s_2_4.root",ver);
  TFile *f2 = new TFile(inname);
  TH2F *pTMassE2 = (TH2F *)f2->Get("dNdpTEMid");
  TH2F *pTMassMu2 = (TH2F *)f2->Get("dNdpTMuMid");

  sprintf(inname,"root/test_%s_4_10.root",ver);
  TFile *f3 = new TFile(inname);
  TH2F *pTMassE3 = (TH2F *)f3->Get("dNdpTEMid");
  TH2F *pTMassMu3 = (TH2F *)f3->Get("dNdpTMuMid");

  int index = -1;
  if(strcmp(ver,"RHIC")==0) {
    index = 0;
  }  else if(strcmp(ver,"SPS")==0) {
    index = 1;
  }

  pTMassE1->Scale(Xsec[index][0]/Xsec[index][1]);
  pTMassMu1->Scale(Xsec[index][0]/Xsec[index][1]);

  pTMassE2->Scale(Xsec[index][1]/Xsec[index][1]);
  pTMassMu2->Scale(Xsec[index][1]/Xsec[index][1]);

  pTMassE3->Scale(Xsec[index][2]/Xsec[index][1]);
  pTMassMu3->Scale(Xsec[index][2]/Xsec[index][1]);

  TH2F pTMassE = (*pTMassE1) + (*pTMassE2) + (*pTMassE3);
  pTMassE.SetName("pTMassE");
  pTMassE.RebinX(5);
  pTMassE.ProfileY();

  TH2F pTMassMu = (*pTMassMu1) + (*pTMassMu2) + (*pTMassMu3);
  pTMassMu.SetName("pTMassMu");
  pTMassMu.RebinX(5);
  pTMassMu.ProfileY();

  
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
   c1->SetGridx();
   c1->SetGridy();
   c1->SetLeftMargin(0.12);
   c1->SetBottomMargin(0.15);
   c1->SetTopMargin(0.02);
   c1->SetRightMargin(0.02);

   double x1 = 0.2;
   double x2 = 9.8;
   double y1 = 0.;
   double y2 = 3.;

   TH1D *h0 = new TH1D("h0","",1,x1,x2);
   h0->SetMaximum(y2);
   h0->SetMinimum(y1);
   h0->GetXaxis()->SetNdivisions(208);
   h0->GetXaxis()->CenterTitle();
   h0->GetXaxis()->SetTitle("M_{inv} (GeV/c^{2})");
   h0->GetXaxis()->SetTitleOffset(0.9);
   h0->GetXaxis()->SetTitleSize(0.07);
   h0->GetXaxis()->SetLabelOffset(0.01);
   h0->GetXaxis()->SetLabelSize(0.05);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetYaxis()->SetNdivisions(204);
   h0->GetYaxis()->SetTitle("<p_{T}> (GeV/c)");
   h0->GetYaxis()->SetTitleOffset(0.8);
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

   pTMassE_pfy->SetFillColor(2);
   pTMassE_pfy->SetLineColor(2);
   pTMassE_pfy->SetLineWidth(10);
   pTMassE_pfy->Draw("e3 same");

   pTMassMu_pfy->SetFillColor(8);
   pTMassMu_pfy->SetLineColor(8);
   pTMassMu_pfy->SetLineWidth(10);
   pTMassMu_pfy->Draw("e3 same");

   char text[100];
   double yt;
   if(strcmp(ver,"RHIC")==0) {
     sprintf(text,"Drell-Yan @ #sqrt{s} = 200 GeV");
     yt = 0.3;
   } else if (strcmp(ver,"SPS")==0) {
     sprintf(text,"Drell-Yan @ #sqrt{s} = 17.2 GeV");
     yt = 2.3;
   }
   TLatex *tex = new TLatex(4.1, yt, text);
   tex->SetTextFont(42);
   tex->SetTextSize(0.05);
   tex->Draw("same");

   TLegend *leg = new TLegend(0.20, 0.64, 0.4, 0.94);
   leg->SetFillColor(10);
   leg->SetFillStyle(0);
   leg->SetLineStyle(4000);
   leg->SetLineColor(10);
   leg->SetLineWidth(0.);
   leg->SetTextFont(42);
   leg->SetTextSize(0.05);
   leg->AddEntry(pTMassE_pfy," e^{+}e^{-}","l");
   leg->AddEntry(pTMassMu_pfy," #mu^{+}#mu^{-}","l");
   leg->Draw();

   c1->Update();
   char outname[100];
   sprintf(outname,"fig/meanPt_%s.eps",ver);
   c1->SaveAs(outname);
   sprintf(outname,"fig/meanPt_%s.gif",ver);
   c1->SaveAs(outname);

}
