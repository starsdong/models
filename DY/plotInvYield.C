Double_t mTExpo(Double_t *x, Double_t *par)
{
  const Double_t Pi = 3.1415927;
  Double_t dNdy = par[0];
  Double_t T = par[1];
  Double_t M0 = par[2];
  Double_t MT = x[0];
  return dNdy/2./Pi/T/(M0+T)*exp(-(MT-M0)/T);
}

void plotInvYield(const char *ver = "RHIC")
{
  gROOT->Reset();

  Double_t mT[100];
  Double_t yy1[50][100], yy1e[50][100];
  Double_t yy2[50][100], yy2e[50][100];

  char title[50][100];

  char inname[100];
  sprintf(inname,"mTInv_slice_%s.root",ver);
  Int_t index = -1;
  if(strcmp(ver,"RHIC")==0) index = 0;
  else if(strcmp(ver,"SPS")==0) index = 1;
  else return;

  TFile *fin = new TFile(inname);
  for(int i=0;i<50;i++) {
    char hisname[100];
    sprintf(hisname,"mTInvEAll_%d",i);
    TH1F *htmp1 = (TH1F *)fin->Get(hisname);
    for(int j=0;j<100;j++) {
      if(i==0) mT[j] = htmp1->GetBinCenter(j+1);
      yy1[i][j] = htmp1->GetBinContent(j+1);
      yy1e[i][j] = htmp1->GetBinError(j+1);
    }

    sprintf(hisname,"mTInvMuAll_%d",i);
    TH1F *htmp2 = (TH1F *)fin->Get(hisname);
    for(int j=0;j<100;j++) {
      yy2[i][j] = htmp2->GetBinContent(j+1);
      yy2e[i][j] = htmp2->GetBinError(j+1);
    }
    sprintf(title[i],"%3.1f<M_{inv}<%3.1f",i*0.2,(i+1)*0.2);
  }

  TGraphErrors *gr1[50];
  TGraphErrors *gr2[50];

  for(int i=0;i<50;i++) {
    gr1[i] = new TGraphErrors(100,mT,yy1[i],0,yy1e[i]);
    gr2[i] = new TGraphErrors(100,mT,yy2[i],0,yy2e[i]);
  }


  TCanvas *c1 = new TCanvas("c1", "c1",0,0,1000,1000);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0.01);
  c1->SetFillColor(10);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetFrameFillColor(10);
  c1->SetFrameBorderMode(0);

  TPad *pad[25];
  double XMIN = 0.13;
  double XMAX = 0.98;
  double YMIN = 0.12;
  double YMAX = 0.98;

  double yMax[5] = {1.e-5, 1.e-6, 1.e-7, 3.e-8, 1.e-8};

  TF1 *fun1 = new TF1("fun1",mTExpo,0.0,5.0,3);

  for(int i=0;i<25;i++) {
    c1->cd();

    int irow = i/5;
    int icol = i%5;
    double xp1 = XMIN + (XMAX-XMIN)/5*icol;
    double xp2 = XMIN + (XMAX-XMIN)/5*(icol+1);
    double yp1 = YMAX - (YMAX-YMIN)/5*(irow+1);
    double yp2 = YMAX - (YMAX-YMIN)/5*irow;
    
    char pname[100];
    sprintf(pname,"pad_pT_%d_pie_%d", irow, icol);
    pad[i] = new TPad(pname,"",xp1,yp1,xp2,yp2);
    pad[i]->SetFillColor(10);
    pad[i]->SetBorderMode(0);
    pad[i]->SetBorderSize(2);
    pad[i]->SetFillStyle(0);
    pad[i]->SetFrameFillColor(10);
    pad[i]->SetFrameFillStyle(0);
    pad[i]->SetFrameBorderMode(0);
    pad[i]->SetTopMargin(0.001);
    pad[i]->SetBottomMargin(0.001);
    pad[i]->SetLeftMargin(0.001);
    pad[i]->SetRightMargin(0.001);
    pad[i]->SetLogy();
    pad[i]->Draw();
    pad[i]->cd();

    double x1 = 0.;
    double x2 = 3.99;
    double y1 = yMax[irow]*1.e-6;
    if(index==1) y1 = y1*1.e-2;
    double y2 = yMax[irow];

    TH1D *h0 = new TH1D("h0","",1,x1,x2);
    h0->SetMinimum(y1);
    h0->SetMaximum(y2);
    h0->GetXaxis()->SetNdivisions(205);
    h0->GetXaxis()->SetTitle("");
    h0->GetXaxis()->SetTitleOffset(999.);
    h0->GetXaxis()->SetTitleSize(0.0001);
    h0->GetXaxis()->SetLabelOffset(999.);
    h0->GetXaxis()->SetLabelSize(0.0001);
    h0->GetXaxis()->SetLabelFont(42);
    h0->GetXaxis()->SetTitleFont(42);
    h0->GetYaxis()->SetNdivisions(105);
    h0->GetYaxis()->SetTitle("");
    h0->GetYaxis()->SetTitleOffset(999.);
    h0->GetYaxis()->SetTitleSize(0.0001);
    h0->GetYaxis()->SetLabelOffset(999.);
    h0->GetYaxis()->SetLabelSize(0.0001);
    h0->GetYaxis()->SetLabelFont(42);
    h0->GetYaxis()->SetTitleFont(42);
    h0->Draw("c");

    gr1[i]->SetMarkerStyle(20);
    gr1[i]->SetMarkerSize(1.5);
    gr1[i]->SetLineWidth(2);

    gr2[i]->SetMarkerStyle(24);
    gr2[i]->SetMarkerSize(1.5);
    gr2[i]->SetLineWidth(2);

    // expo fit
    /*
    if(i==6) {
      gr1[i]->Draw("p");
      fun1->SetParameters(1e-5,0.4,1.3);
      fun1->FixParameter(2,1.3);
      fun1->SetRange(0.0, 3.0);
      gr1[i]->Fit("fun1","R");
    } else {
      continue;
    }
    */
    
    
    double Tslope[30], Terr[30];
    if(i!=0) {
      gr1[i]->Draw("p");      
      gr2[i]->Draw("p");

      fun1->SetParameters(1e-8,0.2,(i+0.5)*0.2);
      fun1->FixParameter(2,(i+0.5)*0.2);
      fun1->SetRange(0.1, 1.);
      gr1[i]->Fit("fun1","R");

      Tslope[i] = fun1->GetParameter(1);
      Terr[i] = fun1->GetParError(1);
      
      TLatex *tex = new TLatex(1.5, yMax[irow]*1.e-2, title[i]);
      tex->SetTextFont(42);
      tex->SetTextSize(0.12);
      tex->Draw("same");
      
    } else {
      char text[100];
      sprintf(text,"p+p @ %s",ver);
      TLatex *tex = new TLatex(0.5, yMax[irow]*1e-2, text);
      tex->SetTextFont(42);
      tex->SetTextSize(0.15);
      tex->Draw("same");

      TLegend *leg = new TLegend(0.30, 0.14, 0.7, 0.54);
      leg->SetFillColor(10);
      leg->SetFillStyle(0);
      leg->SetLineStyle(4000);
      leg->SetLineColor(10);
      leg->SetLineWidth(0.);
      leg->SetTextFont(42);
      leg->SetTextSize(0.15);
      leg->AddEntry(gr1[i]," e^{+}e^{-}","p");
      leg->AddEntry(gr2[i]," #mu^{+}#mu^{-}","p");
      leg->Draw();

    }


    TLine *l1 = new TLine(x1,y1,x2,y1);
    l1->SetLineWidth(3);
    if(irow==4) l1->Draw("same");
    TLine *l2 = new TLine(x1,y2,x2,y2);
    l2->SetLineWidth(3);
    if(irow==0) l2->Draw("same");
    TLine *l3 = new TLine(x1,y1,x1,y2);
    l3->SetLineWidth(3);
    if(icol==0) l3->Draw("same");
    TLine *l4 = new TLine(x2,y1,x2,y2);
    l4->SetLineWidth(3);
    if(icol==4) l4->Draw("same");
    
    pad[i]->Modified();
  }

  c1->cd();

  ////////////////////////////////
  // put on the labels on axis
  ////////////////////////////////
  for(int i=0;i<5;i++) {
    double xp1 = XMIN + (XMAX-XMIN)/5*i;
    double xp2 = XMIN + (XMAX-XMIN)/5*(i+1);

    double x0 = xp1 + ( 0 - 0 )/4.*(xp2-xp1) - 0.005;
    double x1 = xp1 + ( 1 - 0 )/4.*(xp2-xp1) - 0.005;
    double x2 = xp1 + ( 2 - 0 )/4.*(xp2-xp1) - 0.005;
    double x3 = xp1 + ( 3 - 0 )/4.*(xp2-xp1) - 0.005;

    TLatex *tex = new TLatex(x0, YMIN - 0.022, "0");
    tex->SetTextSize(0.025);
    tex->SetTextFont(42);
    tex->Draw("same");

    TLatex *tex = new TLatex(x1, YMIN - 0.022, "1");
    tex->SetTextSize(0.025);
    tex->SetTextFont(42);
    tex->Draw("same");

    TLatex *tex = new TLatex(x2, YMIN - 0.022, "2");
    tex->SetTextSize(0.025);
    tex->SetTextFont(42);
    tex->Draw("same");

    TLatex *tex = new TLatex(x3, YMIN - 0.022, "3");
    tex->SetTextSize(0.025);
    tex->SetTextFont(42);
    tex->Draw("same");
  }

  for(int i=0;i<5;i++) {
    double yp1 = YMAX - (YMAX-YMIN)/5*(i+1);
    double yp2 = YMAX - (YMAX-YMIN)/5*i;

    double ymax = yMax[i];
    double ymin = yMax[i]*1.e-6;
    if(index==1) ymin *= 1.e-2;

    int yL1 = (int)(log10(ymin))+1;
    int yL2 = (int)(log10(ymax))-1;

    double yLp1 = yp1 + (yL1 - log10(ymin)) / (log10(ymax)-log10(ymin)) * (yp2 - yp1) - 0.005;
    double yLp2 = yp1 + (yL2 - log10(ymin)) / (log10(ymax)-log10(ymin)) * (yp2 - yp1) - 0.005;


    char label[10];
    sprintf(label,"10^{%d}",yL1);
    TLatex *tex = new TLatex(XMIN-0.05, yLp1, label);
    tex->SetTextSize(0.025);
    tex->SetTextFont(42);
    tex->Draw("same");
    
    sprintf(label,"10^{%d}",yL2);
    TLatex *tex = new TLatex(XMIN-0.05, yLp2, label);
    tex->SetTextSize(0.025);
    tex->SetTextFont(42);
    tex->Draw("same");
    
  }
  
  TLatex *tex = new TLatex(XMIN+0.2, YMIN-0.09, "m_{T} - m_{inv} (GeV/c^{2})");
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->Draw("same");

  TLatex *tex = new TLatex(XMIN-0.08, YMAX/2, "d^{2}N/m_{T}dm_{T}dy (c^{4}/GeV^{2})");
  tex->SetTextFont(42);
  tex->SetTextSize(0.05);
  tex->SetTextAngle(90);
  tex->Draw("same");

  c1->Update();
  char outfig[100];
  sprintf(outfig,"fig/invYield_%s.eps",ver);
  c1->SaveAs(outfig);
  sprintf(outfig,"fig/invYield_%s.gif",ver);
  c1->SaveAs(outfig);

  for(int i=0;i<25;i++) {
    cout << " T = " << Tslope[i] << " +/- " << Terr[i] << endl;
  }
  
  
}
