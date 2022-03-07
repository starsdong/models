#if (!defined(__CINT__) || defined(__MAKECINT__))
//#ifndef __CINT__
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TPythia6/TPythia6.h"
#include <iostream>
using std::cout;
using std::endl;
#endif

//---------------------------------------------------------------------
// Full function -- any specified pT, y spectrum
//     different D hadron fractions are fixed in the code.
//---------------------------------------------------------------------
void DrellYan(const Int_t Nevts=1000, const Char_t *OutName="test",const Char_t *Energy="RHIC",const Double_t CKIN1=0.0, const Double_t CKIN2=2.0)
{
  gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("$OPTSTAR/alt/lib/libPythia6");
  gSystem->Load("libEGPythia6");

  /////////////////////////////
  // Book the histogram
  ////////////////////////////
  TH1F *hInvMassEMid = new TH1F("invMassEMid","",1000,0.,10.);
  TH1F *hInvMassMuMid = new TH1F("invMassMuMid","",1000,0.,10.);
  TH1F *hInvMassE = new TH1F("invMassE","",1000,0.,10.);
  TH1F *hInvMassMu = new TH1F("invMassMu","",1000,0.,10.);
  TH2F *hPtE2D = new TH2F("PtE2D","",1000,0.,10.,1000,0.,10.);
  TH2F *hPtMu2D = new TH2F("PtMu2D","",1000,0.,10.,1000,0.,10.);
  TH1F *hEvent = new TH1F("event","",10,0.,10.);

  TH2F *hdNdpTEMid = new TH2F("dNdpTEMid","",1000,0.,10.,100,0.,10.);
  TH2F *hdNdpTMuMid = new TH2F("dNdpTMuMid","",1000,0.,10.,100,0.,10.);
  TH2F *hdNdmTEMid = new TH2F("dNdmTEMid","",1000,0.,10.,100,0.,10.);
  TH2F *hdNdmTMuMid = new TH2F("dNdmTMuMid","",1000,0.,10.,100,0.,10.);


  /////////////////////////////
  // Initialize the PYTHIA
  ////////////////////////////
  TPythia6 *pythia = TPythia6::Instance();

  pythia->SetMSEL(0);
  pythia->SetMSUB(1,1);
  pythia->SetMSTP(43,1);
  // switch off all channels
  for(int i=174;i<=189;i++) {
    pythia->SetMDME(i,1,0);
  }
  pythia->SetMDME(182,1,1);
  pythia->SetMDME(184,1,1);
  pythia->SetCKIN(1,CKIN1);
  pythia->SetCKIN(2,CKIN2);

  if(strcmp(Energy,"RHIC")==0) {
    pythia->Initialize("CMS","p","p",200.);
  } else if(strcmp(Energy,"SPS")==0) {
    pythia->Initialize("CMS","p","p",17.2);
  } else {
    cout << " Unknown collision energy! EXIT!!!" << endl;
    return;
  }


//  pythia->Pylist(12);
  /////////////////////////
  // start the event loop
  /////////////////////////

  for (Int_t i_evt=0; i_evt<Nevts; i_evt++) {
    pythia->GenerateEvent();
    hEvent->Fill(0);

    if(i_evt<10) pythia->Pylist(1);

    Int_t n_part=pythia->GetN();

    bool findEP = false;
    bool findEM = false;
    bool findMuP = false;
    bool findMuM = false;

    TLorentzVector pos(0.,0.,0.);
    TLorentzVector neg(0.,0.,0.);
    for (Int_t i_part=1; i_part<=n_part; i_part++) {
      if(pythia->GetK(pythia->GetK(pythia->GetK(i_part,3),3),2)!=23) continue;
      if(pythia->GetK(i_part,1)!=1) continue;

      if(pythia->GetK(i_part,2)==11) {        
        neg.SetXYZM(pythia->GetP(i_part,1), pythia->GetP(i_part,2),
pythia->GetP(i_part,3),0.000511);
        findEM = true;
      } else if(pythia->GetK(i_part,2)==-11) {
        pos.SetXYZM(pythia->GetP(i_part,1), pythia->GetP(i_part,2),
pythia->GetP(i_part,3),0.000511);
        findEP = true;
      }
      if(pythia->GetK(i_part,2)==13) {        
        neg.SetXYZM(pythia->GetP(i_part,1), pythia->GetP(i_part,2),
pythia->GetP(i_part,3),0.10566);
        findMuM = true;
      } else if(pythia->GetK(i_part,2)==-13) {
        pos.SetXYZM(pythia->GetP(i_part,1), pythia->GetP(i_part,2),
pythia->GetP(i_part,3),0.10566);
        findMuP = true;
      }
    }

    TLorentzVector pair = neg + pos;

    if(findEM && findEP) {
        hPtE2D->Fill(neg.Pt(), pos.Pt());
        hInvMassE->Fill(pair.M());
        if(fabs(neg.Eta())<1.&&fabs(pos.Eta())<1.) hInvMassEMid->Fill(pair.M());

        if(fabs(pair.Rapidity())<0.5) {
           hdNdpTEMid->Fill(pair.Perp(),pair.M());
           hdNdmTEMid->Fill(pair.Mt()-pair.M(),pair.M(),1./pair.Mt());
        }
        hEvent->Fill(1);
    }
    if(findMuM && findMuP) {
        hPtMu2D->Fill(neg.Pt(), pos.Pt());
        hInvMassMu->Fill(pair.M());
        if(fabs(neg.Eta())<1.&&fabs(pos.Eta())<1.) hInvMassMuMid->Fill(pair.M());
        if(fabs(pair.Rapidity())<0.5) {
           hdNdpTMuMid->Fill(pair.Perp(),pair.M());
           hdNdmTMuMid->Fill(pair.Mt()-pair.M(),pair.M(),1./pair.Mt());
        }
        hEvent->Fill(2);
    }


    if (i_evt%1000==0)
      cout << "Event " << i_evt << endl;
  }

  pythia->Pystat(1);

  //////////////////////////
  // output the histograms
  //////////////////////////
  char outname[512];
  sprintf(outname,"%s_%s_%d_%d.root",OutName,Energy,(Int_t)CKIN1,(Int_t)CKIN2);
  TFile fout(outname,"RECREATE");
  fout.cd();
  hEvent->Write();
  hInvMassE->Write();
  hInvMassMu->Write();
  hInvMassEMid->Write();
  hInvMassMuMid->Write();
  hPtE2D->Write();
  hPtMu2D->Write();

  hdNdpTEMid->Write();
  hdNdmTEMid->Write();
  hdNdpTMuMid->Write();
  hdNdmTMuMid->Write();
  fout.Close();
  return;
}
