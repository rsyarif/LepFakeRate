#include "interface/utilities.h"
#include "interface/CMSStyle.C"
#include <TH2.h>
#include <TH1F.h>
#include "TFile.h"
#include <iostream>


void fakeratio(){
    TFile* f = TFile::Open("fakeeventloop_Data2lep.root");
    f->cd();

    //histograms are loose but not tight. 
    TH1F* h_tightloose_CR1_ele = (TH1F*)f->Get("h_tightloose_CR1_ele");
    TH1F* h_lepEta_Loose_CR1_ele = (TH1F*)f->Get("h_lepEta_Loose_CR1_ele");
    TH1F* h_lepEta_Tight_CR1_ele = (TH1F*)f->Get("h_lepEta_Tight_CR1_ele");
    TH1F* h_lepPt_Loose_CR1_ele = (TH1F*)f->Get("h_lepPt_Loose_CR1_ele");
    TH1F* h_lepPt_Tight_CR1_ele = (TH1F*)f->Get("h_lepPt_Tight_CR1_ele");
    TH2F* h_lepEta_vs_lepPt_Loose_CR1_ele = (TH2F*)f->Get("h_lepEta_vs_lepPt_Loose_CR1_ele");
    TH2F* h_lepEta_vs_lepPt_Tight_CR1_ele = (TH2F*)f->Get("h_lepEta_vs_lepPt_Tight_CR1_ele");
    TH2F* h_lepEta_vs_coneCorPt_Loose_CR1_ele = (TH2F*)f->Get("h_lepEta_vs_coneCorPt_Loose_CR1_ele");
    TH2F* h_lepEta_vs_coneCorPt_Tight_CR1_ele = (TH2F*)f->Get("h_lepEta_vs_coneCorPt_Tight_CR1_ele");

    TH1F* h_tightloose_CR1_mu = (TH1F*)f->Get("h_tightloose_CR1_mu");
    TH1F* h_lepEta_Loose_CR1_mu = (TH1F*)f->Get("h_lepEta_Loose_CR1_mu");
    TH1F* h_lepEta_Tight_CR1_mu = (TH1F*)f->Get("h_lepEta_Tight_CR1_mu");
    TH1F* h_lepPt_Loose_CR1_mu = (TH1F*)f->Get("h_lepPt_Loose_CR1_mu");
    TH1F* h_lepPt_Tight_CR1_mu = (TH1F*)f->Get("h_lepPt_Tight_CR1_mu");
    TH2F* h_lepEta_vs_lepPt_Loose_CR1_mu = (TH2F*)f->Get("h_lepEta_vs_lepPt_Loose_CR1_mu");
    TH2F* h_lepEta_vs_lepPt_Tight_CR1_mu = (TH2F*)f->Get("h_lepEta_vs_lepPt_Tight_CR1_mu");
    TH2F* h_lepEta_vs_coneCorPt_Loose_CR1_mu = (TH2F*)f->Get("h_lepEta_vs_coneCorPt_Loose_CR1_mu");
    TH2F* h_lepEta_vs_coneCorPt_Tight_CR1_mu = (TH2F*)f->Get("h_lepEta_vs_coneCorPt_Tight_CR1_mu");


    h_lepEta_Loose_CR1_ele->Add( h_lepEta_Tight_CR1_ele);
    h_lepPt_Loose_CR1_ele->Add( h_lepPt_Tight_CR1_ele);
    h_lepEta_vs_lepPt_Loose_CR1_ele->Add( h_lepEta_vs_lepPt_Tight_CR1_ele);
    h_lepEta_vs_coneCorPt_Loose_CR1_ele->Add( h_lepEta_vs_coneCorPt_Tight_CR1_ele);

    h_lepEta_Loose_CR1_mu->Add( h_lepEta_Tight_CR1_mu);
    h_lepPt_Loose_CR1_mu->Add( h_lepPt_Tight_CR1_mu);
    h_lepEta_vs_lepPt_Loose_CR1_mu->Add( h_lepEta_vs_lepPt_Tight_CR1_mu);
    h_lepEta_vs_coneCorPt_Loose_CR1_mu->Add( h_lepEta_vs_coneCorPt_Tight_CR1_mu);
    enum tightloose{loose = 1, tight = 2};
    //tightlooose has bin 0 = Loose but not tight, bin 1 = tight. 
    //We change this to bin 0 = inclusive loose and bin 1 = tight. 

    cout<<"******************* Electron Fake Rates **********************"<<endl;
    h_tightloose_CR1_ele->SetBinContent(loose,h_tightloose_CR1_ele->GetBinContent(loose) + h_tightloose_CR1_ele->GetBinContent(tight) );
    float fe = h_tightloose_CR1_ele->GetBinContent(tight)/h_tightloose_CR1_ele->GetBinContent(loose);
    float ufe = fe*sqrt( (1.0/h_tightloose_CR1_ele->GetBinContent(tight)) + (1.0/h_tightloose_CR1_ele->GetBinContent(loose)) );//sqrt(N) stat uncert.
    //uncert(A/B) = (A/B)*sqrt( 1/A + 1/B);
    printf("Flat Electron Fake Rate: %.3f +- %.3f. From %i tight and %i inclusive loose\n",fe,ufe,(int)h_tightloose_CR1_ele->GetBinContent(tight),(int)h_tightloose_CR1_ele->GetBinContent(loose));

    TH1F* h_lepEta_CR1_ele = DivideHists( h_lepEta_Tight_CR1_ele, h_lepEta_Loose_CR1_ele);
    TH1F* h_lepPt_CR1_ele  = DivideHists( h_lepPt_Tight_CR1_ele , h_lepPt_Loose_CR1_ele );
    TH2F* h_lepEta_vs_lepPt_CR1_ele = DivideHists( h_lepEta_vs_lepPt_Tight_CR1_ele, h_lepEta_vs_lepPt_Loose_CR1_ele);
    TH2F* h_lepEta_vs_coneCorPt_CR1_ele = DivideHists( h_lepEta_vs_coneCorPt_Tight_CR1_ele, h_lepEta_vs_coneCorPt_Loose_CR1_ele);

    cout<<"******************* Muon Fake Rates **********************"<<endl;
    h_tightloose_CR1_mu->SetBinContent(loose,h_tightloose_CR1_mu->GetBinContent(loose) + h_tightloose_CR1_mu->GetBinContent(tight) );
    float fm = h_tightloose_CR1_mu->GetBinContent(tight)/h_tightloose_CR1_mu->GetBinContent(loose);
    float ufm = fm*sqrt( (1.0/h_tightloose_CR1_mu->GetBinContent(tight)) + (1.0/h_tightloose_CR1_mu->GetBinContent(loose)) );//sqrt(N) stat uncert.
    //uncert(A/B) = (A/B)*sqrt( 1/A + 1/B);
    printf("Flat Muon Fake Rate: %.3f +- %.3f. From %i tight and %i inclusive loose\n",fm,ufm,(int)h_tightloose_CR1_mu->GetBinContent(tight),(int)h_tightloose_CR1_mu->GetBinContent(loose));

    TH1F* h_lepEta_CR1_mu = DivideHists( h_lepEta_Tight_CR1_mu, h_lepEta_Loose_CR1_mu);
    TH1F* h_lepPt_CR1_mu  = DivideHists( h_lepPt_Tight_CR1_mu , h_lepPt_Loose_CR1_mu );
    TH2F* h_lepEta_vs_lepPt_CR1_mu = DivideHists( h_lepEta_vs_lepPt_Tight_CR1_mu, h_lepEta_vs_lepPt_Loose_CR1_mu);
    TH2F* h_lepEta_vs_coneCorPt_CR1_mu = DivideHists( h_lepEta_vs_coneCorPt_Tight_CR1_mu, h_lepEta_vs_coneCorPt_Loose_CR1_mu);

    TFile* out = new TFile("fakeevent_result.root","RECREATE");
    out->cd();
    h_tightloose_CR1_ele->Write();
    h_tightloose_CR1_mu->Write();
    h_lepEta_CR1_ele->Write();
    h_lepPt_CR1_ele->Write();
    h_lepEta_vs_lepPt_CR1_ele->Write();
    h_lepEta_vs_coneCorPt_CR1_ele->Write();
    h_lepEta_CR1_mu->Write();
    h_lepPt_CR1_mu->Write();
    h_lepEta_vs_lepPt_CR1_mu->Write();
    h_lepEta_vs_coneCorPt_CR1_mu->Write();
    out->Close();
    f->Close();
}

/*  	h_tightloose_CR1_ele
	h_lepEta_Loose_CR1_ele
	h_lepEta_Tight_CR1_ele
	h_lepPt_Loose_CR1_ele
	h_lepPt_Tight_CR1_ele
	h_lepEta_vs_lepPt_Loose_CR1_ele
	h_lepEta_vs_lepPt_Tight_CR1_ele
	h_lepEta_vs_coneCorPt_Loose_CR1_ele
	h_lepEta_vs_coneCorPt_Tight_CR1_ele
	h_tightloose_CR1_mu
	h_lepEta_Loose_CR1_mu
	h_lepEta_Tight_CR1_mu
	h_lepPt_Loose_CR1_mu
	h_lepPt_Tight_CR1_mu
	h_lepEta_vs_lepPt_Loose_CR1_mu
	h_lepEta_vs_lepPt_Tight_CR1_mu
	h_lepEta_vs_coneCorPt_Loose_CR1_mu
	h_lepEta_vs_coneCorPt_Tight_CR1_mu*/
