#define fakeratecalc_cxx
#include "interface/fakeratecalc.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include "TProfile.h"
#include <map>
#include <set>
#include <cmath>
#include "TMath.h"
#include <algorithm>
#include <utility>
//my additions
#include "TAxis.h"
#include <TH3.h>
#include "TFile.h"
#include "TTree.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include <assert.h>
#include "TF1.h"
#include <iostream>
#include <fstream>
#include <string.h>
#include "TVector3.h"
//#include <algorithm>
#include "TLorentzVector.h"
#include "interface/CMSStyle.C"
#include <time.h>
#include "interface/stringmap.h"
#include "interface/DMCblock_fakerates.h"
#include "interface/DMCblock_list_fakerates.h"
#include "interface/controlpannel.h"
#include "interface/KinematicVar.h"
#include "eventRegistry.C"
#include "interface/EleIDCuts.h"
#include "check_dup.C"
//#include "EventTopology.h"

//   In a ROOT session, you can do:
//      root> .L fakeratecalc.C
//      root> fakeratecalc* t = new fakeratecalc("blockname")
//      root> t->Loop();       // Loop on all entries
//	Loops over a DMC block. 
//	
//	Origionally:
//      root> .L fakeratecalc.C
//      root> fakeratecalc t
//      root> t.Loop();       // Loop on all entries

using namespace std;

template<typename T> bool EtGreater(const T* p1, const T* p2) {
    return (p1->momentum.Et() > p2->momentum.Et());
}

typedef stringmap<TH1F*> LabelHist;
typedef stringmap<LabelHist> Label2Hist;
typedef stringmap<TH2F*> LabelHist2D;
typedef stringmap<LabelHist2D> Label2Hist2D;

class counters: public stringmap<int>{
    public:
	vector <string> order;
	void addCounter(string Label){
	    set(Label,0);
	    order.push_back(Label);
	}
	int increment(string Label){
	    if(tmap.find(Label) != tmap.end()) return ++tmap[Label];
	    else{
		cerr<<"Warning! in counters.increment received Invalid key string "<<Label<< endl;
		return -1;
	    }
	}
};

TH1F* init_TH1F_special(string label, string var, LabelKinVars &KinVars ); //cant throw  <std::string,int> errorpair
TH1F* init_TH1F_special(string label, KinematicVar* k); 
TH2F* init_TH2F_special(string label, KinematicVar* kx, KinematicVar* ky);
void MakeTriggerIndicies(std::vector<int> & TriggerIndicies, int NTriggers, string Triggers[], vector<string>  * vsSelTriggers);
std::map<int,int> Make_HLT_map_selIndex_to_stringIndex(int NTriggers, string Triggers[], vector<string>  * vsSelTriggers);
void Fill_HLT_study(TH2F* h, std::map<int,int> indexmap, int nTrig, std::vector<int> TriggerIndicies, vector<int> *viSelTriggers, bool verbose);
bool checkTriggerList(std::vector<int> & TriggerIndicies, int NTriggers, string Triggers[], vector<string>  * vsSelTriggers);
void indexPtSort(vector<int>& indicies, int n, vector<double> pt);
//bool pairCompare(const std::pair<int, double>& firstElem, const std::pair<int, double>& secondElem);
struct mypair{ int index; double pt; };
bool pairCompare(mypair* firstElem, mypair* secondElem);

//      __                    
//     / /   ____  ____  ____ 
//    / /   / __ \/ __ \/ __ \
//   / /___/ /_/ / /_/ / /_/ /
//  /_____/\____/\____/ .___/ 
//                   /_/      
//
void fakeratecalc::Loop(eventRegistry* EventRegistry,eventRegistry* BadEventRegistry)
{
    cout<<"Loop lives"<<endl;
    ofstream firealarm;
    bool onfire = false; //set to true in case of alarming error.
    firealarm.open("tmpfirealarm.txt");
    firealarm << "fakeratecalc is on FIRE"<<endl;
    firealarm.close();

    const bool isSignal = dmcblock->type == 1 or dmcblock->type == 2;
    ///////////////////////////////////////////////////////
    //control pannel; dashboard
    static const int printlevel=3;//print level printlevel printlevel
    //0 = Stealth
    //1 = Annouce birth, death, marriage, and explosions
    //2 = Christmas and easter only; turn on 1000th event heart monitor, Normal for running
    //3 = Casual:         Tick off major blocks and probe suspicious lines
    //***Begin compartmentalizing the verbosity with print blocks***
    //4 = Currious:	 Tick off minor blocks
    //5 = Probe loops 1 level deep and chosen significant lines
    //6 = Probe loops 2 levels deep
    //7 = Probe loops 3+ levels deep
    //8 = Paranoid:        Probe every significant operation that looks like it could fail
    //9 = FBI Colonoscopy: Probe every operation, even if it's unlikely that it'll fail
    //10= Hyperbarf:       Print a probe after every line
    //11= ALICE: 	       Blast appart its nuclei and examine the entrails 

    ///////////////////////////////////////////////////////
    if(printlevel > 3) std::cout << "Initialize event counters." << std::endl;	
    //vector <string> Counter_order;
    counters Counters;
    Counters.addCounter("no cuts");
    Counters.addCounter("pre HLT");
    Counters.addCounter("post HLT");
    Counters.addCounter("Has 2 leptons");
    Counters.addCounter("Has exactly 1 lepton");
    namedbool* KinVarSwitches = MakeKinVarSwitches(); //out of all the kinvars you've imagined, which are to be used? 

    if(printlevel >= 5) std::cout << "counters started" << std::endl;	

    //namedbool* AlltopoCut = new namedbool();
    namedint AlltopoNpass;
    //namedbool* topoCut = new namedbool();
    if(printlevel >= 9) std::cout << "made namedbool" << std::endl;	

    for (int iTopo = 0; iTopo<nEventTopologies; iTopo++) {
	//topoCut->set(s_EventTopology[iTopo], false);
	AlltopoNpass.set(s_EventTopology[iTopo], 0);
    }
    if(printlevel >= 9) std::cout << "done setting topoCuts" << std::endl;	


    string s_forTopo[nEventTopologies];
    for (int iTopo=0; iTopo<nEventTopologies; iTopo++){
	s_forTopo[iTopo]=string(" for ")+s_EventTopology[iTopo];
    }

    typedef stringmap<KinematicVar*>* LabelKinVars;
    LabelKinVars KinVars = setupKinematicVar();

    if(printlevel >= 4) cout<<"initializing specials"<<endl;
    ///////////////////Initalize Specials//////////////////////
    TH1F* h_tightloose_CR1_ele;
    TH1F* h_lepEta_Loose_CR1_ele;
    TH1F* h_lepEta_Tight_CR1_ele;
    TH1F* h_lepPt_Loose_CR1_ele;
    TH1F* h_lepPt_Tight_CR1_ele;
    TH2F* h_lepEta_vs_lepPt_Loose_CR1_ele;
    TH2F* h_lepEta_vs_lepPt_Tight_CR1_ele;
    TH2F* h_lepEta_vs_coneCorPt_Loose_CR1_ele;
    TH2F* h_lepEta_vs_coneCorPt_Tight_CR1_ele;

    TH1F* h_tightloose_CR1_mu;
    TH1F* h_lepEta_Loose_CR1_mu;
    TH1F* h_lepEta_Tight_CR1_mu;
    TH1F* h_lepPt_Loose_CR1_mu;
    TH1F* h_lepPt_Tight_CR1_mu;
    TH2F* h_lepEta_vs_lepPt_Loose_CR1_mu;
    TH2F* h_lepEta_vs_lepPt_Tight_CR1_mu;
    TH2F* h_lepEta_vs_coneCorPt_Loose_CR1_mu;
    TH2F* h_lepEta_vs_coneCorPt_Tight_CR1_mu;

    TH1F* h_tightloose_CR_T53_ele;
    TH1F* h_lepEta_Loose_CR_T53_ele;
    TH1F* h_lepEta_Tight_CR_T53_ele;
    TH1F* h_lepPt_Loose_CR_T53_ele;
    TH1F* h_lepPt_Tight_CR_T53_ele;
    TH2F* h_lepEta_vs_lepPt_Loose_CR_T53_ele;
    TH2F* h_lepEta_vs_lepPt_Tight_CR_T53_ele;
    TH2F* h_lepEta_vs_coneCorPt_Loose_CR_T53_ele;
    TH2F* h_lepEta_vs_coneCorPt_Tight_CR_T53_ele;

    TH1F* h_tightloose_CR_T53_mu;
    TH1F* h_lepEta_Loose_CR_T53_mu;
    TH1F* h_lepEta_Tight_CR_T53_mu;
    TH1F* h_lepPt_Loose_CR_T53_mu;
    TH1F* h_lepPt_Tight_CR_T53_mu;
    TH2F* h_lepEta_vs_lepPt_Loose_CR_T53_mu;
    TH2F* h_lepEta_vs_lepPt_Tight_CR_T53_mu;
    TH2F* h_lepEta_vs_coneCorPt_Loose_CR_T53_mu;
    TH2F* h_lepEta_vs_coneCorPt_Tight_CR_T53_mu;

    try{
	KinematicVar* ktightloose = KinVars->get_throwable("tightloose",3);
	h_tightloose_CR1_ele = init_TH1F_special("h_tightloose_CR1_ele", ktightloose); 
	h_tightloose_CR1_mu = init_TH1F_special("h_tightloose_CR1_mu", ktightloose); 

        KinematicVar* klepEta = KinVars->get_throwable("lepEta",3);
	h_lepEta_Loose_CR1_ele = init_TH1F_special("h_lepEta_Loose_CR1_ele", klepEta); 
	h_lepEta_Tight_CR1_ele = init_TH1F_special("h_lepEta_Tight_CR1_ele", klepEta); 
	h_lepEta_Loose_CR1_mu = init_TH1F_special("h_lepEta_Loose_CR1_mu", klepEta); 
	h_lepEta_Tight_CR1_mu = init_TH1F_special("h_lepEta_Tight_CR1_mu", klepEta); 

        KinematicVar* klepPt = KinVars->get_throwable("lepPt",3);
	h_lepPt_Loose_CR1_ele = init_TH1F_special("h_lepPt_Loose_CR1_ele", klepPt); 
	h_lepPt_Tight_CR1_ele = init_TH1F_special("h_lepPt_Tight_CR1_ele", klepPt); 
	h_lepPt_Loose_CR1_mu = init_TH1F_special("h_lepPt_Loose_CR1_mu", klepPt); 
	h_lepPt_Tight_CR1_mu = init_TH1F_special("h_lepPt_Tight_CR1_mu", klepPt); 

	h_lepEta_vs_lepPt_Loose_CR1_ele = init_TH2F_special("h_lepEta_vs_lepPt_Loose_CR1_ele", klepPt, klepEta);
	h_lepEta_vs_lepPt_Tight_CR1_ele = init_TH2F_special("h_lepEta_vs_lepPt_Tight_CR1_ele", klepPt, klepEta);
	h_lepEta_vs_lepPt_Loose_CR1_mu = init_TH2F_special("h_lepEta_vs_lepPt_Loose_CR1_mu", klepPt, klepEta);
	h_lepEta_vs_lepPt_Tight_CR1_mu = init_TH2F_special("h_lepEta_vs_lepPt_Tight_CR1_mu", klepPt, klepEta);

        KinematicVar* kconeCorrPt = KinVars->get_throwable("coneCorrPt",3);
	h_lepEta_vs_coneCorPt_Loose_CR1_ele = init_TH2F_special("h_lepEta_vs_coneCorPt_Loose_CR1_ele", kconeCorrPt, klepEta);
	h_lepEta_vs_coneCorPt_Tight_CR1_ele = init_TH2F_special("h_lepEta_vs_coneCorPt_Tight_CR1_ele", kconeCorrPt, klepEta);
	h_lepEta_vs_coneCorPt_Loose_CR1_mu = init_TH2F_special("h_lepEta_vs_coneCorPt_Loose_CR1_mu", kconeCorrPt, klepEta);
	h_lepEta_vs_coneCorPt_Tight_CR1_mu = init_TH2F_special("h_lepEta_vs_coneCorPt_Tight_CR1_mu", kconeCorrPt, klepEta);



	h_tightloose_CR_T53_ele = init_TH1F_special("h_tightloose_CR_T53_ele", ktightloose); 
	h_tightloose_CR_T53_mu = init_TH1F_special("h_tightloose_CR_T53_mu", ktightloose); 

	h_lepEta_Loose_CR_T53_ele = init_TH1F_special("h_lepEta_Loose_CR_T53_ele", klepEta); 
	h_lepEta_Tight_CR_T53_ele = init_TH1F_special("h_lepEta_Tight_CR_T53_ele", klepEta); 
	h_lepEta_Loose_CR_T53_mu = init_TH1F_special("h_lepEta_Loose_CR_T53_mu", klepEta); 
	h_lepEta_Tight_CR_T53_mu = init_TH1F_special("h_lepEta_Tight_CR_T53_mu", klepEta); 

	h_lepPt_Loose_CR_T53_ele = init_TH1F_special("h_lepPt_Loose_CR_T53_ele", klepPt); 
	h_lepPt_Tight_CR_T53_ele = init_TH1F_special("h_lepPt_Tight_CR_T53_ele", klepPt); 
	h_lepPt_Loose_CR_T53_mu = init_TH1F_special("h_lepPt_Loose_CR_T53_mu", klepPt); 
	h_lepPt_Tight_CR_T53_mu = init_TH1F_special("h_lepPt_Tight_CR_T53_mu", klepPt); 

	h_lepEta_vs_lepPt_Loose_CR_T53_ele = init_TH2F_special("h_lepEta_vs_lepPt_Loose_CR_T53_ele", klepPt, klepEta);
	h_lepEta_vs_lepPt_Tight_CR_T53_ele = init_TH2F_special("h_lepEta_vs_lepPt_Tight_CR_T53_ele", klepPt, klepEta);
	h_lepEta_vs_lepPt_Loose_CR_T53_mu = init_TH2F_special("h_lepEta_vs_lepPt_Loose_CR_T53_mu", klepPt, klepEta);
	h_lepEta_vs_lepPt_Tight_CR_T53_mu = init_TH2F_special("h_lepEta_vs_lepPt_Tight_CR_T53_mu", klepPt, klepEta);

	h_lepEta_vs_coneCorPt_Loose_CR_T53_ele = init_TH2F_special("h_lepEta_vs_coneCorPt_Loose_CR_T53_ele", kconeCorrPt, klepEta);
	h_lepEta_vs_coneCorPt_Tight_CR_T53_ele = init_TH2F_special("h_lepEta_vs_coneCorPt_Tight_CR_T53_ele", kconeCorrPt, klepEta);
	h_lepEta_vs_coneCorPt_Loose_CR_T53_mu = init_TH2F_special("h_lepEta_vs_coneCorPt_Loose_CR_T53_mu", kconeCorrPt, klepEta);
	h_lepEta_vs_coneCorPt_Tight_CR_T53_mu = init_TH2F_special("h_lepEta_vs_coneCorPt_Tight_CR_T53_mu", kconeCorrPt, klepEta);
        //                h_ST_SSDLsansS  = init_TH1F_special( "h_STSSDL_sansST", kST);//change this; this name makes no sense
        //                        h_ST_OSDL1sansS = init_TH1F_special("h_ST_OSDL1sansS", kST);
        //                                h_ST_OSDL1sansSlb       = init_TH1F_special("h_ST_OSDL1sansSlb", kST);
        //                                        h_ST_OSDL1sansSH        = init_TH1F_special("h_ST_OSDL1sansSH", kST);
        //                                                h_ST_OSDL2sansS = init_TH1F_special("h_ST_OSDL2sansS", kST);
        //                                                        h_ST_OSDL2sansSH        = init_TH1F_special("h_ST_OSDL2sansSH", kST);
        //                                                                h_ST_SSDLsansSH = init_TH1F_special("h_ST_SSDLsansSH", kST);
        //
        //
   }//end try
    catch(std::pair <std::string,int> errorpair){
        cerr<<"Error while initalizing special histograms; Invalid kinematic variable key string \""<<errorpair.first<<"\" sought in KinVars"<<endl;
        std::terminate();
    } //end catch failed stringmap gets

    //define an array of trigger names.
    //Ele-mu triggers live in both lists; it doesn't matter which list you put them in. 
    //but be careful to only put double electron triggers in the electron list and muon triggers in the muon list else weird failures. 
    const int NEleTriggers = 5;
    //string EleTriggers[NEleTriggers ] = {"HLT_Ele17_CaloIdL_TrackIdL_IsoVL_v"};
    string EleTriggers[NEleTriggers ] = {
     "HLT_Ele8_CaloIdM_TrackIdM_PFJet30_v", "HLT_Ele12_CaloIdM_TrackIdM_PFJet30_v", "HLT_Ele18_CaloIdM_TrackIdM_PFJet30_v", "HLT_Ele23_CaloIdM_TrackIdM_PFJet30_v", "HLT_Ele33_CaloIdM_TrackIdM_PFJet30_v"}; //your main.


//        "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_v", "HLT_Ele18_CaloIdL_TrackIdL_IsoVL_PFJet30_v", "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_v", "HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_v"};


    const int NMuTriggers = 4;
    string MuTriggers[NMuTriggers ] = {
        "HLT_Mu8_v", "HLT_Mu17_v", "HLT_Mu24_v", "HLT_Mu34_v"};
//        "HLT_Mu8_TrkIsoVVL_v", "HLT_Mu24_TrkIsoVVL_v", "HLT_Mu34_TrkIsoVVL_v"};

	//Initialize the HLT study histograms and set their bin labels. 
    TH2F* HLT_coincidience_ele_CR1 = new TH2F("HLT_coincidience_ele_CR1","HLT_coincidience_ele_CR1",NEleTriggers+1,0,NEleTriggers+1,NEleTriggers,0,NEleTriggers);
    TH2F* HLT_coincidience_mu_CR1 = new TH2F("HLT_coincidience_mu_CR1","HLT_coincidience_mu_CR1",NMuTriggers+1,0,NMuTriggers+1,NMuTriggers,0,NMuTriggers);
    TH2F* HLT_coincidience_ele = new TH2F("HLT_coincidience_ele","HLT_coincidience_ele",NEleTriggers+1,0,NEleTriggers+1,NEleTriggers,0,NEleTriggers);
    TH2F* HLT_coincidience_mu = new TH2F("HLT_coincidience_mu","HLT_coincidience_mu",NMuTriggers+1,0,NMuTriggers+1,NMuTriggers,0,NMuTriggers);
    for (int i=0;i<NEleTriggers;i++){ 
	HLT_coincidience_ele_CR1->GetXaxis()->SetBinLabel(i+1,EleTriggers[i].c_str());
	HLT_coincidience_ele_CR1->GetYaxis()->SetBinLabel(i+1,EleTriggers[i].c_str());
	HLT_coincidience_ele->GetXaxis()->SetBinLabel(i+1,EleTriggers[i].c_str());
	HLT_coincidience_ele->GetYaxis()->SetBinLabel(i+1,EleTriggers[i].c_str());
    }
    HLT_coincidience_ele_CR1->GetXaxis()->SetBinLabel(NEleTriggers+1,"uniquely");
    HLT_coincidience_ele->GetXaxis()->SetBinLabel(NEleTriggers+1,"uniquely");
    for (int i=0;i<NMuTriggers;i++){ 
	HLT_coincidience_mu_CR1->GetXaxis()->SetBinLabel(i+1,MuTriggers[i].c_str());
	HLT_coincidience_mu_CR1->GetYaxis()->SetBinLabel(i+1,MuTriggers[i].c_str());
	HLT_coincidience_mu->GetXaxis()->SetBinLabel(i+1,MuTriggers[i].c_str());
	HLT_coincidience_mu->GetYaxis()->SetBinLabel(i+1,MuTriggers[i].c_str());
    }
    HLT_coincidience_mu_CR1->GetXaxis()->SetBinLabel(NMuTriggers+1,"uniquely");
    HLT_coincidience_mu->GetXaxis()->SetBinLabel(NMuTriggers+1,"uniquely");

    int nPassDilepCut = 0;
    int nDuplicates = 0;

    int nPos = 0;
    int nNeg = 0;

    //     __  ______   _____  __  __   ____  ____  ___ 
    //    /  |/  / _ | /  _/ |/ / / /  / __ \/ __ \/ _ \
    //   / /|_/ / __ |_/ //    / / /__/ /_/ / /_/ / ___/
    //  /_/  /_/_/ |_/___/_/|_/ /____/\____/\____/_/    
    //                                                  
    if(printlevel >= 5) std::cout << "Approach Main Loop" << std::endl;	
    //MAIN EVENT LOOP
    if (fChain == 0) return;
    bool first_event = true;
    vector<int> MuTriggerIndicies;
    vector<int> EleTriggerIndicies;
    std::map<int,int> HLT_map_ele;
    std::map<int,int> HLT_map_mu;

    if(printlevel >= 6) std::cout << "try GetEntriesFast" << std::endl;	
    Long64_t nentries = fChain->GetEntriesFast();//this segfaults. 
    if(printlevel >= 6) std::cout << "did it" << std::endl;	
    //for (Long64_t jentry=0; jentry<1000;jentry++) {
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
	if( jentry%1000 == 0 && printlevel >=3 ) printf("%i out of %i events processing\n",int(jentry), int(nentries));

	if(printlevel >= 6) std::cout << "Attempt LoadTree" << std::endl;	
	Long64_t ientry = LoadTree(jentry);
	if(printlevel >= 6) std::cout << "fin LoadTree, attempt GetEntry" << std::endl;	
	if (ientry < 0) break;
	if(GetEntry(jentry) == 0) break;  
	//fChain->GetEntry(jentry);   
	if(printlevel >= 6) std::cout << "fin attempt GetEntry" << std::endl;	


	Counters.increment("no cuts");

	Counters.increment("pre HLT");

	if(first_event){ //search through the string of trigger indicies 
	    if(dmcblock->isMC){
		MakeTriggerIndicies(EleTriggerIndicies, NEleTriggers, EleTriggers, vsSelMCTriggersEl_singleLepCalc);
		MakeTriggerIndicies(MuTriggerIndicies, NMuTriggers, MuTriggers, vsSelMCTriggersMu_singleLepCalc);

		//cout<<"Running the map maker MC"<<endl;
		HLT_map_ele = Make_HLT_map_selIndex_to_stringIndex(NEleTriggers, EleTriggers, vsSelMCTriggersEl_singleLepCalc);
		HLT_map_mu = Make_HLT_map_selIndex_to_stringIndex(NMuTriggers, MuTriggers, vsSelMCTriggersMu_singleLepCalc);
	    }
	    else{
		MakeTriggerIndicies(EleTriggerIndicies, NEleTriggers, EleTriggers, vsSelTriggersEl_singleLepCalc);
		MakeTriggerIndicies(MuTriggerIndicies, NMuTriggers, MuTriggers, vsSelTriggersMu_singleLepCalc);

		//cout<<"Running the map maker"<<endl;
		HLT_map_ele = Make_HLT_map_selIndex_to_stringIndex(NEleTriggers, EleTriggers, vsSelTriggersEl_singleLepCalc);
		HLT_map_mu =  Make_HLT_map_selIndex_to_stringIndex(NMuTriggers,  MuTriggers,  vsSelTriggersMu_singleLepCalc);
	    }
	} //end compile trigger indicies

	//occasionally check that the trigger location really is consistent. 
	if( jentry%1000 == 0){ 
	    //if( jentry%1000 == 0 && !first_event ){ 

	    //look at the trigger list
	    bool trigs_match = true;
	    if(dmcblock->isMC){
		trigs_match &= checkTriggerList(EleTriggerIndicies, NEleTriggers, EleTriggers, vsSelMCTriggersEl_singleLepCalc);
		trigs_match &= checkTriggerList(MuTriggerIndicies, NMuTriggers, MuTriggers, vsSelMCTriggersMu_singleLepCalc);
	    }
	    else{
		trigs_match &= checkTriggerList(EleTriggerIndicies, NEleTriggers, EleTriggers, vsSelTriggersEl_singleLepCalc);
		trigs_match &= checkTriggerList(MuTriggerIndicies, NMuTriggers, MuTriggers, vsSelTriggersMu_singleLepCalc);
	    }
	    if(!trigs_match){
		cout<<"ERROR! trigger list no longer matches what it started as. Your trigger selection is probably broken!"<<endl;
		std::terminate();
	    }

	} //end occasionally
	first_event = false;

	bool passHLT = false;
	if(dmcblock->isMC){
	    for(std::vector<int>::iterator i = EleTriggerIndicies.begin(); i != EleTriggerIndicies.end() && !passHLT; i++)
		passHLT |= (*viSelMCTriggersEl_singleLepCalc)[*i];
	    for(std::vector<int>::iterator i = MuTriggerIndicies.begin(); i != MuTriggerIndicies.end() && !passHLT; i++)
		passHLT |= (*viSelMCTriggersMu_singleLepCalc)[*i];
	}
	else{
	    for(std::vector<int>::iterator i = EleTriggerIndicies.begin(); i != EleTriggerIndicies.end() && !passHLT; i++)
		passHLT |= (*viSelTriggersEl_singleLepCalc)[*i];
	    for(std::vector<int>::iterator i = MuTriggerIndicies.begin(); i != MuTriggerIndicies.end() && !passHLT; i++)
		passHLT |= (*viSelTriggersMu_singleLepCalc)[*i];
	}


	if(!passHLT) continue;

	Counters.increment("post HLT");

	//cout<< run_CommonCalc <<"  "<< lumi_CommonCalc<<"  "<< event_CommonCalc<<endl;//xxx
	//if(isIn1LepDataset( run_CommonCalc , lumi_CommonCalc, event_CommonCalc))
	//	cout<<"Single Lep Dataset Event Found: "<< run_CommonCalc <<"  "<< lumi_CommonCalc<<"  "<< event_CommonCalc<<endl;//xxx
		
	if(printlevel >= 6) std::cout << "incremented some coutners" << std::endl;	
	//make vectors of objects like electrons, muons, jets, and b-jets

	/////////////////////////////////////////////////////////////////
	//////////////////////// FIND THE PRIMARY VERTEX ////////////////
	/////////////////////////////////////////////////////////////////
	//This is already done by LJmet

	/////////////////////////////////////////////////////////////////
	//////////////////////// FIND MUONS  ////////////////////////////
	/////////////////////////////////////////////////////////////////
	//    __  __                       
	//   |  \/  |_   _  ___  _ __  ___ 
	//   | |\/| | | | |/ _ \| '_ \/ __|
	//   | |  | | |_| | (_) | | | \__ \
	//   |_|  |_|\__,_|\___/|_| |_|___/
	//                                 
	//xxmu
	if(printlevel >= 5) cout<<"Begind muons"<<endl;
	float MuT = 0;
	TLorentzVector VLep(0,0,0,0);
	/*if(nLooseMuons_CommonCalc != (int)muEnergy_singleLepCalc->size() ){
	  cout<<"neither tight nor nLooseMuons_CommonCalc matches muon vector size "<<endl;//soemtimes
	  else cout<<"muon vector size matches tight"<<endl; //this happens a lot
	  }*/

	vector<int> selMuonIndicies;
	int nSelMuons = 0;
	vector<int> TightMuonIndicies;
	int nTightMuons = 0;
	vector<int> LooseMuonIndicies;
	int nLooseMuons = 0;
	//vector<int> LooseNotTightMuonIndicies;
	//int nLooseNotTightMuons = 0;
	//for(int imu = 0; imu < nLooseMuons_CommonCalc; ++imu)
	for(int imu = 0; imu < (int)muEnergy_singleLepCalc->size(); ++imu){
	    if((*muPt_singleLepCalc)[imu] < 20.) continue;

	    if((*muIsLoose_singleLepCalc)[imu] and (*muMiniIso_singleLepCalc)[imu] < 0.4){
		nLooseMuons++;
		LooseMuonIndicies.push_back(imu);
		if((*muIsTight_singleLepCalc)[imu] and (*muMiniIso_singleLepCalc)[imu] < 0.2){
		
		    nTightMuons++;
		    TightMuonIndicies.push_back(imu);

		    TLorentzVector Vmu;
		    Vmu.SetPtEtaPhiE ( (*muPt_singleLepCalc)[imu], (*muEta_singleLepCalc)[imu], (*muPhi_singleLepCalc)[imu], (*muEnergy_singleLepCalc)[imu]); 
		    VLep = VLep + Vmu;
		}//end mu is tight
		//else{
		    //nLooseNotTightMuons++;
		    //LooseNotTightMuonIndicies.push_back(imu); //xxx why is this commented out??? Should be using this for bkgs. 
		//}
	    }//end mu is loose

	}//end muon for loop
	selMuonIndicies = LooseMuonIndicies;
	//selMuonIndicies = makeBkgs?LooseMuonIndicies:TightMuonIndicies;
	nSelMuons = nLooseMuons;
	//nSelMuons = makeBkgs?nLooseMuons:nTightMuons;
	indexPtSort(selMuonIndicies,nSelMuons, (*muPt_singleLepCalc));

	for(int imu = 0; imu<nSelMuons; ++imu) MuT+= (*muPt_singleLepCalc)[selMuonIndicies[imu]];




	/////////////////////////////////////////////////////////////////
	//////////////////////// FIND pfElectrons  //////////////////////
	/////////////////////////////////////////////////////////////////
	//    _____ _           _                       
	//   | ____| | ___  ___| |_ _ __ ___  _ __  ___ 
	//   |  _| | |/ _ \/ __| __| '__/ _ \| '_ \/ __|
	//   | |___| |  __/ (__| |_| | | (_) | | | \__ \
	//   |_____|_|\___|\___|\__|_|  \___/|_| |_|___/
	//                                              
	//xxele
	//if(nSelElectrons_CommonCalc != (int)elEnergy_singleLepCalc->size() )
	//cout<<"nSelElectrons_CommonCalc does not match the size of the el vectors"<<endl;

	if(printlevel >= 5) cout<<"Begind electrons"<<endl;
	vector<int> selEleIndicies;
	int nSelEle = 0;
	vector<int> LooseEleIndicies;
	int nLooseEle = 0;
	vector<int> TightEleIndicies;
	int nTightEle = 0;
	//vector<int> LooseNotTightEleIndicies;
	//int nLooseNotTightEle = 0;
	int nVetoEle= 0;
	float EleT = 0;
	if(printlevel >= 6) cout<<"Begind electrons loop"<<endl;
	for(int iele = 0; iele<(int)elEnergy_singleLepCalc->size(); ++iele){
	    if(printlevel >= 6) cout<<"Start ele loop for ele "<<iele<<endl;

	    if((*elPt_singleLepCalc)[iele] < 25.) continue; //online uct is 23.
	    if(printlevel >= 7) cout<<"ele survived pt cut"<<endl;

	    //dR 0.3 cut between electrons and muons with preference to muons 
	    bool eleFailsMuDRcut = false;
	    TLorentzVector Vele, Vmu;
	    Vele.SetPtEtaPhiE ((*elPt_singleLepCalc)[iele],(*elEta_singleLepCalc)[iele], (*elPhi_singleLepCalc)[iele], (*elEnergy_singleLepCalc)[iele]);
	    if(printlevel >= 7) cout<<"approach muon dR cutwork"<<endl;
	    for(int imu = 0; imu<nSelMuons && !eleFailsMuDRcut; ++imu){
		Vmu.SetPtEtaPhiE ( (*muPt_singleLepCalc)[selMuonIndicies[imu]], (*muEta_singleLepCalc)[selMuonIndicies[imu]], (*muPhi_singleLepCalc)[selMuonIndicies[imu]], (*muEnergy_singleLepCalc)[selMuonIndicies[imu]]); 
		eleFailsMuDRcut |= isSameObject(Vele, Vmu, 0.1);
	    } //end dR cut. 
	    if(eleFailsMuDRcut) continue;
	    if(printlevel >= 7) cout<<"survived muon dR cut"<<endl;

	    //int eleIdLevel = T50ns_F25ns?(*elIDLevel50ns_singleLepCalc)[iele]:(*elIDLevel25ns_singleLepCalc)[iele]; //eleIdLevel is generated by ele_ID_level_2015 in EleIDCuts.h
	    int eleIdLevel = 0;
	    switch(EleSelectionToUse){
		case cutbased:
		    eleIdLevel = ele_ID_level_2015(
			    T50ns_F25ns,
			    (*elEta_singleLepCalc)[iele],
			    (*elSihih_singleLepCalc)[iele],
			    (*elDeta_singleLepCalc)[iele],
			    (*elDphi_singleLepCalc)[iele],
			    (*elHoE_singleLepCalc)[iele],
			    (*elRelIso_singleLepCalc)[iele], //I'm guessing this is reliso with effective area. need to check. 
			    (*elOoemoop_singleLepCalc)[iele],
			    (*elD0_singleLepCalc)[iele],
			    (*elDZ_singleLepCalc)[iele],
			    (*elMHits_singleLepCalc)[iele],
			    (*elVtxFitConv_singleLepCalc)[iele]);
		    break;
		case MVA_nontrigger:
		    eleIdLevel = ele_ID_level_2015_MVA_nontrig( 
			    (*elPt_singleLepCalc)[iele],
			    (*elEta_singleLepCalc)[iele],
			    (*elMVAValue_singleLepCalc)[iele],
			    (*elMiniIso_singleLepCalc)[iele]);
		    break;
		case MVA_trigger:
		    cout<<"Error! Electron MVA selection is not yet possible"<<endl;
		    std::terminate();
		    /*eleIdLevel = ele_ID_level_2015_MVA_trig(
		      (*elMVAValue_singleLepCalc)[iele],
		      T50ns_F25ns,
		      (*elPt_singleLepCalc)[iele],
		      (*elEta_singleLepCalc)[iele],
		      (*elSihih_singleLepCalc)[iele],
		      (*elDeta_singleLepCalc)[iele],
		      (*elDphi_singleLepCalc)[iele], 
		      (*elHoE_singleLepCalc)[iele],
		      float EcalRelIso, //yet to be included in ntuple
		      float HcalRelIso,//yet to be included in ntuple
		      float Track03RelIso); //yet to be included in ntuple
		      break; */
	    }//end switch

	    //elID does not exist in data -> segfault.
	    if(printlevel >= 8) cout<<"got eleIDLevel: "<<eleIdLevel<<endl;
	    //ele ID level single lept calc. 

	    if(eleIdLevel >= 1) { //if veto or tighter
		nVetoEle++;
		if(eleIdLevel >= 2){ // if loose
		    LooseEleIndicies.push_back(iele);
		    ++nLooseEle;
		    //EleT += (*elPt_singleLepCalc)[iele];//with the tt cut, same as if you had asked for it with tight. 
		    if(eleIdLevel >= 4){ //if tight
			TightEleIndicies.push_back(iele);
			nTightEle++;
			if(printlevel >= 9) cout<<"try to get elPt"<<endl;
			//EleT += (*elPt_singleLepCalc)[iele];
			VLep = VLep + Vele;
			if(printlevel >= 9) cout<<"gotit"<<endl;
		    }
		    //LooseNotTightEleIndicies.push_back(iele);
		    //++nLooseNotTightEle;
		}//end if loose ele;
		else {
		    if(printlevel >= 8) cout<<"running continue"<<endl;
		    continue;
		}
	    }//end if veto electron
	    if(printlevel >= 6) cout<<"at end of ele loop"<<iele<<endl;

	}//end electron loop
	selEleIndicies = LooseEleIndicies;
	//selEleIndicies = makeBkgs?LooseEleIndicies:TightEleIndicies;
	nSelEle = nLooseEle;
	//nSelEle = makeBkgs?nLooseEle:nTightEle;
	indexPtSort(selEleIndicies,nSelEle,  (*elPt_singleLepCalc));

	for(int iele = 0; iele<nSelEle; ++iele) EleT+= (*elPt_singleLepCalc)[selEleIndicies[iele]];

	/////////////////////////////////////////////////////////////////
	//////////////////////// Lepton Cut /////////////////////////////
	/////////////////////////////////////////////////////////////////

	if(printlevel >= 5) cout<<"Make lepton Cut"<<endl;
	//then make a dielectron cut.
	int nLooseLep = nLooseMuons + nLooseEle;
	int nTightLep = nTightMuons + nTightEle;

	if(!( nLooseLep == 1 and nTightLep<=1)) continue;

	//if(nLooseEle + nLooseMuons <2 ) continue; //clean up all of this. xxx
	//Counters.increment("Has 2 leptons");

//	if(nVetoEle + nLooseMuons >2 ) continue;
	Counters.increment("Has exactly 1 lepton");

	if((nLooseLep >= 2)) cerr<<"Error! somehow the mono lepton cut is failing, multi-leptons!"<<endl;
	if((nLooseLep <= 0)) cerr<<"Error! somehow the mono lepton cut is failing, no leptons!"<<endl;
	assert(nTightLep < 2);
	assert(nLooseLep < 2);
	assert(nLooseLep > 0);

	float LepT = EleT + MuT;

	bool is_el = (bool) nSelEle;
	bool is_mu = (bool) nSelMuons;
	assert(is_el ^ is_mu);

	//if( (is_elel && is_mumu) || (is_mumu && is_elmu) || (is_elmu && is_elel)) 
	//cout<<"What the hell? Check the trilepton cut."<<endl;
	//assert(!((is_elel && is_mumu) || (is_mumu && is_elmu) || (is_elmu && is_elel)) );

	TLorentzVector Lep;
	if(is_el){
	    int a = selEleIndicies[0];
	    Lep.SetPtEtaPhiE((*elPt_singleLepCalc)[a],(*elEta_singleLepCalc)[a], (*elPhi_singleLepCalc)[a], (*elEnergy_singleLepCalc)[a]);
	}
	else if(is_mu){
	    int a = selMuonIndicies[0];
	    Lep.SetPtEtaPhiE((*muPt_singleLepCalc)[a],(*muEta_singleLepCalc)[a], (*muPhi_singleLepCalc)[a], (*muEnergy_singleLepCalc)[a]);
	}

	//cut out DY MC overlap
	/*	if(dmcblock->type == 71 && //if on DY
		((dmcblock->processID == 1 && Mll > 100) || //if DYJetsToLL50toINFf and M > 100
		(dmcblock->processID == 2 && Mll < 100))  ) //or if DYJetsToLL100to200f and M<100
		continue;*/

	/*if(isSpecial( run_CommonCalc , lumi_CommonCalc, event_CommonCalc) ){
		cout<<"Special Event Found after dilepton cut: "<< run_CommonCalc <<"  "<< lumi_CommonCalc<<"  "<< event_CommonCalc<<endl;//xxx
		if(is_elel) cout<<"double electron event"<<endl;
		if(is_elmu) cout<<"electron muon event"<<endl;
		if(is_mumu) cout<<"double muon event"<<endl;
	}*/
	nPassDilepCut++;
	if(!dmcblock->isMC){
	    if(   RegisterEvent(run_CommonCalc , lumi_CommonCalc, event_CommonCalc, EventRegistry)  ){
		//If this passes, this event is a duplicate. 
		//cout<<"Duplicate Event Found after dilepton cut: "<< run_CommonCalc <<"  "<< lumi_CommonCalc<<"  "<< event_CommonCalc<<endl;//xxx
		nDuplicates++;
		continue;
		}
	    if(   CheckEvent(run_CommonCalc , lumi_CommonCalc, event_CommonCalc, BadEventRegistry)  ){
                //n_rejected_BadEvents++;
                continue;
                }
	}//end if isData

	/////////////////////////////////////////////////////////////////
	//////////////////////// FIND JETS  /////////////////////////////
	/////////////////////////////////////////////////////////////////
	//        _      _       
	//       | | ___| |_ ___ 
	//    _  | |/ _ \ __/ __|
	//   | |_| |  __/ |_\__ \
	//    \___/ \___|\__|___/
	//                       


	bool there_is_an_opposing_jet = false;

	if(printlevel >= 5) cout<<"Begind jets"<<endl;
	vector<int> JetAK4Indicies;
	int nJetAK4= 0;
	//float HT = 0; //HT INVALID
	//TLorentzVector VHad(0,0,0,0);
	int closest_jet_index = -1;
	double minJetDR = 999.;
	for(int ijet = 0; ijet< (int)AK4JetPt_singleLepCalc->size(); ++ijet){ //jets are approxomately sorted in Pt order. 
	    //but you saw an instance of
	    //JetPt: 189.88 25.9926 80.4774 
	    //strange thing: we always have at least two jets. 
	    //these are born with pt > 25 cut and |eta| < 2.5 cut
	    if((*AK4JetPt_singleLepCalc)[ijet] < 30 || fabs((*AK4JetEta_singleLepCalc)[ijet]) > 2.4) continue;
	    //if((*AK4JetPt_singleLepCalc)[ijet] < 40 || fabs((*AK4JetEta_singleLepCalc)[ijet]) > 2.4) continue;


	    TLorentzVector vJ;
	    vJ.SetPtEtaPhiE((*AK4JetPt_singleLepCalc)[ijet],(*AK4JetEta_singleLepCalc)[ijet], (*AK4JetPhi_singleLepCalc)[ijet], (*AK4JetEnergy_singleLepCalc)[ijet]);
	    float dR_jet_lep = dR(vJ,Lep);
	    if(dR_jet_lep > 1.0) there_is_an_opposing_jet = true;

	    if(dR_jet_lep < minJetDR and dR_jet_lep < 0.4){
		minJetDR = dR_jet_lep;
		closest_jet_index = ijet;
	    }

	    //seperate the leptons from the jets using cuts listed in B2G-12-015
	    /*
 		//This is removed so that it's possible to look at jets overlapping your leptons. 
 		//The absence of this code invalidates your HT. See lines marked HT INVALID
	    bool jet_passes_DR = true;
	    TLorentzVector vJ, vEl, vMu;
	    vJ.SetPtEtaPhiE((*AK4JetPt_singleLepCalc)[ijet],(*AK4JetEta_singleLepCalc)[ijet], (*AK4JetPhi_singleLepCalc)[ijet], (*AK4JetEnergy_singleLepCalc)[ijet]);
	    for(int iel = 0; iel<nSelEle && jet_passes_DR; ++iel){
		int iele = selEleIndicies[iel];
		vEl.SetPtEtaPhiE((*elPt_singleLepCalc)[iele],(*elEta_singleLepCalc)[iele], (*elPhi_singleLepCalc)[iele], (*elEnergy_singleLepCalc)[iele]);
		jet_passes_DR &= !isSameObject(vJ, vEl, 0.3);
	    }	
	    if(!jet_passes_DR) continue;
	    for(int iMu = 0; iMu<nSelMuons && jet_passes_DR; ++iMu){
		int imu = selMuonIndicies[iMu];
		vMu.SetPtEtaPhiE((*muPt_singleLepCalc)[imu], (*muEta_singleLepCalc)[imu], (*muPhi_singleLepCalc)[imu], (*muEnergy_singleLepCalc)[imu]); 
		jet_passes_DR &= !isSameObject(vJ, vMu, 0.4);
	    }	
	    if(!jet_passes_DR) continue;
		*/

	    //jet passes:
	    JetAK4Indicies.push_back(ijet);
	    nJetAK4++;
	    //HT += (*AK4JetPt_singleLepCalc)[ijet]; //HT INVALID
	    //VHad = VHad + vJ;
	}//end for all jets
	indexPtSort(JetAK4Indicies, nJetAK4, (*AK4JetPt_singleLepCalc));
	//
	//      AK4HT_singleLepCalc
	//      theJetHT_JetSubCalc

	float coneCorPt = Lep.Pt();
	if(closest_jet_index > -1){
		//FIX
		//elMiniIso_singleLepCalc
		//muMiniIso_singleLepCalc
	}




	//float ST = LepT + HT + corr_met_singleLepCalc;  //HT INVALID
	//TLorentzVector VSum = VHad + VLep;
	//float MSum = VSum.M();
	//float MtSum = Mt(VSum, corr_met_singleLepCalc, corr_met_phi_singleLepCalc);

	/*float leadjetDR = -1;
	if(nJetAK4 >=2){
	    TLorentzVector vJ1, vJ2;
	    int ijet = JetAK4Indicies[0];
	    vJ1.SetPtEtaPhiE((*AK4JetPt_singleLepCalc)[ijet],(*AK4JetEta_singleLepCalc)[ijet], (*AK4JetPhi_singleLepCalc)[ijet], (*AK4JetEnergy_singleLepCalc)[ijet]);
	    ijet = JetAK4Indicies[1];
	    vJ2.SetPtEtaPhiE((*AK4JetPt_singleLepCalc)[ijet],(*AK4JetEta_singleLepCalc)[ijet], (*AK4JetPhi_singleLepCalc)[ijet], (*AK4JetEnergy_singleLepCalc)[ijet]);
	    leadjetDR = dR(vJ1,vJ2);
	}*///end calculate leadjetDR

	/////////////////////////////////////////////////////////////////
	//////////////////////// FIND B-JETS  ///////////////////////////
	/////////////////////////////////////////////////////////////////
	//    ____            _      _       
	//   | __ )          | | ___| |_ ___ 
	//   |  _ \ _____ _  | |/ _ \ __/ __|
	//   | |_) |_____| |_| |  __/ |_\__ \
	//   |____/       \___/ \___|\__|___/
	//                                   
	// xxbtag
	float BTl = 0;
	float BTm = 0;
	//theJetBTag_JetSubCalc
	//theJetFlav_JetSubCalc;
	//theJetCSV_JetSubCalc;
	//int AK4JetBTag_singleLepCalc;  = int value representing a bool = isJetTagged
	//double AK4JetBDisc_singleLepCalc; = pfCombinedInclusiveSecondaryVertexV2BJetTags
	//CISVv2L
	//CISVv2M
	//CISVv2T 

	if(printlevel >= 5) cout<<"Begind Bjets"<<endl;
	vector<int> BJetAK4CISVmIndicies; //these include all the efficiencies and scale factors
	int nBJetAK4CISVm= 0;

	vector<int> BJetAK4CISVtIndicies; //dumb CISV cut, no efficiencies or SF's
	int nBJetAK4CISVt= 0;

	vector<int> BJetAK4CISVmSimpleIndicies;//dumb CISV cut, no efficiencies or SF's
	int nBJetAK4CISVmSimple= 0;

	vector<int> BJetAK4CISVlIndicies;//dumb CISV cut, no efficiencies or SF's
	int nBJetAK4CISVl= 0;

	for(int ijet = 0; ijet<nJetAK4; ++ijet){
	    if((*AK4JetBTag_singleLepCalc)[JetAK4Indicies[ijet]] == 1){ //these include all the scale factors and efficiencies.
		BJetAK4CISVmIndicies.push_back(JetAK4Indicies[ijet]);
		nBJetAK4CISVm++;
	    }
	    float cisv = (*AK4JetBDisc_singleLepCalc)[JetAK4Indicies[ijet]];


	    if(cisv > CISVv2L){
		BJetAK4CISVlIndicies.push_back(JetAK4Indicies[ijet]);
		nBJetAK4CISVl++;
		BTl += (*AK4JetPt_singleLepCalc)[JetAK4Indicies[ijet]];

		if(cisv > CISVv2M){
		    BJetAK4CISVmSimpleIndicies.push_back(JetAK4Indicies[ijet]);
		    nBJetAK4CISVmSimple++;
		    BTm += (*AK4JetPt_singleLepCalc)[JetAK4Indicies[ijet]];

		    if(cisv > CISVv2T){
			BJetAK4CISVtIndicies.push_back(JetAK4Indicies[ijet]);
			nBJetAK4CISVt++;
		    }//end if tight

		}//end if medium
	    }//end if loose

	}//end for all used jets

	indexPtSort(BJetAK4CISVmIndicies, nBJetAK4CISVm, (*AK4JetPt_singleLepCalc));
	indexPtSort(BJetAK4CISVtIndicies, nBJetAK4CISVt, (*AK4JetPt_singleLepCalc));
	indexPtSort(BJetAK4CISVmSimpleIndicies, nBJetAK4CISVmSimple, (*AK4JetPt_singleLepCalc));
	indexPtSort(BJetAK4CISVlIndicies, nBJetAK4CISVl, (*AK4JetPt_singleLepCalc));

	/*float leadBjetDR = -1;
	if(nBJetAK4CISVm >=2){
	    TLorentzVector vJ1, vJ2;
	    int ijet = BJetAK4CISVmIndicies[0];
	    vJ1.SetPtEtaPhiE((*AK4JetPt_singleLepCalc)[ijet],(*AK4JetEta_singleLepCalc)[ijet], (*AK4JetPhi_singleLepCalc)[ijet], (*AK4JetEnergy_singleLepCalc)[ijet]);
	    ijet = BJetAK4CISVmIndicies[1];
	    vJ2.SetPtEtaPhiE((*AK4JetPt_singleLepCalc)[ijet],(*AK4JetEta_singleLepCalc)[ijet], (*AK4JetPhi_singleLepCalc)[ijet], (*AK4JetEnergy_singleLepCalc)[ijet]);
	    leadBjetDR = dR(vJ1,vJ2);
	}//end calculate leadjetDR
	*/

	if(printlevel >= 5) cout<<"Ready for Physics"<<endl;

	//                   __    __        __          __    __        __          __    __
	//                  (//    \\)    __(//   __    (//    \\)    __(//   __    (//    \\)  
	//                  /"      / __  \\)"    \\)_  /"      / __  \\)"    \\)_  /"      / __
	//                '|-..__..-''\_''-.\__..-''  '|-..__..-''\_''-.\__..-''  '|-..__..-''\
	//                (\\  \_    _(\\      _/     (\\  \_    _(\\      _/     (\\  \_    //)
	//                 ""  (\\  //)""     //)      ""  (\\  //)""     //)      ""  (\\   ""
	//                      ""  ""        ""            ""  ""        ""            ""
	//////////////////////////////////////////////////////////////////////
	//////////////////////// READY FOR PHYSICS ///////////////////////////
	//////////////////////////////////////////////////////////////////////
	//xxphysics
	
	//calculate SS
	/*bool is_SSDL = false; //else OSDL
	if(is_elel)      is_SSDL= (*elCharge_singleLepCalc)[selEleIndicies[0]] == (*elCharge_singleLepCalc)[selEleIndicies[1]];
	else if(is_mumu) is_SSDL= (*muCharge_singleLepCalc)[selMuonIndicies[0]] == (*muCharge_singleLepCalc)[selMuonIndicies[1]];
	else if(is_elmu) is_SSDL= (*elCharge_singleLepCalc)[selEleIndicies[0]] == (*muCharge_singleLepCalc)[selMuonIndicies[0]];
	else{
	    cerr<<"WTF dilepton pair that is neither ee, em, or mm."<<endl;
	    assert(0);
	}*/
	bool there_is_a_Mlj_on_Z = false;
	
        for(int ijet = 0; ijet < nJetAK4 and not there_is_a_Mlj_on_Z; ijet++){
            TLorentzVector vJ;
            vJ.SetPtEtaPhiE((*AK4JetPt_singleLepCalc)[ijet],(*AK4JetEta_singleLepCalc)[ijet], (*AK4JetPhi_singleLepCalc)[ijet], (*AK4JetEnergy_singleLepCalc)[ijet]);
            float Mlj = (vJ+Lep).M();
	    there_is_a_Mlj_on_Z |= (Mlj>71. and Mlj<91.);
        }//end calculate leadjetDR




	bool is_CR1 = (corr_met_singleLepCalc < 20.) and 
		(nJetAK4 >= 1) and
		there_is_an_opposing_jet and
		Mt(Lep) < 20.;

	bool is_CR_T53 = (corr_met_singleLepCalc < 25.) and
                (nJetAK4 >= 1) and
                there_is_an_opposing_jet and
		Mt(Lep) < 25. and
		!there_is_a_Mlj_on_Z;



    //void Fill_HLT_study(TH2F* h, std::map<int,int> indexmap, int nTrig, std::vector<int> TriggerIndicies, vector<int> *viSelTriggers){
    //
	Fill_HLT_study(HLT_coincidience_ele, HLT_map_ele, NEleTriggers, EleTriggerIndicies, viSelTriggersEl_singleLepCalc,false);
	Fill_HLT_study(HLT_coincidience_mu,  HLT_map_mu , NMuTriggers,  MuTriggerIndicies,  viSelTriggersMu_singleLepCalc,false);
	if(is_CR1){//FIX
	    Fill_HLT_study(HLT_coincidience_ele_CR1, HLT_map_ele, NEleTriggers, EleTriggerIndicies, viSelTriggersEl_singleLepCalc,false);
	    Fill_HLT_study(HLT_coincidience_mu_CR1,  HLT_map_mu , NMuTriggers,  MuTriggerIndicies,  viSelTriggersMu_singleLepCalc,false);
	}

	//////////////////////////////////////////////////////////////////////
	////////////////////////// Fill Variables ////////////////////////////
	//////////////////////////////////////////////////////////////////////

	float weight = (MCWeight_singleLepCalc<0.)? -1. : 1.;
	if(weight >0) nPos++;
	else nNeg++;


	if(is_CR1){
	    if(is_mu){
		h_tightloose_CR1_mu->Fill(nTightLep);
		if(nTightLep > 0){
		    h_lepEta_Tight_CR1_mu->Fill(Lep.Eta());
		    h_lepPt_Tight_CR1_mu->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Tight_CR1_mu->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Tight_CR1_mu->Fill(coneCorPt,Lep.Eta());
		}
		else{

		    h_lepEta_Loose_CR1_mu->Fill(Lep.Eta());
		    h_lepPt_Loose_CR1_mu->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Loose_CR1_mu->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Loose_CR1_mu->Fill(coneCorPt,Lep.Eta());
		}
	    }//end is muon
	    if(is_el){
		h_tightloose_CR1_ele->Fill(nTightLep); //h_tightloose_CR1_ele = new TH1F("h_tightloose_CR1","0 loose, 1 tight, 2 ratio",3,0,3);

		if(nTightLep > 0){
		    h_lepEta_Tight_CR1_ele->Fill(Lep.Eta());
		    h_lepPt_Tight_CR1_ele->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Tight_CR1_ele->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Tight_CR1_ele->Fill(coneCorPt,Lep.Eta());
		}
		else{

		    h_lepEta_Loose_CR1_ele->Fill(Lep.Eta());
		    h_lepPt_Loose_CR1_ele->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Loose_CR1_ele->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Loose_CR1_ele->Fill(coneCorPt,Lep.Eta());
		}
	    }//end is_ele
	}//end is control region 1

	if(is_CR_T53){
	    if(is_mu){
		h_tightloose_CR_T53_mu->Fill(nTightLep);
		if(nTightLep > 0){
		    h_lepEta_Tight_CR_T53_mu->Fill(Lep.Eta());
		    h_lepPt_Tight_CR_T53_mu->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Tight_CR_T53_mu->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Tight_CR_T53_mu->Fill(coneCorPt,Lep.Eta());
		}
		else{

		    h_lepEta_Loose_CR_T53_mu->Fill(Lep.Eta());
		    h_lepPt_Loose_CR_T53_mu->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Loose_CR_T53_mu->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Loose_CR_T53_mu->Fill(coneCorPt,Lep.Eta());
		}
	    }//end is muon
	    if(is_el){
		h_tightloose_CR_T53_ele->Fill(nTightLep); //h_tightloose_CR_T53_ele = new TH1F("h_tightloose_CR_T53","0 loose, 1 tight, 2 ratio",3,0,3);

		if(nTightLep > 0){
		    h_lepEta_Tight_CR_T53_ele->Fill(Lep.Eta());
		    h_lepPt_Tight_CR_T53_ele->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Tight_CR_T53_ele->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Tight_CR_T53_ele->Fill(coneCorPt,Lep.Eta());
		}
		else{

		    h_lepEta_Loose_CR_T53_ele->Fill(Lep.Eta());
		    h_lepPt_Loose_CR_T53_ele->Fill(Lep.Pt());
		    h_lepEta_vs_lepPt_Loose_CR_T53_ele->Fill(Lep.Pt(),Lep.Eta());
		    h_lepEta_vs_coneCorPt_Loose_CR_T53_ele->Fill(coneCorPt,Lep.Eta());
		}
	    }//end is_ele
	}//end is control region T53

	}//end main event loop
	//////////////////////////////////////////////////////////////
	//////////////////////// FINISH UP ///////////////////////////
	///////////////////////////////////////////////////////////////


	TFile* fout = new TFile(dmcblock->string_meta["EventLoopOutRoot"].c_str(),"RECREATE");
	fout->cd(); //write out the files



	/////////////////////Write Kinvar Histograms//////////////////
	//h_MinMlb_OSDL1sanslb->Write();

	h_tightloose_CR1_ele->Write();
	h_lepEta_Loose_CR1_ele->Write();
	h_lepEta_Tight_CR1_ele->Write();
	h_lepPt_Loose_CR1_ele->Write();
	h_lepPt_Tight_CR1_ele->Write();
	h_lepEta_vs_lepPt_Loose_CR1_ele->Write();
	h_lepEta_vs_lepPt_Tight_CR1_ele->Write();
	h_lepEta_vs_coneCorPt_Loose_CR1_ele->Write();
	h_lepEta_vs_coneCorPt_Tight_CR1_ele->Write();

	h_tightloose_CR1_mu->Write();
	h_lepEta_Loose_CR1_mu->Write();
	h_lepEta_Tight_CR1_mu->Write();
	h_lepPt_Loose_CR1_mu->Write();
	h_lepPt_Tight_CR1_mu->Write();
	h_lepEta_vs_lepPt_Loose_CR1_mu->Write();
	h_lepEta_vs_lepPt_Tight_CR1_mu->Write();
	h_lepEta_vs_coneCorPt_Loose_CR1_mu->Write();
	h_lepEta_vs_coneCorPt_Tight_CR1_mu->Write();

	HLT_coincidience_ele_CR1->Write();
	HLT_coincidience_mu_CR1->Write();
	HLT_coincidience_ele->Write();
	HLT_coincidience_mu->Write();

	h_tightloose_CR_T53_ele->Write();
	h_lepEta_Loose_CR_T53_ele->Write();
	h_lepEta_Tight_CR_T53_ele->Write();
	h_lepPt_Loose_CR_T53_ele->Write();
	h_lepPt_Tight_CR_T53_ele->Write();
	h_lepEta_vs_lepPt_Loose_CR_T53_ele->Write();
	h_lepEta_vs_lepPt_Tight_CR_T53_ele->Write();
	h_lepEta_vs_coneCorPt_Loose_CR_T53_ele->Write();
	h_lepEta_vs_coneCorPt_Tight_CR_T53_ele->Write();

	h_tightloose_CR_T53_mu->Write();
	h_lepEta_Loose_CR_T53_mu->Write();
	h_lepEta_Tight_CR_T53_mu->Write();
	h_lepPt_Loose_CR_T53_mu->Write();
	h_lepPt_Tight_CR_T53_mu->Write();
	h_lepEta_vs_lepPt_Loose_CR_T53_mu->Write();
	h_lepEta_vs_lepPt_Tight_CR_T53_mu->Write();
	h_lepEta_vs_coneCorPt_Loose_CR_T53_mu->Write();
	h_lepEta_vs_coneCorPt_Tight_CR_T53_mu->Write();

	/////////////////End Write Kinvar Histograms//////////////////

	fout->Close();


	if(!onfire){
	    ofstream firealarm2;
	    firealarm2.open("tmpfirealarm.txt");
	    firealarm2 << "SusyAna_MC_BTagEff_arg is ok"<<endl;
	    firealarm2.close();
	}
    }//end Loop

    bool fakeratecalc::isSameObject(TLorentzVector& p1, TLorentzVector& p2, const float dRmax){
	float dEta = p1.Eta() - p2.Eta();
	float dPhi = TVector2::Phi_mpi_pi(p1.Phi() - p2.Phi());
	return std::sqrt(dEta*dEta + dPhi*dPhi) < dRmax;
    }//dR cut

//       ____________________________________________________________
//	   /\    /\    /\    /\    /\    /\    /\    /\    /\    /\
//	  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \
//	 /    \/    \/    \/    \/    \/    \/    \/    \/    \/    \
//	 \    /\    /\    /\    /\    /\    /\    /\    /\    /\    /
//	  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /  \  /
//	 __\/____\/____\/____\/____\/____\/____\/____\/____\/____\/__
//	 


    TH1F* init_TH1F_special(string label, string var, LabelKinVars &KinVars ){ //cant throw  <std::string,int> errorpair
	KinematicVar* k= KinVars->get_throwable(var,3);
	return init_TH1F_special(label, k); 
    }
    TH1F* init_TH1F_special(string label, KinematicVar* k){ 
	TH1F* h;
	if(k->useCustomBinning) h = new TH1F( label.c_str(),
		(char*)(k->titles+";"+k->xlabels).c_str(),
		k->nbins, k->CustomBinning);

	else h = new TH1F( label.c_str(),
		(char*)(k->titles+";"+k->xlabels).c_str(),
		k->nbins, k->hmin, k->hmax);
	return h;
    }

    TH2F* init_TH2F_special(string label, KinematicVar* kx, KinematicVar* ky){ 
	TH2F* h;
	if(kx->useCustomBinning){
	    if(ky->useCustomBinning){
		h = new TH2F( label.c_str(),
			(char*)(kx->titles+";"+kx->xlabels+";"+ky->xlabels).c_str(),
			kx->nbins, kx->CustomBinning,
			ky->nbins, ky->CustomBinning);
	    }
	    else{
		h = new TH2F( label.c_str(),
			(char*)(kx->titles+";"+kx->xlabels+";"+ky->xlabels).c_str(),
			kx->nbins, kx->CustomBinning,
			ky->nbins, ky->hmin, ky->hmax);
	    }
	}
	else{
	    if(ky->useCustomBinning){
		h = new TH2F( label.c_str(),
			(char*)(kx->titles+";"+kx->xlabels+";"+ky->xlabels).c_str(),
			kx->nbins, kx->hmin, kx->hmax,
			ky->nbins, ky->CustomBinning);
	    }
	    else{
		h = new TH2F( label.c_str(),
			(char*)(kx->titles+";"+kx->xlabels+";"+ky->xlabels).c_str(),
			kx->nbins, kx->hmin, kx->hmax,
			ky->nbins, ky->hmin, ky->hmax);
	    }
	}
	return h;
    }//end init_TH2F_special

    std::map<int,int> Make_HLT_map_selIndex_to_stringIndex(int NTriggers, string Triggers[], vector<string>  * vsSelTriggers){
	std::map<int,int> m;
	int i= 0;
	//cout<<"in Make_HLT_map_selIndex_to_stringIndex"<<endl;
	for( std::vector<string>::iterator itrig = vsSelTriggers->begin(); itrig != vsSelTriggers->end(); itrig++){ 
		//cout<<"consider vsSelTriggers::"<<*itrig<<endl;
	    for(int jtrig = 0; jtrig<NTriggers; jtrig++){
		//cout<<"  against "<<Triggers[jtrig]<<endl;
		if((*itrig).compare(0,Triggers[jtrig].size(),Triggers[jtrig]) ==0  ){ //if the vector on the lists starts with Trigger[jtrig]
		    m[i] = jtrig; //map to the first string named. 
			//cout<<"Mapping LJMET trigger index "<<i<<" to string index "<<jtrig<<endl;
		    break;
		}
	    }//end for every trigger that you actually want
	    i++;
	}//for every trigger in the tuple's vector
	return m;
    }

    void MakeTriggerIndicies(std::vector<int> & TriggerIndicies, int NTriggers, string Triggers[], vector<string>  * vsSelTriggers){
	int i= 0;
	for( std::vector<string>::iterator itrig = vsSelTriggers->begin(); itrig != vsSelTriggers->end(); itrig++){ 
	    for(int jtrig = 0; jtrig<NTriggers; jtrig++){
		if((*itrig).compare(0,Triggers[jtrig].size(),Triggers[jtrig]) ==0  ){ //if the vector on the lists starts with Trigger[jtrig]
		    TriggerIndicies.push_back(i);
		    break;
		}
	    }//end for every trigger that you actually want
	    i++;
	}//for every trigger in the tuple's vector
    }
    bool checkTriggerList(std::vector<int> & TriggerIndicies, int NTriggers, string Triggers[], vector<string>  * vsSelTriggers){
	bool all_trigs_are_ok = true;
	for(std::vector<int>::iterator i=TriggerIndicies.begin(); i!=TriggerIndicies.end();i++){ //for every index that ought to be there.
	    bool trig_is_ok = false;
	    for(int jtrig = 0; jtrig<NTriggers && !trig_is_ok; jtrig++){
		trig_is_ok |= (*vsSelTriggers)[*i].compare(0,Triggers[jtrig].size(),Triggers[jtrig]) ==0;
	    }
	    all_trigs_are_ok &= trig_is_ok;
	    if(!trig_is_ok) cout<<"Warning! trigger at index "<<*i<< " has moved and is not on the trigger list. The trig at that index now reads "<<(*vsSelTriggers)[*i]<<endl;
	}
	return all_trigs_are_ok;
    }


    void indexPtSort(vector<int>& indicies, int n, vector<double> pt){
	//given an array length n of indecies of pt, sort the elements of indicies[] so that pt[indicies[]] is in decending order.

	if(n<2) return;
	else if(n==2){
	    if(pt[indicies[0]]<pt[indicies[1]]){ 
		//swap the indicies so the higher pt one is first
		//note: the xor swap trick works here, but would take 70% longer to run.
		int temp = indicies[0];
		indicies[0] = indicies[1];
		indicies[1] = temp;
	    }
	}//end if 2
	else{
	    //check if sorted
	    bool is_sorting = true;
	    double last = n>0?pt[indicies[0]]:0;
	    for(int i= 1; i<n; ++i){
		double now = pt[indicies[i]];
		if(now>last){ 
		    is_sorting = false;
		    break;
		}
		last = now;
	    }
	    if(is_sorting) return;

	    //bundle the index and the pt
	    std::vector<mypair*> pairs;
	    for(int i= 0; i<n; ++i){
		mypair *apair = new mypair();
		apair->index = indicies[i];
		apair->pt = pt[indicies[i]];
		pairs.push_back(apair);
	    }
	    std::sort(pairs.begin(), pairs.end(), pairCompare);
	    for(int i= 0; i<n; ++i) indicies[i] = pairs[i]->index;//put the reults back.
	}
    }//end sort by pt

    bool pairCompare(mypair* firstElem, mypair* secondElem) {
	return firstElem->pt > secondElem->pt;
    }

    void Fill_HLT_study(TH2F* h, std::map<int,int> indexmap, int nTrig, std::vector<int> TriggerIndicies, vector<int> *viSelTriggers, bool verbose){
	for(std::vector<int>::iterator i = TriggerIndicies.begin(); i != TriggerIndicies.end(); i++){
	    if(!(*viSelTriggers)[*i]) continue;
	    int I = indexmap[*i];
	    h->Fill(I,I); //diag is how many passed this trigger in general. 
		if(verbose) cout<<"Filling diagonal at "<<I<<endl;
	    bool something_else_passed = false;
	    for(std::vector<int>::iterator j = TriggerIndicies.begin(); j != TriggerIndicies.end(); j++){
		if(*j == *i) continue;
		if(!(*viSelTriggers)[*j]) continue;
		something_else_passed = true;
		int J = indexmap[*j];
		h->Fill(I,J);
		if(verbose) cout<<"Filling off- diagonal at "<<I<<" "<<J<<endl;
	    }//end second loop;
	    if(!something_else_passed){
		h->Fill(nTrig,I);
		if(verbose) cout<<"Filling unique at "<<I<<endl;
	    }
	}//end first loop
    }//end Fill_HLT_study

