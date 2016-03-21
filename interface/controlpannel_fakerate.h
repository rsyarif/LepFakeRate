#ifndef _CONTROLPANNEL_H // if not defined
#define _CONTROLPANNEL_H 

#include <iostream>
#include <vector>
using namespace std; 

namespace controlpannel{

//data
const float Integrated_Luminosity_Data = 2.215;//0.552;//in inverse femptobarns
/*
string plotsroot= "raw_plots_data_GLL12ABC.root"; 
string logplotsroot= "raw_plots_data_GLL12ABC.log";
string formatedplotsroot= "formatted_plots_dataGLL12.root";
bool useMETFilter_data = true;
*/

//MC
//int type_to_run = 10;
//string which_MC_to_use2 = "FSst_185_mu_150"; //this should become "which_DMCblock_to_use"
const bool makeEventsList = false;
const bool T50ns_F25ns = false;
const bool fix_BR_event_by_event = true; //turn this off when doing the real thing. 

enum EleSelection {cutbased, MVA_nontrigger, MVA_trigger};
EleSelection EleSelectionToUse = MVA_nontrigger;

//Data+MC
//string formatedplotsroot_combined= "formatted_plots_dataMC.root"; //output from format_plots.C

//Limit Setting
//bool remake_Limit_packages = true;

long int probeevent = 10496;

//shift jet energy corrections and all up/down by one sigma. 
int bumpJEC = 0;//is allowed to be +- 1 and 0	//assuems 2% uncertainty.  
int bumpBtagEff = 0;//is allowed to be +- 1 and 0 

	///LISTS FOR INDEXING
const int nEventTopologies = 36; // the number of types of cuts selected, like 1JB...
string s_EventTopology[nEventTopologies] = {
 "DL","OSDL1","OSDL2","SSDL",//4
 "OSDLraw", "SSDLraw",//2
 "OSDLm", "OSDLMll20", "OSDLmMll20",//3
 "SSDLm", "SSDLMll20", "SSDLbasemMll20",//3
 "OSDL1mMll", "OSDL1basemMllJ", "OSDL1basemMllJB",//3
 "OSDL2mMllJ4", "OSDL2mMllJ4B", "OSDL2basemMllJ5", "OSDL2basemMllJ5B",//4
 "SSDLmMllJ2", "SSDLmMllJ2B", "SSDLbasemMllJ", "SSDLbasemMllJB",//4
 "OSDL1sansm", "OSDL1sansBlb", "OSDL1sansMll", "OSDL1sansmlb", "OSDL1sansMlllb", "OSDL1sanslb",//6
 "OSDL2sansm", "OSDL2sansMll", "OSDL2sansJ5",//3
 "SSDLsansm", "SSDLsansMll", "SSDLsansJ3", "SSDLsansB"};//4




const int nKinemVars_all = 25;
string s_KinemVars_all[nKinemVars_all]={
				"nBl", "nJ",
                                "nBm", "BTm", "BTl",
                                "LHT", "HT", "ST",
                                "LepT", //"EleTT", "MuT",
                                "MET", "Mll", "MinMlb",
				"lepPt","lepEta","lepDPhi", "lepDR",
				"jetPt","jetEta","leadjetDR",
				"leadBjetDR", "lepJetDR",
                                //"yield", "yieldsum"};
				//"BLovHT", "BMovHT", "LepTovHT", 
				//"HTovST", "LepTovST", "METovST", "MSumovST", 
				"MSum", "MtSum",
				"tightloose", "coneCorrPt"};
				//nJ
				//Ask for JetPtThreshl

//const int nEventTopologies_limit = 5; // the number of types of cuts selected, like 1JB...
//string s_EventTopology_limit[nEventTopologies_limit] = {
	//"NULL","gbar2",
	//"2JewkMjjgbar2lep0notH"};

		//Will be used if either bumpBtagEff != 0 or bumpJEC != 0. 
		//list of kin vars to make limits on. These must be singly-filled per event and otherwise make sense.
//const int nKinemVars_limit = 2;
//string s_KinemVars_limit[nKinemVars_limit] = {
//"MET", "PtGGPtBB"};//, "PtGGPtVjj",


bool IsTopo(string topo);
bool IsKinvar(string kinvar);
//bool IsTopoLim(string topo);
//bool IsKinvarLim(string kinvar);

bool IsTopo(string topo){
                //      true iff topo is in s_EventTopology.
        for (int iTopo = 0; iTopo<nEventTopologies; iTopo++)
                if( topo.compare(s_EventTopology[iTopo]) == 0)
		    return true;
        return false;
} //end IsTopo
bool IsKinvar(string kinvar){
                //      true iff kinvar is in s_KinemVars_all.
        for (int kKinVar = 0; kKinVar<nKinemVars_all; kKinVar++)
                if( kinvar.compare(s_KinemVars_all[kKinVar]) == 0) //I think this is a "contains and starts at 0, not equality"
		    return true;
        return false;
}// end IsKinvar
/*bool IsTopoLim(TString topo){
                //      true iff topo is in s_EventTopology_limit.
        for (int iTopo = 0; iTopo<nEventTopologies_limit; iTopo++)
                if( topo.compare(s_EventTopology_limit[iTopo]) == 0)
		    return true;
        return false;
}*/ //end IsTopoLim
/*bool IsKinvarLim(TString kinvar){
                //      true iff kinvar is in s_KinemVars_limit.
        for (int kKinVar = 0; kKinVar<nKinemVars_limit; kKinVar++)
                if( kinvar.compare(s_KinemVars_limit[kKinVar]) == 0) //I think this is a "contains and starts at 0, not equality"
                                return true;
        return false;
}//end IsKinvarLim
*/
	const float CISVv2L = 0.605;
	const float CISVv2M = 0.890;
	const float CISVv2T = 0.907;

	//Controls for post.C
	const bool showData = true; 
	const bool preliminary = true;
	const bool saveImages = true;
        const bool writeRootFile = true;
        const string savewhat = "gpdf";//ggifpdfeps
        const float SignalInflationFactor = 1000; //scale up the signal by this much. 
        const bool divideOutBinWidths = 1;

	string plotsdirC   = "plots/C/";
	string plotsdirroot ="plots/root/";
	string plotsdirjpg = "plots/jpg/";
	string plotsdirpdf = "plots/pdf/";
	string plotsdirpng = "/uscms_data/d2/abarker/public_html/plots/";
	string plotsdirgif = "plots/gif/";
	string plotsdireps = "plots/eps/";
	string htmldir   = "/uscms_data/d2/abarker/public_html/";

	const int nmodes = 7;
	static const int nBByields = 6;
	static const int nSysYields = 1+//nominal
					2+//PU corrections
					2+//renorm envelope
					2+//pdf 1 sigmas. 
					2;//Julie Jet SFs 
	static const int nRenorm = 6;
	static const int nPDF = 100;
	static const int nPolySysYields = 1+//nominal
					2+//PU corrections
					nRenorm+//renorm envelope
					nPDF+//pdf 1 sigmas. 
					2;//Julie Jet SFs
	//mode 0 = nominal. 1 = fakerate plus, 2 = fakerate minus, 3 = passrate plus, 4 = passrate minus 5=qg fake plus, 6=qg fake minus.
	//see feModeBehavior in fakerate.h
	bool useFakeQflips = true;

} //end namespace params
#endif

