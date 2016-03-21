#define runone_cxx
#include "fakeratecalc.C"
#include "interface/DMCblock_fakerates.h"
#include "interface/DMCblock_list_fakerates.h"
#include "interface/DMCclass_list_fakerates.h"
#include <iostream>
#include <time.h>
#include <TH1F.h>
#include "eventRegistry.C"
#include "populate_BadEventRegistry.C"
using namespace std;
//nohup root -l -b -q runall.C+ > runall.log
void runonefake(){
	clock_t startfor = clock();
	TH1::AddDirectory(kFALSE);
	//DMCclass * dc =setupDMCclass("DYJetsToLLf");
	//DMCclass * dc =setupDMCclass("ttjf");
	DMCclass * dc =setupDMCclass("Data2lepf");
	//DMCclass * dc =setupDMCclass("Data1lepf");
	//DMCclass * dc =setupDMCclass("TpTp800f");
	//DMCclass * dc =setupDMCclass("QCDf");
	eventRegistry EventRegistry; //for data only, registers all events
	eventRegistry BadEventRegistry; //for data only, a list of known bad events to reject. 
	bool there_is_data = false;
	for(std::vector<DMCblock*>::iterator dit = dc->blocks.begin(); dit != dc->blocks.end(); dit++) there_is_data |= (*dit)->type <0;

	if(there_is_data) populate_BadEventRegistry(&BadEventRegistry,true);
	for(std::vector<DMCblock*>::iterator dit = dc->blocks.begin(); dit != dc->blocks.end(); dit++){
		if((*dit)->Tis50ns_Fis25ns) continue; //don't run 50ns swill
		//if((*dit)->type <0 ) continue; //don't data

		cout<<endl<<"RUNNING BLOCK "<<(*dit)->name<<endl;
	
		clock_t startloop= clock();
		fakeratecalc* t = new fakeratecalc(*dit);

		t->Loop(&EventRegistry,&BadEventRegistry);
		//if((*dit)->type == 71 && (*dit)->processID == 1) cout<<"Found DY 50? "<<(*dit)->name<<endl;
		clock_t endloop= clock();
		cout << "Block "<<(*dit)->name<<" ran in "<<((float) (endloop - startloop))/CLOCKS_PER_SEC<<" seconds"<<endl;
	}//end for all blocks
	clock_t endfor = clock();
	cout << "All blocks ran in "<<((float) (endfor - startfor))/CLOCKS_PER_SEC<<" seconds"<<endl;
}//end runall

