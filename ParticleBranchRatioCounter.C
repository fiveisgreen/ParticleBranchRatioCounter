#include <iostream> 
#include "ParticleBranchRatioCounter.h"

/***************************************************************************
 * This is a calculator for considering all on-shell decay branching
 * fractions of complex decays, such as di-Higgs or Tprime pair produceion. 
 * The main calculator, ParticleBranchRatioCounter, will decay all 
 * combinations of on shell decay's. 
 *
 * ParticleBranchRatioCounter takes as input 1 to 4 particle pid's 
 * which will be the inital particles to be decayed. 
 *
 * Each combination of stable particle is supplied to 
 * analyzer::consider and tallied. This should be the main point 
 * of interest for editing since all sorts of things can be tallied about 
 * physics siganls. Plug in different analyzers using polymorphism 
 * to shape "consider" and "report" to the question at hand.
 *
 * At the end of the calculation, analyzer::report gives a report
 * of what was tellied. 
 *
 * SM Particle branching ratios are based on the 2016 PDG. 
 * Higgs decays are figured for M(Higgs) = 125.5 GeV 
 *
 * Created by Dr. Anthony Barker, January 2017. 
 * *************************************************************************/

using namespace std;

//HERE IS AN EXAMPLE ANALYSER IMPLEMTATION. 
//MAKE YOUR OWN AND PLUG IT INTO 
//
/*class dilepAna : public analyzer{
	//required elements:
       	//constructor, destructor 
        //void consider(stage* astage) //function to be run on every leaf. 
        //void report() //end of computation print out report

	private: 
		long n_leaves;
		long n_leavesDL;
		double net_br;
		double br_0lep;
		double br_1lep;
		double br_2lep;
		double br_3lep;
		double br_4lep;
		double br_OSSF;
		double br_SSSF;
		double br_OSOF;
		double br_SSOF;

		//0, 1, 2 top

	public:
		void consider(stage* astage);
		void report(){
			cout<<"Number of leaves: "<<n_leaves<<" net br: "<<net_br<<" N dilepton leaves "<< n_leavesDL<< " br 0 lep " << br_0lep<<endl;
			cout<< "BRs: 2lep      OSDL      OSSF      OSOF      SSDL      SSSF      SSOF      1lep      3lep      >3lep"<<endl;
			printf("   %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", br_2lep,
				br_OSSF+br_OSOF, br_OSSF, br_OSOF,
				br_SSSF+br_SSOF, br_SSSF, br_SSOF,
				br_1lep, br_3lep, br_4lep);

		}
		dilepAna(){
			n_leaves=0;
			n_leavesDL=0;
			net_br=0.;
			br_0lep=0.;
			br_1lep=0.;
			br_2lep=0.;
			br_3lep=0.;
			br_4lep=0.;
			br_OSSF=0.;
			br_SSSF=0.;
			br_OSOF=0.;
			br_SSOF=0.;
		}
		~dilepAna(){}
};//end class analyzer

void dilepAna::consider(stage* astage){
		#if verbose==1
		//cout<<"$$GOT LEAF$$ "; astage->print(); cout<<endl;
		#endif
    ++n_leaves;
    net_br += astage->branching_fraction;
	//basics

    int nlep = 0;
    for(auto inode = astage->list.begin(), endnode = astage->list.end(); inode< endnode ;inode++){
	int pid = abs( (*inode)->p->pid ) % 100;
	if(pid == 11 or pid == 13) ++nlep; 
    }

    if(nlep==0) br_0lep += astage->branching_fraction;
    else if(nlep==1) br_1lep += astage->branching_fraction;
    else if(nlep ==2){// and nB >= 1 and nJ >= 2 ){
	br_2lep += astage->branching_fraction;
	int pid1 = 0, pid2=0;
	for(auto inode = astage->list.begin(), endnode = astage->list.end(); inode< endnode ;inode++){
	    int pid = abs( (*inode)->p->pid ) % 100;
	    if(pid == 11 or pid == 13){
		if(pid1==0) pid1= (*inode)->p->pid;
		else pid2= (*inode)->p->pid;
	    }
	}//end second loop
	if(pid1*pid2 > 0){ //same sign
	    ++n_leavesDL;
	    //astage->print();
		if(pid1%100 == pid2%100) br_SSSF += astage->branching_fraction;
		else br_SSOF += astage->branching_fraction;
	}
	else{ //opposite sign
		if(abs( pid1 ) % 100 == abs( pid2 ) % 100) br_OSSF += astage->branching_fraction;
		else br_OSOF += astage->branching_fraction;
	}
    }//end dielpton
    else if(nlep ==3){// and nB >= 1 and nJ >= 2 ){
	br_3lep += astage->branching_fraction;
    }
    else if(nlep >=4){// and nB >= 1 and nJ >= 2 ){
	br_4lep += astage->branching_fraction;
    }
}//end coinsider
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
enum whichAnalyzer:int{nlep, signflavor};

void ParticleBranchRatioCounter(whichAnalyzer WhichAnalyzer, int pid1,int pid2=0, int pid3=0,int pid4=0){
    //the main UI function. 
    //not performance critical, just be careful how much work you give to the analyzer. 


    //specify an alayser. Plug in custom analyzers here. 
    analyzer* Analyzer; //collects data about the leaves of the decay tree. //delete4

    switch(WhichAnalyzer){
	case signflavor:
	case default: 
	    Analyzer = new dilepAna(); 
    }


    //Perform the main analysis: 
    ParticleBranchRatioCounter_core(Analyzer, pid1,pid2, pid3,pid4);

    //report what you found. 
    Analyzer->report(); 

    //cleanup
    delete Analyzer;//delete4
}//end ParticleBranchRatioCounter


