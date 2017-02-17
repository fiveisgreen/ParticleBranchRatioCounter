#include <iostream> 
#include <stack> 
#include <vector>
#include <array>
#include "interface/intmap.h"
//#include <TH1F.h>

//some precompiler controls 
#define verbose 0   //turns on all the print statements to give a running history. Don't use this for complex decays.
#define WATCH_THE_CLOCK 1   //Turns on a clock to monitor how long the computation takes

/***************************************************************************
 * This is a calculator for considering all on-shell decay branching
 * fractions of complex decays, such as di-Higgs or Tprime pair produceion. 
 * The main calculator, ParticleBranchRatioCounter, will decay all 
 * combinations of on shell decay's. 
 *
 * ParticleBranchRatioCounter takes as input 1 or 2 particle pid's 
 * which will be the inital particles to be decayed. 
 *
 * Each combination of stable particle is supplied to 
 * analyser::consider and tallied. This should be the main point 
 * of interest for editing since all sorts of things can be tallied about 
 * physics siganls. Edit "consider" and the coutners in analyser 
 * to taste. 
 *
 * At the end of the calculation, analyser::report gives a report
 * of what was tellied. 
 *
 * Created by Dr. Anthony Barker, January 2017. 
 * *************************************************************************/

#define LEN_D_ARR 4
#if WATCH_THE_CLOCK== 1   //Turns on a clock to monitor how long the computation takes
#include <time.h>
#endif
using namespace std;

class particle{
	//Here lies all the physics data. 
	//particles store what particle this is, some of it's pedigree, what it can decay into, and correspoing branching fractions. 
    public:
	const int pid; 
	/*
	   a numeric name for this particle. Negative pid are anti-particles or W-. 
	These generally follow the Pythia naming convention, with a variety of differences. 
	Most notably, we include some of the particle's pedigree and some particles don't have antiparticles right now. 
	particle 1 represents all light flavor quark jets: u,d,s
	particle 12 represents all neutrinos.
	The left two digits are used to store particle identity: you can use abs(pid%100) to strip the pedigree and charge information. 
	If digit third from the left is 0, then the particle came from a W or Z resulting from a more interesting grandmother particl described by the left 2 digits.
	If digit third from the left is 5, then the particle came from a tau resulting from a more interesting grandmother particl described by the left 2 digits.
	The left two digits encode what mother or grandmother particle generated this particle. 
	Examples: 
	-13 is an anti-muon. 
	-2513 is an antimuon whose mother is a higgs boson (pid 25 = higgs). 
	-25013 is an higgs->V->antimuon. 
	-6013 is an anti_top->W->antimuon (pid 6 = top). 
	See the constructor for the full list. 
	*/

	int nbranchings; 
	//number of decay channels for this particle. Stable particles have nbranchings = 0

	vector< array<int,LEN_D_ARR>> daughters; // Stores pID's of all possible daughter particles. 
	//vector (lenght = nbranchings) of arrays of pID's 
	//Except the first element of each array is the acutal length of the array, followed by up to 3 pID's 

	vector<double> branching_fraction; //length = nbranchings

	//Particle physics data about particles decays is hard coded into the constructor and these set functions: 
	particle(int _pid); 
    private:
	void setW(int pid); //used to clean up the constructor a bit. 
	void setZ(int pid);
	void setTau(int pid,bool decay_taus);
	void setH();
}; // end class particle

class pIDmap: public intmap<particle*>{
    public:
	pIDmap(){}
	virtual ~pIDmap(){ for( std::map<int,particle*>::iterator ip = tmap.begin(); ip != tmap.end(); ip++) delete ip->second; } //delete1
};

particle* mapfetch(int pid); 
//mapfetch efficiently delivers a particle* corresponding to a provided pID. 
//If you give it a bad pid the program does a controled termination and tells you that you gave it a bad pid. 

struct node{
    //this is a wrapper for a particle that keeps track of which one we're on
    //particles are global, nodes are spicific to a moment in the decay chain. Many nodes point to the same particle. 
    int which_branch;//keeps track of which decay channel we're on.
    particle * p; //is not destroyed when the node is destroyed. 

    node(particle* _p): p(_p), which_branch(0) {}
    node(int pid): p(mapfetch(pid)), which_branch(0) {}
    ~node(){}; //deletes nodes without deleting the particle. 
    bool increment();//Returns false iff overflow, if overlfow which_branch = 0
};

struct stage{
	//A stage is a stage in the decay sequence. It is a list of particles, potentially a mix of stable and unstable, together 
	//with the branching fraction for arriving there. 
    double branching_fraction;
    vector<node*> list; //must be extensible, must allow repete objects, 
	//stage(double br_previous = 1.0):branching_fraction(br_previous){};
	stage(double br_previous = 1.0):branching_fraction(br_previous){list.reserve(8);} //make space for 8 pointers initially.
	bool increment(); //returns true iff overflow
	bool is_stable();
	void print();
	~stage(){ for(auto inode = list.begin(), endnode = list.end(); inode< endnode ;inode++) delete *inode; } //delete2
};

class analyser{
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
		analyser(){
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
};//end class analyser

void analyser::consider(stage* astage){
		#if verbose==1
		//cout<<"$$GOT LEAF$$ "; astage->print(); cout<<endl;
		#endif
    ++n_leaves;
    net_br += astage->branching_fraction;
	//basics

    int nlep = 0;
    //int nB = 0;
    //int nJ = 0; //light jets
    for(auto inode = astage->list.begin(), endnode = astage->list.end(); inode< endnode ;inode++){
	int pid = abs( (*inode)->p->pid ) % 100;
	if(pid == 11 or pid == 13) ++nlep; 
//	if(pid == 5) nB++;
//	if(pid == 1 or pid==4 or pid == 21 ) nJ++;
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool stage::is_stable(){
    //This is a performance critical function. 
	bool top_is_stable = true;
        for(auto inode = list.begin(), endnode = list.end(); top_is_stable and inode< endnode ;inode++)
		top_is_stable &= (*inode)->p->nbranchings == 0;
        return top_is_stable;
}//end is_stable

void stage::print(){ 
	//This prints out some data about a stage
	//pid:which_branch pid:which_branch ... stage's_branching_fraction
	//turn on verbose to make it go. 
        for(auto inode = list.begin(), endnode = list.end(); inode< endnode ;inode++){
		cout<<(*inode)->p->pid<<':'<<(*inode)->which_branch<<' '; //main line 
	}
	cout<<"\t\t\t BR: "<<branching_fraction<<endl;
}//end print

inline int sign(int a){return a<0?-1:1;}

particle::particle(int _pid):pid(_pid){
	/* Anti particles, and W-  have pid < 0;
	   pid % 100 (the first two digits) encode the particle. The most significant digits encode the parent.
	   so 2411 is a matter electron from a W+
	   if the digit third from the left is a 0, that means the particle came from a parent by way of an intermediate boson. 
	   if the digit third from the left is a 5, that means the particle came from a parent by way of an intermediate tau. 
	   so 25014 is a muon from a higgs by way of a vector boson. 2514 is a muon directly from the higgs. 
	
	   If any changes are made here they must by synchrnoized by hand with makemap.
	   This is not a performance critical function. 
	   Particles only created at the start of the calculation, and only one of each type. 
	 */

	int i_sign = sign(pid);
	int apid = abs(pid);
	switch(apid){
		//here define what particles are and their branchings. 
		case     1://q
		case     4: //c,
		case     5: //b,
		case   605: //tb
		case  2305: //Zb
		case  2505: //Hb
		case    11://ele 
		case  6011://twe
		case  2311://ze
		case  2411://we
		case 25011://HVe
		case   511://tau ele
		case  6511://top->tau->ele
		case 23511://Z->tau->ele excluding Z from H
		case 24511://W->tau->ele excluding W from H
		case 25511://H->tau->ele
		case    12://nu
		case    13://mu , 
		case  6013://twm
		case  2313://zm
		case  2413://wm
		case 25013://HVm
		case  2513://Hm
		case   513://tau mu
		case  6513://top->tau->mu
		case 23513://Z->tau->mu, excluding Z from H
		case 24513://W->tau->mu, excluding W from H
		case 25513://H->tau->mu, 
		case    21://glue
		case    22://gamma 
			nbranchings = 0; //stable
			branching_fraction.push_back(1.);
			break;
		case    15://tau //listed here as stable for convenience; should be changed. 
		case  6015://twtau
		case  2315://ztau
		case  2415://wtau
		case  2515://Htau
		case 25015://HVtau
			setTau(pid,true);
			break;
		case     6://t
			nbranchings = 1;
			{
				array<int,LEN_D_ARR> tbw = {2,-i_sign*605,i_sign*624};//b and W-
				daughters.push_back(tbw);
			}
			branching_fraction.push_back(1.);
			break;
		case    23://Z
		case  2523://HZ
			setZ(pid);
			break;
		case    24://W
		case   624://tW
		case  2524://HW
			setW(pid);
			break;
		case    25://H
			setH();
			break;
		case    55://B'
			nbranchings = 3;
			{
				array<int,LEN_D_ARR> tw = {2,-i_sign*6,i_sign*24};//t and W-
				array<int,LEN_D_ARR> bZ = {2,i_sign*5,23};//b and Z
				array<int,LEN_D_ARR> bH = {2,i_sign*5,25};//b and H
				daughters.push_back(tw);
				branching_fraction.push_back(0.5);
				daughters.push_back(bZ);
				branching_fraction.push_back(0.25);
				daughters.push_back(bH);
				branching_fraction.push_back(0.25);
			}
			break;
		case    66://T'
			nbranchings = 3;
			{
				array<int,LEN_D_ARR> bw = {2,-i_sign*5,i_sign*24};//b and W-
				array<int,LEN_D_ARR> tZ = {2,i_sign*6,23};//t and Z
				array<int,LEN_D_ARR> tH = {2,i_sign*6,25};//t and H
				daughters.push_back(bw);
				branching_fraction.push_back(0.5);
				daughters.push_back(tZ);
				branching_fraction.push_back(0.25);
				daughters.push_back(tH);
				branching_fraction.push_back(0.25);
			}
			break;
		default:
			nbranchings = 0;
			branching_fraction.push_back(1.);
	}//end switch
}//end constructor 

pIDmap* makemap(){
	//All particles are born in make map and live in the resulting map.
	//so you only ever make the particle once, then point to it. 
	//particles indexed by their pid via this map.
	pIDmap* out = new pIDmap();
	array<int,6>uncharged = {25, 23, 2523, 21, 22, 12};
	array<int,38> charged = {1, 4, 5, 605, 2305, 2505, 11, 6011, 2311, 2411, 25011,511, 6511, 23511, 24511, 25511, 513, 6513, 23513, 24513, 25513, 13, 6013, 2313, 2413, 25013, 2513, 15, 6015, 2315, 2415, 2515, 25015, 6, 24, 624, 2524, 66};
		   
	//int uncharged[] = {25, 23, 2523, 21, 22, 12};
	//int uncharged_size = 6;
	//int charged[] = {1, 4, 5, 605, 2305, 2505, 11, 6011, 2311, 2411, 25011, 13, 6013, 2313, 2413, 25013, 2513, 15, 6015, 2315, 2415, 2515, 25015, 6, 24, 624, 2524, 66};
	//int charged_size = 28;

	for(unsigned int i=0;i<uncharged.size();i++) out->set(uncharged[i], new particle(uncharged[i]) ); //delete1 //doesn't compile
	for(unsigned int i=0;i<charged.size();i++){
		out->set(charged[i],  new particle( charged[i]) ); //delete1
		out->set(-charged[i], new particle(-charged[i]) ); //delete1
	}
	return out;
}//end makemap

particle* mapfetch(int pid){
	//safely serves the particle pointers given the PID. 
	//This is a performance critical function. 
    static pIDmap* pmap = makemap(); //only done once at startup. 

    particle* p;
    try{      p = pmap->get_throwable(pid,1);      }

    catch(std::pair <int,int> errorpair){
	cerr<<"Error in mapfetch: Invalid pID "<<errorpair.first<< " Check makemap and particle constructor correspondence"<<endl;
	std::terminate();
    }
    return p;
}//end mapfetch.

void particle::setTau(int pid,bool decay_taus){
    if(not decay_taus){
	nbranchings = 0; //stable
	branching_fraction.push_back(1.);
	return;
    }

    int i_sign = sign(pid);
    int prefix;

    switch(abs(pid)){
	case    15://tau 
	    prefix = 500; //lep from tau
	    break;
	case  6015://twtau
	    prefix = 6500; //lep from tau from W from top
	case  2315://ztau
	    prefix = 23500; //lep from tau from Z
	    break;
	case  2415://wtau
	    prefix = 24500; //lep from tau from W
	    break;
	case  2515://Htau
	case 25015://HVtau
	    prefix = 25500; //lep from tau from H
	    break;
	default:
	    prefix = 0; //ele from W
    }

	//Set Branches; 
	//	Has any one ever told you tau's are messy? 
	//	Taus tend to have multi-partile decays, split this up for whatever number of decays are being modeled 

	
	//Leptonic Decays
			#if LEN_D_ARR==3 //only 2 particle decays allowed, l nu nu becomes l nu; mu nu nu gamma is neglected.
	//e- nn	0.1781 as e nu
	array<int,LEN_D_ARR> en = {2,i_sign*(prefix+11), 12}; //ele nu representing ele nu nu
	daughters.push_back(en);
	branching_fraction.push_back(0.1781);

	//mu- nn	0.1707
	//mu- nng	0.003 
	array<int,LEN_D_ARR> mn = {2,i_sign*(prefix+13), 12}; //mu nu representing mu nu nu and mu nu nu gamma (rare)
	daughters.push_back(mn);
	branching_fraction.push_back(0.1737);

			#elif LEN_D_ARR>=4 //3 or more particles allowed in decay
	//e- nn	0.1781
	array<int,LEN_D_ARR> enn = {3,i_sign*(prefix+11), 12, 12}; //had nu
	daughters.push_back(enn);
	branching_fraction.push_back(0.1781);

	//mu- nn	0.1707
	array<int,LEN_D_ARR> mnn = {3,i_sign*(prefix+13), 12, 12}; //had nu
	daughters.push_back(mnn);
	branching_fraction.push_back(0.1707);

	//mu- nng	0.003 
				#if LEN_D_ARR==4
	array<int,LEN_D_ARR> mng = {3,i_sign*(prefix+13), 12, 22}; //mu nu nu gamma represented as mu nu gamma
	daughters.push_back(mng);
	branching_fraction.push_back(0.003);
				#else
	array<int,LEN_D_ARR> mnng = {4,i_sign*(prefix+13), 12, 12, 22}; //mu nu nu gamma fully represented. 
	daughters.push_back(mnng);
	branching_fraction.push_back(0.003);
				#endif
			#endif  //end lepton compiler if


	//Hadronic Decays
			#if LEN_D_ARR<=3 //only 2 daughters allowed, model all hadronic decays as hnu
	array<int,LEN_D_ARR> hnu = {2,i_sign*1, 12}; //had nu accounts for all 
	daughters.push_back(hnu);
	branching_fraction.push_back(0.6482);

			#else //more than 2 daughters allowed, break up hadronic decays into single hadron and multi-hadron.
	//1 had- nu	0.4952
	array<int,LEN_D_ARR> hnu = {2,i_sign*1, 12}; //had nu
	daughters.push_back(hnu);
	branching_fraction.push_back(0.4952);
			#endif

			#if LEN_D_ARR==4 //3 decay products allowed
	array<int,LEN_D_ARR> hhnu = {3,i_sign*1,i_sign*1, 12}; //3 had nu modeled as h nu
	daughters.push_back(hhnu);
	branching_fraction.push_back(0.1530); //4 decay products allowed, can fully model tau to hhhnu

			#elif LEN_D_ARR==5 //5 particles allowed in the decay.
	array<int,LEN_D_ARR> hhhnu = {4,i_sign*1,i_sign*1, -i_sign*1, 12}; //fully modeled. 
	daughters.push_back(hhhnu);
	branching_fraction.push_back(0.1530); //5h nu subsumed in 3h nu

			#elif LEN_D_ARR==6 //6 particles allowed in the decay.
	array<int,LEN_D_ARR> hhhnu = {4,i_sign*1,i_sign*1, -i_sign*1, 12}; //fully modeled. 
	daughters.push_back(hhhnu);
	branching_fraction.push_back(0.15203); 

	array<int,LEN_D_ARR> h5nu = {5,i_sign*1,i_sign*1, -i_sign*1,i_sign*1, 12}; //model 5h nu as 4h nu
	daughters.push_back(h5nu);
	branching_fraction.push_back(0.00097);

			#elif LEN_D_ARR>=7 //7 or more particles allowed in the decay.
	array<int,LEN_D_ARR> hhhnu = {4,i_sign*1,i_sign*1, -i_sign*1, 12}; //fully modeled. 
	daughters.push_back(hhhnu);
	branching_fraction.push_back(0.15203); 

	array<int,LEN_D_ARR> h5nu = {6,i_sign*1,i_sign*1, -i_sign*1,i_sign*1, -i_sign*1, 12}; //fully modeled. 
	daughters.push_back(h5nu);
	branching_fraction.push_back(0.00097);
			#endif
		
	nbranchings = daughters.size();//update later
}//end setTau

void particle::setW(int pid){
	//This sets the branches, daughters, and branching fractions of W
	//while preserving the grandmother information. 
	//This is NOT a performance critical function. 
	int prefix = 0;
	int i_sign = sign(pid);
	switch(abs(pid)){
		case   624://tW+
			prefix = 6000; //ele from W from t
			break;
		case  2524://HW+
			prefix = 25000; //ele from W from H
			break;
		default:
			prefix = 2400; //ele from W
	}
	array<int,LEN_D_ARR> qq = {2,1,1};
	daughters.push_back(qq);
	branching_fraction.push_back(0.6775);
	array<int,LEN_D_ARR> enu = {2,i_sign*(prefix + 11),12};
	daughters.push_back(enu);
	branching_fraction.push_back(0.1085);
	array<int,LEN_D_ARR> mnu = {2,i_sign*(prefix + 13),12};
	daughters.push_back(mnu);
	branching_fraction.push_back(0.1015);
	array<int,LEN_D_ARR> tnu = {2,i_sign*(prefix + 15),12};
	daughters.push_back(tnu);
	branching_fraction.push_back(0.1125);
	 //pi+ &  photon rare decay. //you may comment out the next three lines if you like, speeds it up. 
	array<int,LEN_D_ARR> pig = {2,1,22};
	daughters.push_back(pig);
	branching_fraction.push_back(0.00022);
	nbranchings = daughters.size();//update later
}//end setW

void particle::setZ(int pid){
	//This sets the branches, daughters, and branching fractions of Z
	//while preserving the grandmother information. 
	//This is NOT a performance critical function. 
	int prefix = 2300;
	if(abs(pid)==2523) prefix = 25000; // z form higgs
	array<int,LEN_D_ARR> ee = {2,prefix+11,-prefix-11};
	array<int,LEN_D_ARR> mm = {2,prefix+13,-prefix-13};
	array<int,LEN_D_ARR> tt = {2,prefix+15,-prefix-15};
	array<int,LEN_D_ARR> nn = {2,12,12};//nunu
	array<int,LEN_D_ARR> qq = {2,1,-1};//light hadrons
	array<int,LEN_D_ARR> cc = {2,4,-4};
	array<int,LEN_D_ARR> bb = {2,2305,-2305};
			#if LEN_D_ARR>=4
	array<int,LEN_D_ARR> ggg = {3,21,21,21}; //glue glue glue. this is the only listed 3 particle decay.
			#else 
	array<int,LEN_D_ARR> ggg = {2,21,21}; //model glue glue glue as 2 gluons 
			#endif
	daughters.push_back(qq);
	branching_fraction.push_back(0.3943);//qq
	daughters.push_back(nn);
	branching_fraction.push_back(0.2001);//nunu
	daughters.push_back(bb);
	branching_fraction.push_back(0.1516);//bb
	daughters.push_back(cc);
	branching_fraction.push_back(0.1240);//cc
	daughters.push_back(ee);
	branching_fraction.push_back(0.03366);//ee
	daughters.push_back(mm);
	branching_fraction.push_back(0.03367);//mm
	daughters.push_back(tt);
	branching_fraction.push_back(0.03360);//tt
	daughters.push_back(ggg);
	branching_fraction.push_back(0.0101);//ggg
	//br's sum to 0.98103
	nbranchings = daughters.size();//update later
}//setZ

void particle::setH(){
	//This sets the branches, daughters, and branching fractions of Z
	//while preserving the grandmother information. 
	//Branching ratios correspond to a 125.5 GeV SM Higgs 
	//This is NOT a performance critical function. 
	array<int,LEN_D_ARR> bb = {2,2505,-2505};
	array<int,LEN_D_ARR> tt = {2,2515,-2515};
	array<int,LEN_D_ARR> mm = {2,2513,-2513};
	array<int,LEN_D_ARR> cc = {2,4,-4};
	//array<int,LEN_D_ARR> ss = {2,1,1};//strange strange, merged with glue glue since both produce light flavor jets.
	array<int,LEN_D_ARR> gg = {2,21,21};//gluons
	array<int,LEN_D_ARR> phopho = {2,22,22};
	array<int,LEN_D_ARR> Zpho = {2,2523,22};
	array<int,LEN_D_ARR> WW = {2,2524,-2524};
	array<int,LEN_D_ARR> ZZ = {2,2523,2523};

	daughters.push_back(bb);
	branching_fraction.push_back(0.569);//bb
	daughters.push_back(WW);
	branching_fraction.push_back(0.223);//WW
	daughters.push_back(gg);
	branching_fraction.push_back(0.0852+2.43E-04);//glue glue + ss
	daughters.push_back(tt);
	branching_fraction.push_back(6.24E-02);//tau tau
	daughters.push_back(cc);
	branching_fraction.push_back(2.87E-02);//cc
	daughters.push_back(ZZ);
	branching_fraction.push_back(2.76E-02);//ZZ
	daughters.push_back(phopho);
	branching_fraction.push_back(2.28E-03);//gamma gamma
	daughters.push_back(Zpho);
	branching_fraction.push_back(1.58E-03);//Zgamma
	daughters.push_back(mm);
	branching_fraction.push_back(2.16E-04);//mu mu 
	nbranchings = daughters.size();//9
}//end setH

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool node::increment(){
	//Normall increments which_branching and returns true. 
	//If which_branching then exceeds nbranchigns, reset which_branching to 0 and return false;
	//relies on stable particles having nbranchings = 0
	//This is a perforance critical function
	if(++which_branch >= p->nbranchings){
		which_branch = 0;
		return false; //overflow
	}
	return true; //normal
}//end increment node

bool stage::increment(){
	//increment the particles in the stage sequentially. Returns false unless we increment off the end, in which case it returns true. . 
	//This is a perforance critical function
	for(auto inode = list.begin(), lend = list.end();inode< lend;inode++) if( (*inode)->increment()) return false; //if node doesn't overflow, stop
	return true; //overflow
}

bool next(stack<stage*>& sstack){ 
	//returns false if next produces an empty stack
	//This is a perforance critical function
	delete sstack.top();//delete3
	sstack.pop();
	if(sstack.empty()) return false;
			#if verbose==1
			cout<<'^'; sstack.top()->print();//for check. 
			#endif
	if(sstack.top()->increment()){
			#if verbose==1
			cout<<"increment overflow"<<endl; //COUT potential infinite loop
			#endif
		return next(sstack); //overflowed, go pop up level. 
	}
        		#if verbose==1
			cout<<"increment says overflow"<<endl;
			#endif
	return true;
}//end next

bool decay_stage(stack<stage*>& sstack){ //the name "decay" is taken somewhere in root. 
	//decays the top of the stack into a new stage, and pushes that onto the stack. 
	//returns true unless the stage was stable. 
	//This is a perforance critical function
	if(sstack.top()->is_stable()) return false;

	stage* s = new stage(sstack.top()->branching_fraction); //delete3
	for(auto inode = sstack.top()->list.begin();inode< sstack.top()->list.end();inode++){
	    if((*inode)->p->nbranchings == 0) s->list.push_back(new node( (*inode)->p ));//delete2
	    else{ 
		particle* P = (*inode)->p;

		for(int i=1,n = P->daughters[(*inode)->which_branch][0];i<=n;++i)
		    s->list.push_back(new node( P->daughters[(*inode)->which_branch][i]));//delete2 //n may be 2 or 3. 

		s->branching_fraction *= P->branching_fraction[(*inode)->which_branch];
	    }
	}//end for every node. 
	sstack.push(s);
	return true;
}//end decay_stage

void Decay_until_stable(stack<stage*> &sstack){ 
	//recursively decays the top of the stack until the result of the decay is stable 
	//This is a perforance critical function
        		#if verbose==1
			cout<<"decaing      "; sstack.top()->print();
			#endif
	if(decay_stage(sstack)){
			#if verbose==1
			cout<<"decay result "; sstack.top()->print();
			#endif
	    Decay_until_stable(sstack); //decay sstack.top if stable, then decay that.
	}
}//end Decay_until_stable

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initalize_stack_with_4_particles(stack<stage*>& sstack, int pid1, int pid2, int pid3, int pid4){ 
	//This is NOT a perforance critical function
	//Use this for two unlike particles and two uncharged particles. 
	//initalizes the first stage of the stack with 2 parciles.
	stage* out = new stage();//delete3
	out->list.push_back(new node(pid1));//delete2
	out->list.push_back(new node(pid2));//delete2
	out->list.push_back(new node(pid3));//delete2
	out->list.push_back(new node(pid4));//delete2
	sstack.push(out);
}

void initalize_stack_with_3_particles(stack<stage*>& sstack, int pid1, int pid2, int pid3){ 
	//This is NOT a perforance critical function
	//Use this for two unlike particles and two uncharged particles. 
	//initalizes the first stage of the stack with 2 parciles.
	stage* out = new stage();//delete3
	out->list.push_back(new node(pid1));//delete2
	out->list.push_back(new node(pid2));//delete2
	out->list.push_back(new node(pid3));//delete2
	sstack.push(out);
}

void initalize_stack_with_2_particles(stack<stage*>& sstack, int pid1, int pid2){ 
	//This is NOT a perforance critical function
	//Use this for two unlike particles and two uncharged particles. 
	//initalizes the first stage of the stack with 2 parciles.
	stage* out = new stage();//delete3
	out->list.push_back(new node(pid1));//delete2
	out->list.push_back(new node(pid2));//delete22
	sstack.push(out);
}

void initalize_stack_with_1_particle(stack<stage*>& sstack, int pid1){ 
	//This is NOT a perforance critical function
	//initalizes the first stage of the stack with a single Tprime 
	stage* out = new stage();//delete3
	out->list.push_back(new node( pid1));//delete2
	sstack.push(out);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ParticleBranchRatioCounter(int pid1,int pid2=0, int pid3=0,int pid4=0){
    //the main function. 
    //The contents of the do-loop are perforance critical.
    //The other sectoins of the function are not performance critical. 
			#if WATCH_THE_CLOCK==1 
			clock_t start_clock = clock();
			#endif

    stack<stage*> sstack;
    if(pid4!=0) initalize_stack_with_4_particles(sstack, pid1, pid2, pid3, pid4);
    else if(pid3!=0) initalize_stack_with_3_particles(sstack, pid1, pid2, pid3 );
    else if(pid2!=0) initalize_stack_with_2_particles(sstack, pid1, pid2); 
    else initalize_stack_with_1_particle( sstack, pid1); 

    //specify an alyzer
    analyser* Analyser = new analyser(); //collects data about the leaves of the decay tree. //delete4

    //go analyze
    do{ 	
	       		#if verbose==1
			cout<<"main loop    ";sstack.top()->print();
			#endif
	Decay_until_stable(sstack);
	Analyser->consider(sstack.top());
    }while(next(sstack));

    //report what you found. 
    Analyser->report(); 

    //cleanup
    delete Analyser;//delete4
			#if WATCH_THE_CLOCK==1 
			clock_t stop_clock = clock();
			cout << "Computation took "<< ((float) (stop_clock - start_clock))/CLOCKS_PER_SEC<<" seconds"<<endl;
			#endif
}//end ParticleBranchRatioCounter
