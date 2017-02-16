#include <iostream> 
#include <stack> 
#include <map>
#include <vector>
#include <array>
//#include <string.h>
//#include <exception>
//#include <TH1F.h>
//#include <time.h>
using namespace std;

class particle{
	public:
	    const int pid;
	    int nbranchings; //stable particles have nbranchings = 0
	    vector< array<int,4>> daughters; //First element is the acutal length of the array. then space for 3 pid's 
					     //list of branchings, each branching is a list of daughters. 
	    vector<double> branching_ratios;
		particle(int _pid);
	private:
		void setW(int pid);
		void setZ(int pid);
		void setH();
};

struct node{
	//this is a wrapper for a particle that keeps track of which one we're on
    int which_branch;//keeps track of which branch we're on.
    particle * p; //I need to have the particle not destroyed when the node is destroyed.  //check
	node(particle* _p): p(_p) {which_branch = 0;}
	node(int pid, map<int,particle*>& pmap): p(pmap[pid]){ which_branch = 0;}
	~node(){}; //do not delete p when deleting node. 
	bool increment();//Returns false iff overflow, if overlfow which_branch = 0

};

struct stage{
    double branching_ratios;
    vector<node*> list; //must be extensible, must allow repete objects, 
	stage(){branching_ratios=1.0;};
	bool increment(); //returns true iff overflow
	//bool my_decay(stack<stage* >& sstack, map<int,particle*> pmap); //returns true unless the stage was stable--thus nothing pushed. 
	bool is_stable();
	void print();
	~stage(){ for(vector<node*>::iterator inode = list.begin(), endnode = list.end(); inode< endnode ;inode++) delete *inode; } //delete2
};

/*class analyzer{
	public:
		virtual void consider(stage* astage)=0; //abstract 
		virtual void report()=0; //abstract 
		virtual ~analyzer();
};*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool stage::is_stable(){
	bool top_is_stable = true;
        for(vector<node*>::iterator inode = list.begin(), endnode = list.end(); top_is_stable and inode< endnode ;inode++)
		top_is_stable &= (*inode)->p->nbranchings == 0;
        return top_is_stable;
}

void stage::print(){
        for(vector<node*>::iterator inode = list.begin(), endnode = list.end(); inode< endnode ;inode++)
		cout<<(*inode)->p->pid<<':'<<(*inode)->which_branch<<' '<<endl;
	cout<<endl;
}

inline int sign(int a){return a<0?-1:1;}

particle::particle(int _pid):pid(_pid){
	/*
	   Anti particles, and W-  have pid < 0;
	   pid % 100 (the first two digits) encode the particle. The most significant digits encode the parent.
	   so 2411 is a matter electron from a W+
	   if the digit third from the left is a 0, that means the particle came from a parent by way of an intermediate. 
	   so 25014 is a muon from a higgs by way of a vector boson. 2514 is a muon directly from the higgs. 
	 */

	int fuckyousgn = sign(pid);
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
		case    12://nu
		case    13://mu , 
		case  6013://twm
		case  2313://zm
		case  2413://wm
		case 25013://HVm
		case  2513://Hm
		case    15://tau //listed here as stable for convenience; should be changed. 
		case  6015://twtau
		case  2315://ztau
		case  2415://wtau
		case  2515://Htau
		case 25015://HVtau
		case    21://glue
		case    22://gamma 
			nbranchings = 0;
			branching_ratios.push_back(1.);
			break;
		case     6://t
			nbranchings = 1;
			{
				array<int,4> tbw = {2,fuckyousgn*605,fuckyousgn*624};//b and W-
				daughters.push_back(tbw);
			}
			branching_ratios.push_back(1.);
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
		case    99://T'
			nbranchings = 3;
			{
				array<int,4> bw = {2,fuckyousgn*5,fuckyousgn*24};//b and W-
				array<int,4> tZ = {2,fuckyousgn*6,23};//t and Z
				array<int,4> tH = {2,fuckyousgn*6,25};//t and H
				daughters.push_back(bw);
				branching_ratios.push_back(0.5);
				daughters.push_back(tZ);
				branching_ratios.push_back(0.25);
				daughters.push_back(tH);
				branching_ratios.push_back(0.25);
			}
			break;
		default:
			nbranchings = 0;
			branching_ratios.push_back(1.);
	}//end switch
}//end constructor 

map<int,particle*> makemap(){
	//All particles are born in make map and live in the resulting map.
	//so you only ever make the particle once, then point to it. 
	//particles indexed by their pid via this map.
	map<int,particle*> out;
	int uncharged[] = {25, 23, 2523, 21, 22, 12};
	int uncharged_size = 6;
	int charged[] = {1, 4, 5, 605, 2305, 2505, 11, 6011, 2311, 2411, 25011, 13, 6013, 2313, 2413, 25013, 2513, 15, 6015, 2315, 2415, 2515, 25015, 6, 24, 624, 2524, 99};
	int charged_size = 28;

	for(int i=0;i<uncharged_size;i++) out[uncharged[i]] = new particle(uncharged[i]); //delete1
	for(int i=0;i<charged_size;i++){
		out[charged[i]] = new particle(charged[i]);//delete1
		out[-charged[i]] = new particle(-charged[i]);//delete1
	}
	return out;
}//end makemap
void deletemap(map<int,particle*> & pmap){
	//properly delete the map made by makemap
	for(map<int, particle*>::iterator ip=pmap.begin(); ip!=pmap.end(); ++ip) delete ip->second; //delete1
}

void particle::setW(int pid){
	int prefix = 0;
	int fuckyousgn = sign(pid);
	switch(abs(pid)){
		case   624://tW+
			prefix = 600; //ele from W from t
			break;
		case  2524://HW+
			prefix = 2500; //ele from W from H
			break;
		default:
			prefix = 2400; //ele from W
	}
	array<int,4> qq = {2,1,1};
	daughters.push_back(qq);
	branching_ratios.push_back(0.6775);
	array<int,4> enu = {2,prefix + fuckyousgn*11,12};
	daughters.push_back(enu);
	branching_ratios.push_back(0.1085);
	array<int,4> mnu = {2,prefix + fuckyousgn*13,12};
	daughters.push_back(mnu);
	branching_ratios.push_back(0.1015);
	array<int,4> tnu = {2,prefix + fuckyousgn*15,12};
	daughters.push_back(tnu);
	branching_ratios.push_back(0.1125);
	/* //pi+ &  photon rare decay. 
	array<int,4> pig = {2,1,22};
	daughters.push_back(pig);
	branching_ratios.push_back(0.00022);
	*/
	nbranchings = daughters.size();//update later
}//end setW

void particle::setZ(int pid){
	int prefix = 2300;
	if(abs(pid)==2523) prefix = 25000; // z form higgs
	array<int,4> ee = {2,prefix+11,-prefix-11};
	array<int,4> mm = {2,prefix+13,-prefix-13};
	array<int,4> tt = {2,prefix+15,-prefix-15};
	array<int,4> nn = {2,12,12};//nunu
	array<int,4> qq = {2,1,1};//light hadrons
	array<int,4> cc = {2,4,-4};
	array<int,4> bb = {2,2304,-2304};
	array<int,4> ggg = {3,21,21,21};
	daughters.push_back(qq);
	branching_ratios.push_back(0.3943);//qq
	daughters.push_back(nn);
	branching_ratios.push_back(0.2001);//nunu
	daughters.push_back(bb);
	branching_ratios.push_back(0.1516);//bb
	daughters.push_back(cc);
	branching_ratios.push_back(0.1240);//cc
	daughters.push_back(ee);
	branching_ratios.push_back(0.03366);//ee
	daughters.push_back(mm);
	branching_ratios.push_back(0.03367);//mm
	daughters.push_back(tt);
	branching_ratios.push_back(0.03360);//tt
	daughters.push_back(ggg);
	branching_ratios.push_back(0.0101);//ggg
	//br's sum to 0.98103
	nbranchings = daughters.size();//update later
}//setZ

void particle::setH(){
	//sets to the branching ratios of a 125.5 SM higgs. 
	array<int,4> bb = {2,2504,-2504};
	array<int,4> tt = {2,2515,-2515};
	array<int,4> mm = {2,2513,-2513};
	array<int,4> cc = {2,4,-4};
	//array<int,4> qq = {2,1,1};//light hadrons
	array<int,4> gg = {2,21,21};//gluons
	array<int,4> phopho = {2,22,22};
	array<int,4> Zpho = {2,2523,22};
	array<int,4> WW = {2,2524,-2524};
	array<int,4> ZZ = {2,2523,2523};

	daughters.push_back(bb);
	branching_ratios.push_back(0.569);//bb
	daughters.push_back(WW);
	branching_ratios.push_back(0.223);//WW
	daughters.push_back(gg);
	branching_ratios.push_back(0.0852+2.43E-04);//glue glue + ss
	daughters.push_back(tt);
	branching_ratios.push_back(6.24E-02);//tau tau
	daughters.push_back(cc);
	branching_ratios.push_back(2.87E-02);//cc
	daughters.push_back(ZZ);
	branching_ratios.push_back(2.76E-02);//ZZ
	daughters.push_back(phopho);
	branching_ratios.push_back(2.28E-03);//gamma gamma
	daughters.push_back(Zpho);
	branching_ratios.push_back(1.58E-03);//Zgamma
	daughters.push_back(mm);
	branching_ratios.push_back(2.16E-04);//mu mu 
	nbranchings = daughters.size();//9
}//end setH



bool node::increment(){//wad
	//normall increments which_branching and returns true. 
	//If which_branching then exceeds nbranchigns, reset which_branching to 0 and return false;
	//relies on stable particles having nbranchings = 0
	//cout<<"node increment for pid "<<p->pid<<" which before increment: "<<which_branch<<" nbranches: "<<p->nbranchings<<endl;
	if(++which_branch >= p->nbranchings){
		which_branch = 0;
	//cout<<"node increment overflowed"<<endl;
		return false; //overflow
	}
	//cout<<"node increment simply"<<endl;
	return true; //normal
}//end increment node


bool stage::increment(){//wad
	//increment the particles in the stage sequentially. Returns false unless we increment off the end, in which case it returns true. . 
	for(vector<node*>::iterator inode = list.begin(), lend = list.end();inode< lend;inode++){
		if( (*inode)->increment()){
			//cout<<"stage increment is done incrementing."<<endl;//debug
			return false; //if node doesn't overflow, stop
		}
	}
	//cout<<"stage increment gets to end of list so everything on the list ought to have overflowed."<<endl;//debug
	return true; //overflow
}

bool next(stack<stage*>& sstack){
	//returns false  if next produces an empty stack
	delete sstack.top();//delete3
	sstack.pop();
	if(sstack.empty()) return false;
	sstack.top()->print();//for check. 
	if(sstack.top()->increment()){
		cout<<"increment overflow"<<endl; //COUT potential infinite loop
		return next(sstack); //overflowed, go pop up level. 
	}
	cout<<"increment says overflow"<<endl;
	return true;
}

bool my_decay(stack<stage*>& sstack, map<int,particle*> pmap){ //the name "decay" is taken somewhere in root. 
	//decays the top of the stack into a new stage, and pushes that onto the stack. 
	//returns true unless the stage was stable. 
	cout<<"in my_decay"<<endl;
	if(sstack.top()->is_stable()) return false;

	cout<<"unstable"<<endl;
	stage* s = new stage(); //delete3
	int i=1;//debug
	for(vector<node*>::iterator inode = sstack.top()->list.begin();inode< sstack.top()->list.end();inode++){
		cout<<"considering node "<<i++<<endl;//debug
		if((*inode)->p->nbranchings == 0){	
			s->list.push_back(new node( (*inode)->p ));//delete2
		}
		else{
			array<int,4> thisbranch = (*inode)->p->daughters[(*inode)->which_branch];//a list of a couple ints naming which particles this decays into .
			for(int i=1,n = thisbranch[0];i<=n;++i) s->list.push_back(new node( thisbranch[i], pmap ));//delete2
			s->branching_ratios *= (*inode)->p->branching_ratios[(*inode)->which_branch];
		}
	}//end for every node. 
	cout<< "we have a new stage: "<<endl;
	s->print();
	sstack.push(s);
	return true;
}//end my_decay

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void initalize_stack_with_TpTp(stack<stage*>& sstack, map<int,particle*>& pmap){ 
	stage* out = new stage();//delete3
	out->list.push_back(new node( 99,pmap));//delete2
	out->list.push_back(new node(-99,pmap));//delete22
	sstack.push(out);
}

void initalize_stack_with_Tp(stack<stage*>& sstack, map<int,particle*>& pmap){ 
	stage* out = new stage();//delete3
	out->list.push_back(new node( 99,pmap));//delete2
	sstack.push(out);
}

void initalize_stack_with_Top(stack<stage*>& sstack, map<int,particle*>& pmap){ 
	stage* out = new stage();//delete3
	out->list.push_back(new node(6,pmap));//delete2
	sstack.push(out);
}

void initalize_stack_with_Higgs(stack<stage*>& sstack, map<int,particle*>& pmap){ 
	stage* out = new stage();//delete3
	out->list.push_back(new node(25,pmap));//delete2
	sstack.push(out);
}

void initalize_stack_with_DiHiggs(stack<stage*>& sstack, map<int,particle*>& pmap){ 
	stage* out = new stage();//delete3
	out->list.push_back(new node(25,pmap));//delete2
	out->list.push_back(new node(25,pmap));//delete2
	sstack.push(out);
}

void Decay_until_stable(stack<stage*> &sstack,map<int,particle*>& pmap){
	//recursively decays the top of the stack until the result of the decay is stable 
	if(my_decay(sstack,pmap)){
		cout<<"decaying "; //COUT potential infinite loop
		sstack.top()->print();

		Decay_until_stable(sstack,pmap); //decay sstack.top if stable, then decay that.
	}
}


class lepton_counter{
//class lepton_counter: public analyzer{
	private: 
		long n_leaves;
		double br_2lep;
		double net_br;

	public:
		void consider(stage* astage);
		void report(){
			cout<<"Number of leaves: "<<n_leaves<<" net br: "<<net_br<<" dilepton br "<<br_2lep<<endl;
		}
		lepton_counter(){
			n_leaves=0;
			br_2lep=0.;
			net_br=0.;
		}
};

void lepton_counter::consider(stage* astage){
	++n_leaves;
	net_br += astage->branching_ratios;
	int nlep = 0;
        for(vector<node*>::iterator inode = astage->list.begin(), endnode = astage->list.end(); inode< endnode ;inode++){
                int pid = abs( (*inode)->p->pid ) % 100;
		if(pid == 11 or pid == 13) ++nlep;
	}
	if(nlep ==2) br_2lep += astage->branching_ratios;
}


void ParticleBranchRatioCounter(){ //the main function. 
	map<int,particle*> pmap = makemap(); 
	stack<stage*> sstack;

	//say what to start with
	//initalize_stack_with_DiHiggs(sstack,pmap);
	//initalize_stack_with_Higgs(sstack,pmap);
	initalize_stack_with_Top(sstack,pmap);
//void initalize_stack_with_Top(stack<stage*>& sstack, map<int,particle*>& pmap){ 
	//initalize_stack_with_Tp(sstack,pmap);
	//initalize_stack_with_TpTp(sstack,pmap);

	//specify an alyzer
	//analyzer* Analyser = new lepton_counter(); //collects data about the leaves of the decay tree. //delete4
	lepton_counter* Analyser = new lepton_counter(); //collects data about the leaves of the decay tree. //delete4

	//go analyzer
	do{
		sstack.top()->print();
		Decay_until_stable(sstack,pmap);
		Analyser->consider(sstack.top());
	}while(next(sstack));

	//say what you found. 
	Analyser->report();

	//cleanup
	deletemap(pmap);
	delete Analyser;//delete4
}//end ParticleBranchRatioCounter
