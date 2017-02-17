#ifndef INTMAP
#define INTMAP
#include <map>
#include <string.h>
#include <iostream>
#include <exception>
using namespace std;

//intmap is a wrapper class for map<int, whatever>
//with guards against segfaults from calling ints not in the map. 
//If the index int is not found you can default or throw an exception. 

//Example use:
//intmap<int> mymapthing();
//mymapthing.set("green",5);
//int whatisfive = mymapthing.getquick("green",-1);
//try{
//   cout<<"dangerous critter "<<mymapthing.get_throwable("donkey",2)<<endl;
//}
//catch(std::pair <int,int> errorpair){
//   cerr<<"Error in intmap: Invalid int key "<<errorpair.first<< " code "<<errorpair.second<<endl;
//   std::terminate();
//}
//
//Example Map Looping 
//intmap<bool>* kinvarSwitches = new intmap<bool>();
//for( std::map<int,KinematicVar *>::iterator thisKinvar = allKinVars->begin(); thisKinvar != allKinVars->end(); thisKinvar++) {
//            kinvarSwitches->set(thisKinvar->second->tag,false); 
//}

template <class t> class intmap{
	public:
		void set(int intIndex, t item){ tmap[intIndex] = item; }
		t get(int intIndex, t default_); //returns default_ if intIndex not found; no throws guarente. 
		t get_throwable(int intIndex, int errorcode = 0); //throws pair<int,int> if intIndex not found

		typename std::map<int,t>::iterator begin();
		typename std::map<int,t>::iterator end();
		virtual ~intmap(){
		    //for( std::map<int,t>::iterator i = tmap.begin(); i != tmap.end(); i++) delete *i;
		    //Cannot use delete here because t is not guarenteed to be a poitner. if it's an int, delete [int] errors.
		}
	//protected:
		std::map<int,t> tmap;
}; //end intmap


template <class t> t intmap<t>::get(int intIndex, t default_){
	if(tmap.find(intIndex) != tmap.end()) return tmap[intIndex];
	else return default_;
} //end get


template <class t> t intmap<t>::get_throwable(int intIndex, int errorcode){
	//Attempts to return the value indexed by intIndex. 
	//if not found, throws pair<int,int>		= <intIndex, errorcode>
	if(tmap.find(intIndex) != tmap.end()) return tmap[intIndex];
	else{
		std::pair <int,int> errorpair (intIndex,errorcode);
		 throw errorpair;
	}
} //end get_throwable

template <class t> typename std::map<int,t>::iterator intmap<t>::begin(){
	return tmap.begin();
}

template <class t> typename std::map<int,t>::iterator intmap<t>::end(){
	return tmap.end();
}

//typedef intmap<bool> indexedbool;
//typedef intmap<int> indexedint;
//typedef intmap<float> indexedfloat;
typedef intmap<int> indexedstring;

class indexedint: public intmap<int>{
        public:
	void touch(int Label, int init = 0){ //if the label doesn't exist, add it to the list with value init
                if(tmap.find(Label) == tmap.end()) tmap[Label] = init;
        }
        int increment(int Label){ //increment it if it exists
                if(tmap.find(Label) != tmap.end()) return ++tmap[Label];
                else{
		    touch(Label,1);
		    return 1;
		}
        }
        int decrement(int Label){
                if(tmap.find(Label) != tmap.end()) return --tmap[Label];
                else{
		    touch(Label,-1);
		    return -1;
		}
        }
        int add(int Label, int x){ //map[label] += x
                if(tmap.find(Label) != tmap.end()) return tmap[Label]+= x;
                else{
		    touch(Label,x);
		    return x;
		}
        }
        int subtract(int Label, int x){ //map[label] -= x
                if(tmap.find(Label) != tmap.end()) return tmap[Label]-= x;
                else{
		    touch(Label,-x);
		    return -x;
		}
        }
	virtual ~indexedint(){}
};// end indexedint

class indexedbool: public intmap<bool>{
        public:
	void touch(int Label, bool init =false){ //if the label doesn't exist, add it to the list with value init
                if(tmap.find(Label) == tmap.end()) tmap[Label] = init;
        }
        bool flip(int Label){
                if(tmap.find(Label) != tmap.end()) return tmap[Label] ^= true; //flip the bool and return the new value
                else{
                        cerr<<"Warning! in indexedbool.flip received Invalid key int "<<Label<< endl;
                        return false;
                }
        }
        bool or_with(int Label, bool x){
                if(tmap.find(Label) != tmap.end()) return tmap[Label] |= x;
                else{
			touch(Label,x); //if label doesn't exist, make it, init to false, then do operation
                        return x;
                }
        }
        bool xor_with(int Label, bool x){
                if(tmap.find(Label) != tmap.end()) return tmap[Label] ^= x;
                else{
                        cerr<<"Warning! in indexedbool.xor_with received Invalid key int "<<Label<< endl;
                        return false;
                }
        }
        bool and_with(int Label, bool x){
                if(tmap.find(Label) != tmap.end()) return tmap[Label] &= x;
                else{
                        touch(Label,x); //if label doesn't exist, make it, init to true, then do operation
                        return x;
                }
        }
	virtual ~indexedbool(){}
};//en indexedbool

class indexedfloat: public intmap<float>{
        public:
	void touch(int Label, float init = 0.){ //if the label doesn't exist, add it to the list with value init
                if(tmap.find(Label) == tmap.end()) tmap[Label] = init;
        }
        float add(int Label, float x){ //map[label] += x
                if(tmap.find(Label) != tmap.end()) return tmap[Label]+= x;
                else{
                        touch(Label,x); //if label doesn't exist, make it, init to 0, then do operation. 
                        return x;
                }
        }
        float subtract(int Label, float x){ //map[label] -= x
                if(tmap.find(Label) != tmap.end()) return tmap[Label]-= x;
                else{
                        touch(Label,-x); //if label doesn't exist, make it, init to 0, then do operation. 
                        return -x;
                }
        }
        float multiply(int Label, float x){ //map[label] *= x
                if(tmap.find(Label) != tmap.end()) return tmap[Label]*= x;
                else{
                        cerr<<"Warning! in indexedfloat.multiply received Invalid key int "<<Label<< endl;
                        return 0.;
                }
        }
        float divide(int Label, float x){ //map[label] /= x
                if(tmap.find(Label) != tmap.end()) return tmap[Label]/= x;
                else{
                        cerr<<"Warning! in indexedfloat.multiply received Invalid key int "<<Label<< endl;
                        return 0.;
                }
        }
	virtual ~indexedfloat(){}
};//end indexedfloat

#endif
