// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

//#include <string>
//#include <vector>
//#include <sstream>
//#include <cstdlib>
//#include <algorithm>
//#include <StringSet.h>
#include <Transcript.h>
#include <Read.h>

using namespace std;

// initialise the number of samples to 0.
// this will be set later on.
int Transcript::samples=0;
int Transcript::groups=0;
int Transcript::min_counts=10;

void Transcript::remove(){
    for(int r=0; r<reads_.size() ; r++)
      reads_.at(r)->remove(this);
};

Transcript::Transcript(string name){
    name_=name;
    //    reached_min_counts_=false;
};

void Transcript::add_read( Read * read ){ 
  //  if(reads_.size()>=min_counts)//{
    //    reads_.clear();
    //   reached_min_counts_=true;
    //} else {
    reads_.push_back(read) ; 
    //  }
};

bool Transcript::reached_min_counts(){
  int counts=0;
  for(int i=0; i<reads_.size(); i++){
    counts+=reads_.at(i)->get_weight();
    cout << min_counts << endl;
    if(counts >= min_counts) return true;
  }
  return false;
}

