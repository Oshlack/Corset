// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/**
 * This is a general container class which is used to store a pointer to an 
 * object along with it's associated name (or ID string). This provides
 * a convenient and fast way of storing the read information from the bam file - 
 * as a Read object pointer with it's ID string as the key. Similarly the
 * transcript information is stored as a kep, value pair of transcript name
 * and a pointer to a Transcript object. These objects are only used in the
 * first half of the program. Once the "super" clusters have been set-up, the
 * the StringSets are destroyed to conserve memory because it is not longer
 * necessary to access the objects by their IDs.
 *
 * This template is used in ReadList.h and ClusterList.h
 *
 * Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 * Last Modified: 3 May 2013
 */

#ifndef SPARSEARRAY_H
#define SPARSEARRAY_H

#include <string>
#include <cstdlib>

#ifndef UNORDEREDMAP
#include <map>
#else
#include <unordered_map>
#endif

using namespace std;

class SparseArray {
  
#ifndef UNORDEREDMAP
  typedef map< float, float > mymap;
#else
  typedef unordered_map< float, float > mymap;
#endif

 private:
  /** The underlying data structure is a map **/
  mymap sparse_map;
  int width_;

 public:  
  typedef typename mymap::iterator iterator;

  SparseArray(int width){ width_=width; } ;
  SparseArray(){};
  void set_width(int width){ width_=width; } ;

  float key(int i, int j){ return width_*j+i ; } ;

  float get(int i, int j){
    cout << "in get" <<endl;
    if(sparse_map.count(key(i,j))==0) return 0;
    return sparse_map[key(i,j)];
  };

  void set(int i, int j, float dist){
    cout << "in set" <<endl;
    if(dist==0){
      if(sparse_map.count(key(i,j))>0)
	sparse_map.erase(key(i,j)); 
      return ;
    }
    sparse_map[key(i,j)]=dist;
  } ;
  
  void erase_all( int i ){
    cout << "erase_all" <<endl;
    for(int j=0; j< width_ ; j++)
      sparse_map.erase(key(i,j));
  };
  
  float get_max_pair(int &max_i, int &max_j){
    cout << "in get_max_pair" <<endl;
    // find the shortest distance
    // stop if we find a 0 because we know this must be the shortest
    double max=-1;
    iterator it=sparse_map.begin(); 
    cout << sparse_map.size() << endl;
    for(;it != sparse_map.end(); ++it ){
      cout << it->first << " " << it->second << endl;
      if(it->second==1){
	max_j = (int) it->first / width_;
	max_i = (int) it->first % width_; 
	return 1;
      }
      if(it->second > max){
	max=it->second ;
	max_j = (int) it->first / width_;
	max_i = (int) it->first % width_;
      }
    }
    return max;
  };

};


#endif
