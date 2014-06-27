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

#ifndef STRINGSET_H
#define STRINGSET_H

#include <string>
#include <cstdlib>

#ifndef UNORDEREDMAP
#include <map>
#else
#include <unordered_map>
#endif

using namespace std;

template <class T > class StringSet {
  
#ifndef UNORDEREDMAP
  typedef map< string, T * > mymap;
#else
  typedef unordered_map< string, T * > mymap;
#endif

 private:
  /** The underlying data structure is a map **/
  mymap set_map;
  
 public:  
  typedef typename mymap::iterator iterator;

  /** Using the key string, look up the pointer to the object. 
      For example use the read ID to get a Read * or a 
      transcript ID to get a Transcript *.
   **/
  T * find(string & name){
    iterator it;
    it = set_map.find(name);
    if(it!=set_map.end())
      return it->second;
    return NULL;
  };

  /** Create a new object, insert it and its ID string into the map,
      then return the pointer to the object. Check to see if the name
      already exisits and if it does, just return the already exisiting
      object pointer.
   **/
  T * insert(string name){
    T * look_up = find(name); 
    if(look_up!=NULL)
      return look_up;
    else{
      T * pt = new T(name);
      set_map[name]=pt;
      return pt;
    }
  };
  
  iterator begin(){
    return set_map.begin();
  };

  iterator end(){
    return set_map.end();
  };

  void clear(){ set_map.clear(); }
};


#endif
