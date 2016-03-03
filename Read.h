// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/***
 ** This file contains the declarations and definitions for two classes:
 ** Read and ReadList. 
 **
 ** Author: Nadia Davidson
 ** Modified: 11th July 2014
 **/

#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <Transcript.h>

using namespace std;

// Read is a basic container object for a read. It can record the read ID, 
// all the transcript alignments for that read and which sample the read belongs to. 
class Read{
  vector < Transcript * > alignments_;
  unsigned char sample_;
  int weight_;
  uintptr_t trans_hash;

 public:
  Read(){ weight_=1; }; 
  Read(string name){ weight_=1; }; //dummy for StringSet
  void set_sample(int sample){sample_=sample;};
  int get_sample(){return sample_;};
  int alignments(){return alignments_.size();};
  void set_weight(int weight){ weight_ = weight; };
  int get_weight(){ return weight_; };
  void add_alignment(Transcript * trans){
    alignments_.push_back(trans);
    trans->add_read(this);
  };

  vector < Transcript * >::iterator align_begin(){return alignments_.begin();};
  vector < Transcript * >::iterator align_end(){return alignments_.end();};

  void sort_alignments(){ sort(alignments_.begin(),alignments_.end()); };

  //function to check whether an alignment has already been recorded.
  bool has(Transcript* t){return (find(align_begin(),align_end(),t) != align_end());};

  void remove(Transcript * trans){ alignments_.erase(find(align_begin(),align_end(),trans)); };

  //check if two reads align to all the same transcripts 
  //this assumed that "sort_alignments" has been called first.
  //it is used in ReadList::compactify_reads
  bool has_same_alignments( Read * r);

  // return the hashed values of all transcript pointers that this read aligns to
  // used to quickly compare alignments in the compactify_reads function
  uintptr_t get_trans_hash(){
    return trans_hash ;
  }

  // sets the hashed values of all transcript pointers that this read aligns to using xor
  // used to quickly compare alignments in the compactify_reads function
  uintptr_t set_trans_hash(){
    trans_hash = 0;
    vector < Transcript * >::iterator t_itr=align_begin();
    for( ; t_itr!=align_end() ; t_itr++){
      // this has function was derived from boost::hash_combine()
      trans_hash ^= hash<uintptr_t>()((uintptr_t)(*t_itr)) + 0x9e3779b9 + (trans_hash<<6) + (trans_hash>>2);
    }
  }

  void print_alignments(){
    vector < Transcript * >::iterator t_itr=align_begin();
    for( ; t_itr!=align_end() ; t_itr++)
      cout << (*t_itr)->get_name() << " ";
    cout << endl;
  }

};


// ReadList is a container for a set of Reads and contains functions for
// inserting and accessing reads.
class ReadList{
 private:
    TranscriptList * transcript_list; 
    StringSet<Read> * reads_map; //warning this is deleted after the reads are read
    vector<Read *> reads_vector;

 public:
    //we need to know all the transcripts before we can build a ReadList.
    ReadList( TranscriptList * transcripts){ transcript_list = transcripts; reads_map = new StringSet<Read> ; };

    //add a new alignment into the list
    //this one is used when reading bam files
    void add_alignment(string read, string trans, int sample);
    //this one is used when reading corset read summary files
    void add_alignment(vector<string> trans, int sample, int weight);

    //saves memory by reducing each read into a set of
    //"compact reads" with a weight.
    //Also the map object is clear and the reads are
    //stored as a vector instead. Read IDs are cleared.
    void compactify_reads(TranscriptList * trans, string outputReadsName=""); 

    vector<Read *>::iterator begin(){return reads_vector.begin();};
    vector<Read *>::iterator end(){return reads_vector.end();};

    void write(string outputReadsName);

};


#endif
