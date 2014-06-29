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
 ** Modified: 3 May 2013
 **/

#ifndef READ_H
#define READ_H

#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <Transcript.h>

#include <iostream>

using namespace std;

// Read is a basic container object for a read. It can record the read ID, 
// all the transcript alignments for that read and which sample the read belongs to. 
class Read{
  //  string name_;
  vector < Transcript * > alignments_;
  unsigned char sample_;
  int weight_;
 public:
  Read(string name){ weight_=1; }; //name_=name; };
  void set_sample(int sample){sample_=sample;};
  int get_sample(){return sample_;};
  int alignments(){return alignments_.size();};
  void set_weight(int weight){ weight_ = weight; };
  int get_weight(){ return weight_; };
  void add_alignment(Transcript * trans){
    alignments_.push_back(trans);
    //    if(!trans->reached_min_counts()){
      trans->add_read(this);
      //    };
  };

  vector < Transcript * >::iterator align_begin(){return alignments_.begin();};
  vector < Transcript * >::iterator align_end(){return alignments_.end();};

  void sort_alignments(){ sort(alignments_.begin(),alignments_.end()); };

  //function to check whether an alignment has already been recorded.
  bool has(Transcript* t){return (find(align_begin(),align_end(),t) != align_end());};

  void remove(Transcript * trans){ alignments_.erase(find(align_begin(),align_end(),trans)); };

  //check if two reads have align to the same transcripts 
  //this assumed that "sort_alignments" has been called first.
  //it is used in ReadList::compactify_reads
  bool has_same_alignments( Read * r);
  
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
    void add_alignment(string read, string trans);

    //save memory by clearing the read IDs
    void compactify_reads(TranscriptList * trans); 

    vector<Read *>::iterator begin(){return reads_vector.begin();};
    vector<Read *>::iterator end(){return reads_vector.end();};

};


#endif
