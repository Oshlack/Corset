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
  string name_;
  int sample_;
  vector < Transcript * > alignments_;
 public:
  Read(){name_="";};
  Read(string name){name_=name; };
  string get_name(){return name_; } ;
  void set_sample(int sample){sample_=sample;};
  int get_sample(){return sample_;};
  int alignments(){return alignments_.size();};
  void add_alignment(Transcript * trans){
    alignments_.push_back(trans);
    if(!trans->reached_min_counts()){
      trans->add_read(this);
    };
  };

  vector < Transcript * >::iterator align_begin(){return alignments_.begin();};
  vector < Transcript * >::iterator align_end(){return alignments_.end();};

  //function to check whether an alignment has already been recorded.
  bool has(Transcript* t){return (find(align_begin(),align_end(),t) != align_end());};

  void remove(Transcript * trans){ 
     alignments_.erase(find(align_begin(),align_end(),trans));
  };

};


// ReadList is a container for a set of Reads and contains functions for
// inserting and accessing reads.
class ReadList{
 private:
    TranscriptList * transcript_list; 
    StringSet<Read> reads; 

 public:
    //we need to know all the transcripts before we can build a ReadList.
    ReadList( TranscriptList * transcripts){ transcript_list = transcripts; };
    //add a new alignment into the list
    void add_alignment(string read, string trans){
      //find the transcript id if it already exists:
      Read * r = reads.insert(read);
      Transcript * t = transcript_list->insert(trans);
      //don't try to insert an alignment if 1. it already exists or
      //2. if the transcript ID is not in the TranscriptList.
      if( !r->has(t) )
	r->add_alignment(t);
    };
    //    void print(); 
    StringSet<Read>::iterator begin(){return reads.begin();};
    StringSet<Read>::iterator end(){return reads.end();};
};


#endif
