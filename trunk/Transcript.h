// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/**
 * A class to store the names of transcripts.
 * This allows transcripts to be idenfitied by a pointer to their "pair" 
 * rather than using the full name (which is a string). This makes the code
 * run faster and require less memory.
 *
 * Author: Nadia Davidson
 * Modified: 3 May 3013
 */

#ifndef TRANSCRIPT_H
#define TRANSCRIPT_H

#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <algorithm>
#include <StringSet.h>

using namespace std;

class Read;
//void Read::remove(Transcript *);

class Transcript{
  string name_;
  int pos_; // used later by Cluster
  vector<Read*> reads_; //temporary vector so we can quickly remove the alignments
                 //for transcripts with less then min_count hits.
  bool reached_min_counts_;

 public:
  Transcript(){name_="";};
  Transcript(string name);

  string get_name(){return name_; } ;
  void pos(int position){pos_=position;};
  int pos(){return pos_;};
  void add_read( Read * read );
  bool reached_min_counts(){ return reached_min_counts_; };
  void remove(); //remove myself from the reads lists .. 

  static int samples;
  static int min_counts;
};

typedef StringSet<Transcript> TranscriptList;

#endif
