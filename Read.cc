// Copyright 2014 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include<Read.h>

bool Read::has_same_alignments( Read * r){
  if(r->get_sample()!=get_sample()) return false;
  if(r->alignments()!=alignments()) return false;
  return equal(align_begin(),align_end(),r->align_begin());
};

void ReadList::add_alignment(string read, string trans, int sample){
  //find the transcript id if it already exists:
  Read * r = reads_map->insert(read);
  Transcript * t = transcript_list->insert(trans);
  //don't try to insert an alignment if 1. it already exists or
  //2. if the transcript ID is not in the TranscriptList.
  if( !r->has(t) ){
    r->add_alignment(t);
    r->set_sample(sample);
  }
};

//save memory by reducing all the reads of a sample into
//a vector of "compact reads". compact reads are stored in
//regular read object, with a weight equal to the number
//of regular reads. This saves a lot of RAM if reads from multiple
//samples are processed. The map obect (StringSet) with read IDs
//is also destroyed to save memory.
void ReadList::compactify_reads(TranscriptList * trans){ 
  // first lets sort the alignments for each read

  StringSet<Read>::iterator itr=reads_map->begin();
  for(; itr!=reads_map->end(); itr++) itr->second->sort_alignments();

  // then loop over the transcripts
  TranscriptList::iterator transItr = trans->begin();
  for(;transItr!=trans->end(); transItr++){
    vector<Read*> * reads = transItr->second->get_reads();
    int reads_size=reads->size();
    for(int i=0; i < (reads_size - 1); i++){
      if(reads->at(i)->get_weight()!=0){
	for(int j=(i+1) ; j < reads_size ; j++){
	  if( reads->at(j)->get_weight()!=0 &&
	      reads->at(i)->has_same_alignments(reads->at(j))){
	    int new_weight = reads->at(i)->get_weight() + reads->at(j)->get_weight() ;
	    reads->at(i)->set_weight( new_weight ) ; //add to the weight
	    reads->at(j)->set_weight( 0 );  //set the weight of duplicates to zero
	  }
	}
      }
    }
    //now remove the reads with weight=0 from each transcript's list of reads
    for(int k=reads->size(); k > 0 ; --k){
      if(reads->at(k-1)->get_weight()==0)
	reads->erase(reads->begin()+k-1 );
    }
  }

  for(itr=reads_map->begin(); itr!=reads_map->end(); itr++){
    //remove zero weight reads at the end...
    Read * r = itr->second;
    if(r->get_weight()!=0 )
      reads_vector.push_back( r );
    else 
      delete r;
    //      if(itr->second->alignments()>1000)
    //        cout << "Found a read with "<<itr->second->alignments() << endl;
    //      else
  }
  reads_map->clear();
  delete reads_map ;
}
