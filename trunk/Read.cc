
#include<Read.h>

bool Read::has_same_alignments( Read * r){
  if(r->get_sample()!=get_sample()) return false;
  if(r->alignments()!=alignments()) return false;
  return equal(align_begin(),align_end(),r->align_begin());
};

void ReadList::add_alignment(string read, string trans){
  //find the transcript id if it already exists:
  Read * r = reads_map->insert(read);
  Transcript * t = transcript_list->insert(trans);
  //don't try to insert an alignment if 1. it already exists or
  //2. if the transcript ID is not in the TranscriptList.
  if( !r->has(t) )
    r->add_alignment(t);
};
//    void print();

void ReadList::compactify_reads(TranscriptList * trans){ //save memory by clearing the read IDs
  // first lets sort the alignments for each read
  //  int reads_before=0;
  StringSet<Read>::iterator itr=reads_map->begin();
  for(; itr!=reads_map->end(); itr++) itr->second->sort_alignments();

  // then loop over the transcripts
  TranscriptList::iterator transItr = trans->begin();
  for(;transItr!=trans->end(); transItr++){
    vector<Read*> * reads = transItr->second->get_reads();
    int reads_size=reads->size();
    // reads_before+=reads_size;
    for(int i=0; i < (reads_size - 1); i++){
      if(reads->at(i)->get_weight()!=0){
	for(int j=(i+1) ; j < reads_size ; j++){
	  if( reads->at(j)->get_weight()!=0 &&
	      reads->at(i)->has_same_alignments(reads->at(j))){
	    int new_weight = reads->at(i)->get_weight() + reads->at(j)->get_weight() ;
	    reads->at(i)->set_weight( new_weight ) ;
	    reads->at(j)->set_weight( 0 );  //set the weight of duplicates to zero
	  }
	}
      }
    }
    //now remove the read from each transcripts list of reads
    for(int k=reads->size(); k > 0 ; --k){
      if(reads->at(k-1)->get_weight()==0)
	reads->erase(reads->begin()+k-1 );
    }
  }

  /**  int reads_after=0;
  for(transItr=trans->begin();transItr!=trans->end(); transItr++){
    vector<Read*> * reads = transItr->second->get_reads();
    reads_after+=reads->size();
  }
  cout << reads_before<< " " << reads_after << " ";**/

  int removed2=0;
  for(itr=reads_map->begin(); itr!=reads_map->end(); itr++){
    //remove zero weight reads at the end...
    Read * r = itr->second;
    if(r->get_weight()!=0 )
      reads_vector.push_back( r );
    else{
      removed2++;
      delete r;
    }
    //      if(itr->second->alignments()>1000)
    //        cout << "Found a read with "<<itr->second->alignments() << endl;
    //      else
  }
  //  cout << removed2 << " " << reads_map->size() << " " << reads_vector.size() << endl;
  reads_map->clear();
  delete reads_map ;
}
