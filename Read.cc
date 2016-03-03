// Copyright 2014 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include<Read.h>

bool Read::has_same_alignments( Read * r){

  if(r->get_trans_hash()!=get_trans_hash()) return false;
  return(r->get_sample()==get_sample());

  //   return (r->get_trans_hash()==get_trans_hash()); //equal(align_begin(),align_end(),r->align_begin());
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

//This is the method called when we read corset files as input.
//reads are already 'compact', so we can add to the read_vector instead of the map.
void ReadList::add_alignment(vector<string> trans_names, int sample, int weight){
  //make the new read
  Read * r = new Read();
  reads_vector.push_back(r);
  r->set_sample(sample);
  r->set_weight(weight); //this read represents mutliple reads in the original bam
  
  //loop over all the transcripts that this read aligns to
  vector<string>::iterator itrTrans = trans_names.begin();
  for(; itrTrans!=trans_names.end(); itrTrans++){
    //find the transcript object with the name  
    Transcript * t = transcript_list->insert(*itrTrans);
    r->add_alignment(t);
  }
};

//save memory by reducing all the reads of a sample into
//a vector of "compact reads". compact reads are stored in
//regular read object, with a weight equal to the number
//of regular reads. This saves a lot of RAM if reads from multiple
//samples are processed. The map obect (StringSet) with read IDs
//is also destroyed to save memory.
void ReadList::compactify_reads(TranscriptList * trans, string outputReadsName){ 

  // first lets sort the alignments for each read
  // and calculate a hash value to be used when comparing alignments
  StringSet<Read>::iterator itr=reads_map->begin();
  for(; itr!=reads_map->end(); itr++){ 
    itr->second->sort_alignments(); 
    itr->second->set_trans_hash();
  }

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
  }
  reads_map->clear();
  delete reads_map ;
}

void ReadList::write(string outputReadsName){
  //now output the read mapping to file if requested
  ofstream readFile;
  readFile.open(outputReadsName);
  vector<Read *>::iterator read;
  for(read=begin(); read!=end(); read++){
    readFile << (*read)->get_weight() ;
    vector < Transcript * >::iterator trans;
    for(trans=(*read)->align_begin(); trans!=(*read)->align_end(); trans++)
      readFile << "\t" << (*trans)->get_name();
    readFile << endl;
  }
  readFile.close();
  cout << "Done writing "<< outputReadsName << endl;
}
