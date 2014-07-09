// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <MakeClusters.h>

using namespace std;


pair< Transcript * const, Cluster * > * MakeClusters::getMapElement(Transcript * trans){
  //  StringSet<Cluster>::iterator it=transMap.find(trans);
  map< Transcript*, Cluster*>::iterator it=transMap.find(trans);
  if(it!=transMap.end())
    return &(*it);
  Cluster * clust = new Cluster();
  clust->add_tran(trans);
  clusterList.push_back(clust);
  pair< map< Transcript *, Cluster*>::iterator,bool> newTrans=transMap.insert( 
				       		      pair< Transcript *, Cluster* >(trans,clust));
  return &(*(newTrans.first));
};

void MakeClusters::checkAgainstCurrentCluster(Transcript * trans){
  //compare the current transcript's cluster to the "current" cluster.
  Cluster * this_cluster=getMapElement(trans)->second;
  //"this_cluster" should match "current_cluster"
  if(this_cluster==current_cluster) return;
  
  //otherwise, change the cluster of all transcripts in "this_cluster" to "current_cluster"
  for(int i=0; i<this_cluster->n_trans(); i++){
    transMap[this_cluster->get_tran(i)]=current_cluster;
    current_cluster->add_tran(this_cluster->get_tran(i));
  }
  //and copy all the reads across too
  for(int i=0; i<this_cluster->n_reads(); i++){
    current_cluster->add_read(this_cluster->get_read(i));
  }
  //clean up
  clusterList.erase(find(clusterList.begin(),clusterList.end(),this_cluster));
  delete this_cluster;
};

void MakeClusters::makeSuperClusters(vector<ReadList*> & readLists){
  //loop through each read: 
  int i=0;
  for(int sample=0; sample<readLists.size(); sample++){
    ReadList reads=*(readLists.at(sample));

    vector< Read* >::iterator rIt;
    for( rIt=reads.begin() ; rIt != reads.end(); rIt++ ){
      Read * r = *rIt ; //->second;
      //r->set_sample(sample);
      vector< Transcript * >::iterator tIt ;
      for( tIt=r->align_begin() ; tIt != r->align_end(); tIt++ ){
	Transcript * trans=*tIt;
	if(tIt==r->align_begin()){
	  setCurrentCluster(trans);
	  current_cluster->add_read(r);
	}
	else
	  checkAgainstCurrentCluster(trans);
      }
      if(i % 100000 == 0)
	cout << float(i)/float(1000000) << " million compact reads read" <<endl;
      i++;
    }
    delete readLists.at(sample);
  }
  //clean up 
  readLists.clear();
  transMap.clear();
}

void MakeClusters::processSuperClusters(map<float,string> & distance_thresholds, vector<int> & groups){
  //now do the hierarchical clustering...
  cout << "Starting hierarchial clustering..." << endl;
  for(int i=0; clusterList.size()>0; i++){
    if(i % 1000 == 0)
      cout << float(i)/float(1000) << " thousand clusters done" <<endl;
    Cluster * back = clusterList.back();
    back->set_id(i);
    back->set_sample_groups(groups);
    back->cluster(distance_thresholds);	 
    //    back->print_alignments();
    clusterList.pop_back();
    delete back;
  }
}

MakeClusters::MakeClusters(vector<ReadList*> & readLists, 
			   map<float,string> & distance_thresholds, 
			   vector<int> & groups){

  //stage 1: process all the reads and looked for shared hits.
  //groups all transcripts which share at least one read
  makeSuperClusters(readLists);
  //stage 2: loop over each of the newly created groups (super clusters)
  //and perform the hierarchical clustering
  processSuperClusters(distance_thresholds,groups);

};


