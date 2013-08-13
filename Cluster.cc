// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#include <Cluster.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

//can be set by user (see corset.cc)
double Cluster::D_cut=0;
string Cluster::file_prefix="";

//constants
const string Cluster::file_counts="counts";
const string Cluster::file_clusters="clusters";
const string Cluster::file_ext=".txt";
const string Cluster::cluster_id_prefix="Cluster-";
const string Cluster::cluster_id_joiner=".";
const string Cluster::cluster_id_prefix_no_reads="NoReadsCluster-";

//getter/setter methods
void Cluster::add_tran(Transcript * trans){ cluster_.push_back(trans); };
void Cluster::add_read(Read * read){ read_.push_back(read); };
Transcript * Cluster::get_tran(int i){ return cluster_.at(i); };
Read * Cluster::get_read(int i){ return read_.at(i); };
int Cluster::n_trans(){ return cluster_.size(); };
int Cluster::n_reads(){ return read_.size(); };

void print(vector<int> x){
  vector<int>::iterator it=x.begin();
  for(;it!=x.end(); it++)
    cout << *it << " ";
  cout << endl;
}



double Cluster::get_dist(int i, int j ){

  //first, work out how many shared reads we have
  int shared_reads = 0; 
  double total_reads_i = 0;
  double total_reads_j = 0;

  vector<int> sample_shared_read;
  sample_shared_read.resize(Transcript::samples,0);

  for(int s=0; s < Transcript::samples; s++){
    vector<int> & readsG1=read_groups.at(s).at(i);
    vector<int> & readsG2=read_groups.at(s).at(j);

    vector<int>::iterator it1=readsG1.begin();
    vector<int>::iterator it2=readsG2.begin();
    
    while(it1!=readsG1.end() && it2!=readsG2.end()){
      if(*it1 < *it2) ++it1;
      else if(*it1 > *it2) ++it2;
      else {
	sample_shared_read.at(s)++;
	it1++;
	it2++;
      }
    }
    shared_reads+=sample_shared_read.at(s);

    total_reads_i+=readsG1.size();
    total_reads_j+=readsG2.size();
  }

  double shared_reads_distance = shared_reads / min( total_reads_i, total_reads_j ) ;
  
  // there is no need to look for differences in expression
  // if the distance is already at the maximum.
  if(shared_reads_distance==0)
    return shared_reads_distance; 
  
  //now look to see if there is different differential expression between
  //sample which would suggest that the two transcripts do not belong to the
  //same gene.
  int ngroups=Transcript::groups;
  vector<double> x (ngroups,1); //start with 1 count to avoid dividing by 0.
  vector<double> y (ngroups,1);
  for(int s=0; s < Transcript::samples; s++){
    int g=sample_groups.at(s);
    int shared=0.5*sample_shared_read.at(s); //subtract half the shared reads
    x.at(g)+=read_groups.at(s).at(i).size()-shared; //for each group so that
    y.at(g)+=read_groups.at(s).at(j).size()-shared; //we're not double counting.
  }

  //pretend that there's only one sample per group for the moment
  //calculate the likelihood when the ratios are different:
  double non_null=0;
  double x_all=0;
  double y_all=0;
  for(int g=0; g < x.size(); g++){
    double r=y.at(g)/x.at(g);
    non_null += y.at(g)*log(r*x.at(g)) - r*x.at(g);
    non_null += x.at(g)*log(x.at(g)) - x.at(g);
    x_all+=x.at(g);
    y_all+=y.at(g);
  }

  //when the ratios are the same.. 
  double null=0;
  double r_all=y_all/x_all;
  for(int g=0; g < x.size(); g++){
    double mean_x = (x.at(g)+y.at(g))/((double)(1+r_all));
    null += y.at(g)*log(r_all*mean_x) - r_all*mean_x;
    null += x.at(g)*log(mean_x) - mean_x;
  } 
  
  double D=2*non_null-2*null;
  if(D>D_cut) //10 is approx. p<0.001
    return 0; //set the distance to the max.

  return shared_reads_distance;
}

void Cluster::merge(int i, int j, int n){

  //update the group information
  int l=n-1;
  groups.at(i).insert(groups.at(i).end(), groups.at(j).begin(), groups.at(j).end());
  groups.at(j).swap(groups.at(l));
  groups.pop_back();

  //update the read alignments
  for(int s=0; s< Transcript::samples ; s++){
    //make a temporary set to store the result. (will be ordered and unique)
    vector<int> temp;
    vector<int>::iterator it1 = read_groups.at(s).at(i).begin();
    vector<int>::iterator it2 = read_groups.at(s).at(j).begin();
    vector<int>::iterator it1_end = read_groups.at(s).at(i).end();
    vector<int>::iterator it2_end = read_groups.at(s).at(j).end();
    while(it1!=it1_end && it2!=it2_end){
      if(*it1 < *it2) { temp.push_back(*it1); ++it1; }
      else if(*it1 > *it2){temp.push_back(*it2); ++it2; }
      else {
	temp.push_back(*it1);
	it1++;
	it2++;
      }
    }
    while(it1!=it1_end){
      temp.push_back(*it1);
      it1++;
    }
    while(it2!=it2_end){
      temp.push_back(*it2);
      it2++;
    }
    read_groups.at(s).at(i).clear();
    read_groups.at(s).at(i).insert(read_groups.at(s).at(i).begin(),temp.begin(),temp.end());
    read_groups.at(s).at(j).swap(read_groups.at(s).at(l));
    read_groups.at(s).pop_back(); 
  }

  //now update the distance matrix
  // i <-> i & j
  // j <-> l
  // i is always > j
  int m=i;
  if(m==l) m=j; //m holds the index of the new merged read_groups
  // if(i==l) i=j;

  //update all the ith row and column to contain the ij pairs total counts
  for(int k=0; k<i; k++){
    if(dist[i][k]==0 && ((k<j && dist[j][k]==0)||(j<k && dist[k][j]==0))){
      dist[i][k]=0; //this stuff speeds up the algorithm
    }
    else if( (dist[i][k]==1) && ((j<k && dist[k][j]==1) || (j>k && dist[j][k]==1))){
      dist[i][k]=1;
    }
    else
      dist[i][k]=get_dist(m,k);
  }
  for(int k=i+1; k<l; k++){
    if(dist[k][i]==0 && dist[k][j]==0) //j is always less than i
      dist[k][i]=0;
    else if (dist[k][i]==1 && dist[k][j]==1)
      dist[k][i]=1;
    else
      dist[k][i]=get_dist(m,k);
  }
  
  //move the last row to j to reduce the size of the matrix
  //if(i!=j){ //but only if we haven't already done this.
    for(int k=0; k<j; k++)
      dist[j][k]=dist[l][k];
    for(int k=j+1; k<n; k++)
      dist[k][j]=dist[l][k];
    //  }
}

double Cluster::find_next_pair(int n, int &max_i, int &max_j){

  // find the shortest distance
  // stop if we find a 0 because we know this must be the shortest
  double max=-1;
  for(int i=1; i<n; i++){
      for(int j=i-1; j>=0; j--){
	if(dist[i][j]==1){
	  max_i=i;
	  max_j=j;
	  return 1;
	}
	if(dist[i][j]>max){
	  max=dist[i][j];
	  max_i=i;
	  max_j=j;
	}
      }
  }
  return max;
}

vector<int> Cluster::get_counts(int s){
  // we need to reorder the sets by reads.
  // start by looping over the reads 
  vector< set<int> > read_groups_inv;
  vector<int> counts;
  
  read_groups_inv.resize(n_reads());
  counts.resize(read_groups.at(s).size());

  for(int g=0; g < read_groups.at(s).size() ; g++){
    vector<int>::iterator itrRG = read_groups.at(s).at(g).begin();
    for( ; itrRG != read_groups.at(s).at(g).end() ; itrRG++){
      read_groups_inv.at(*itrRG).insert(g);
    }
  }
  
  //now loop over the reads
  for(int rd=0; rd < read_groups_inv.size() ; rd++){
    int n_align =  read_groups_inv.at(rd).size() ;
    if(n_align>0){
      int which_group_pos = rand() % n_align ;
      set<int>::iterator itr=read_groups_inv.at(rd).begin();
      advance(itr,which_group_pos); 
      counts.at(*itr)++; 
    }
  } //randomly choose a transcript for the moment

  return counts;
}

//vector<Cluster*> Cluster::process(double h, string filename){
void Cluster::cluster( map<double,string> & thresholds){

  //nothing to do if we only have one transcript in the cluster
  if(n_trans()>1000)
    cout << "cluster with "<<n_trans() << " transcripts.. this might take a while"<<endl;

  //prepare a few data objects
  initialise_matrix();

  //start clustering
  map< double,string >::iterator itr_d_cut=thresholds.begin(); //which element of d_cut are we pointing to.
  for(int n=n_trans(); n>1 ; n-- ){
    int i=0;
    int j=0;
    double distance = 1 - find_next_pair(n,i,j);
    if( n>1000 & n % 200==0 )
      cout << "down to "<<n<< " clusters. dist=" << distance << endl;
    while( distance > itr_d_cut->first & itr_d_cut != thresholds.end() ){ 
      output_clusters( itr_d_cut->second ); //print the results
      itr_d_cut++;
    }
    if(distance==1) break;
    
    //or continue on...
    merge(i,j,n);
  }
  //report the final clustering.
  while( itr_d_cut != thresholds.end() ){ 
    output_clusters(itr_d_cut->second);
    itr_d_cut++;
  }
}


void Cluster::output_clusters(string threshold){
      //do the counts
      vector< vector<int> > counts;
      for(int s=0; s<Transcript::samples; s++)
	counts.push_back(get_counts(s));

      ofstream countsFile;
      countsFile.open((file_prefix+string(file_counts+threshold+file_ext)).c_str(),ios_base::app);
      for(int g=0; g < counts.at(0).size(); g++){
	countsFile << cluster_id_prefix << get_id() << cluster_id_joiner << g ;
	for(int s=0; s < Transcript::samples; s++)
	  countsFile << "\t"<< counts.at(s).at(g) ;
	countsFile << endl;
      }
      countsFile.close();

      //do the clusters
      ofstream clusterFile;
      clusterFile.open((file_prefix+string(file_clusters+threshold+file_ext)).c_str(),ios_base::app);
      for(int g=0; g < groups.size(); g++)
	for(int t=0; t < groups.at(g).size(); t++){
	  clusterFile << get_tran(groups.at(g).at(t))->get_name() 
		      << "\t" << cluster_id_prefix << get_id() << cluster_id_joiner << g << endl;
	}
      clusterFile.close();
}


void Cluster::clean_up(){
  for(int i=1; i<n_trans(); i++){
    delete dist[i];
  }
  delete dist;
}

void Cluster::initialise_matrix(){

  //set the position variable for the transcripts
  //loop over transcripts:
  for(int i=0; i<n_trans(); i++)
    get_tran(i)->pos(i);

       
  //allocate space for the distance matrix
  dist = new double*[n_trans()];
  for(int i=1; i<n_trans(); i++){
    dist[i]=new double[i+1];
    for(int j=0; j<i; j++){
      dist[i][j] = 0;
    }
  }
  
  read_groups.resize(Transcript::samples);
  for(int s=0; s < Transcript::samples ; s++)
    read_groups.at(s).resize(n_trans());

  vector< Transcript *>::iterator t1;
  vector< Transcript *>::iterator t2;
  for(int r=0; r<n_reads(); r++){
    Read * read = get_read(r);
    int alignments = read->alignments();
    int sample = read->get_sample();

    for(t1=read->align_begin(); t1!=read->align_end(); t1++){
      int i=(*t1)->pos();
      for(t2=read->align_begin(); t2!=t1; t2++){ //flag the elements in the distance matrix
	int j=(*t2)->pos();                      //which need to be calculated properly.
	if(j<i)
	  dist[i][j]=1;
	else 
	  dist[j][i]=1;
      }
      read_groups.at(sample).at(i).push_back(r); 
    }
  }
  
  groups.resize(n_trans());
  for(int n=0; n < n_trans(); n++)
    groups.at(n).push_back(n);
  
  //now set the distances
  for(int i=1; i<n_trans(); i++){
    for(int j=0; j<i; j++){
      if(dist[i][j]==1)
	dist[i][j] = get_dist(i,j);
    }
  }
};


