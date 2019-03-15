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
float Cluster::D_cut=0;
string Cluster::file_prefix="";

//constants
const string Cluster::file_counts="counts";
const string Cluster::file_clusters="clusters";
const string Cluster::file_ext=".txt";
const string Cluster::cluster_id_prefix="Cluster-";
const string Cluster::cluster_id_joiner=".";
const string Cluster::cluster_id_prefix_no_reads="NoReadsCluster-";
//const unsigned char Cluster::distMAX=255;

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


float Cluster::get_dist(int i, int j ){

  //first, work out how many shared reads we have
  int shared_reads = 0; 
  float total_reads_i = 0;
  float total_reads_j = 0;

  vector<int> sample_shared_read( Transcript::samples,0 );

  for(int s=0; s < Transcript::samples; s++){
    vector<int> & readsG1=read_groups.at(s).at(i);
    vector<int> & readsG2=read_groups.at(s).at(j);

    vector<int>::iterator it1=readsG1.begin();
    vector<int>::iterator it2=readsG2.begin();
    
    while(it1!=readsG1.end() && it2!=readsG2.end()){
      if(*it1 < *it2){
	++it1;
      }
      else if(*it1 > *it2){ 
	++it2;
      }
      else {
	sample_shared_read[s]+=get_read(*it1)->get_weight();
	it1++;
	it2++;
      }
    }
    shared_reads+=sample_shared_read[s];
    total_reads_i+=read_group_sizes[i][s];
    total_reads_j+=read_group_sizes[j][s];
  }

  float shared_reads_distance = shared_reads / min( total_reads_i, total_reads_j ) ;

  // there is no need to look for differences in expression
  // if the distance is already at the maximum.
  if(shared_reads_distance==0)
    return shared_reads_distance; 
  
  //now look to see if there is different differential expression between
  //samples which would suggest that the two transcripts do not belong to the
  //same gene.
  int ngroups=Transcript::groups;
  vector<float> x (ngroups,1); //start with 1 count to avoid dividing by 0.
  vector<float> y (ngroups,1);
  for(int s=0; s < Transcript::samples; s++){
    int g=sample_groups[s];
    int shared=0.5*sample_shared_read[s]; //subtract half the shared reads
    x[g]+=((read_group_sizes[i][s])-shared); //for each group so that
    y[g]+=((read_group_sizes[j][s])-shared); //we're not double counting.
  }

  //pretend that there's only one sample per group for the moment
  //calculate the likelihood when the ratios are different:
  float non_null=0;
  float x_all=0;
  float y_all=0;
  for(int g=0; g < x.size(); g++){
    float r=y[g]/x[g];
    non_null += y[g]*log(r*x[g]) - r*x[g];
    non_null += x[g]*log(x[g]) - x[g];
    x_all+=x[g];
    y_all+=y[g];
  }

  //when the ratios are the same.. 
  float null=0;
  float r_all=y_all/x_all;
  for(int g=0; g < x.size(); g++){
    float mean_x = (x[g]+y[g])/((float)(1+r_all));
    null += y[g]*log(r_all*mean_x) - r_all*mean_x;
    null += x[g]*log(mean_x) - mean_x;
  } 
  
  float D=2*non_null-2*null;

  if(D>D_cut) //10 is approx. p<0.001?
    return 0; //set the distance to the max.

  return shared_reads_distance;
}

void Cluster::merge(const int i, const int j){

  //update the group information
  groups.at(i).insert(groups.at(i).end(), groups.at(j).begin(), groups.at(j).end());
  groups.at(j).clear();

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
      else if(*it1 > *it2){ temp.push_back(*it2); ++it2; }
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
    read_groups.at(s).at(j).clear();
  }

  //recalculate the number of reads for cluster, i, sample, s. 
  for(int s=0; s< Transcript::samples ; s++){
    read_group_sizes.at(i).at(s)=0;
    for(int r=0; r< read_groups.at(s).at(i).size(); r++)
      read_group_sizes.at(i).at(s)+=get_read(read_groups.at(s).at(i).at(r))->get_weight();
  }
  read_group_sizes.at(j).clear();

  //now update the distance matrix
  // i <-> i & j
  // j <-> clear
  // i is always > j

  //update all the ith row and column to contain the ij pairs total counts
  for(int k=0; k<i; k++){
    if(dist.no_link(i,k) && ((k<j && dist.no_link(j,k)) || (j<k && dist.no_link(k,j)) )){
      continue;
    }
    else if( (dist.get(i,k)==UCHAR_MAX) && ((j<k && dist.get(k,j)==UCHAR_MAX) || (j>k && dist.get(j,k)==UCHAR_MAX))){
      continue; 
    }
    else
      dist.set(i,k,get_dist(i,k)*UCHAR_MAX);
  }
  for(int k=i+1; k<n_trans(); k++){
    if(dist.no_link(k,i) && dist.no_link(k,j)) //j is always less than i
      continue;
    else if (dist.get(k,i)==UCHAR_MAX && dist.get(k,j)==UCHAR_MAX)
      continue; 
    else
      dist.set(k,i,get_dist(i,k)*UCHAR_MAX);
  }

  //move the last row to j to reduce the size of the matrix
  for(int k=0; k<j; k++)
      dist.remove(j,k); 
  for(int k=j+1; k<n_trans(); k++)
      dist.remove(k,j); 

}


unsigned char Cluster::find_next_pair(int &max_i, int &max_j){

  dist_iterator ditr=dist.begin();
  unsigned char max_value=0;
  //loop through non-zero distances and search for the maximum
  for(; ditr!=dist.end(); ditr++){
    if(ditr->second > max_value){
      max_i = dist.get_i(ditr->first);
      max_j = dist.get_j(ditr->first);
      max_value = ditr->second;
    }
    if(max_value==UCHAR_MAX) //maximum value, no need to look any further
      break;
  }
  return(max_value);

}

//returns a vector of counts (one for each cluster), for a given sample
vector<int> Cluster::get_counts(int s){
  // we need to reorder the sets by reads.
  // start by looping over the reads 
  vector< set<int> > read_groups_inv;
  vector<int> counts;
  
  read_groups_inv.resize(n_reads());
  counts.resize(read_groups.at(s).size());
  //s=sample, g=cluster
  for(int g=0; g < read_groups.at(s).size() ; g++){
    vector<int>::iterator itrRG = read_groups.at(s).at(g).begin();
    for( ; itrRG != read_groups.at(s).at(g).end() ; itrRG++){
      read_groups_inv.at(*itrRG).insert(g);
    }
  } //we now have a 2D array - dim 1: reads, dim 2: cluster
  
  //now loop over the reads (rd=read)
  for(int rd=0; rd < read_groups_inv.size() ; rd++){
    int n_align =  read_groups_inv.at(rd).size() ;
    if(n_align>0){
       int weight = get_read(rd)->get_weight();
       for( ; weight!=0 ; weight--){
          int which_group_pos = rand() % n_align ;
          set<int>::iterator itr=read_groups_inv.at(rd).begin();
          advance(itr,which_group_pos); 
          counts.at(*itr)++; 
       }
    }
  } //randomly choose a cluster is a read maps to more than one cluster

  return counts;
}

//vector<Cluster*> Cluster::process(float h, string filename){
void Cluster::cluster( map<float,string> & thresholds){

  //nothing to do if we only have one transcript in the cluster
  if(n_trans()>1000)
    cout << "cluster with "<<n_trans() << " transcripts.. this might take a while"<<endl;

  //prepare a few data objects
  initialise_matrix();

  //start clustering
  map< float,string >::iterator itr_d_cut=thresholds.begin(); //which element of d_cut are we pointing to.
  for(int n=n_trans(); n>1 ; n-- ){
    int i=0;
    int j=0;
    float distance = 1 - ((float)find_next_pair(i,j))/UCHAR_MAX;
    if( n>1000 & n % 200==0 )
      cout << "down to "<<n<< " clusters. dist=" << distance << endl;
    while( distance > itr_d_cut->first & itr_d_cut != thresholds.end() ){ 
      output_clusters( itr_d_cut->second ); //print the results
      itr_d_cut++;
    }
    //stop when we go over the final threshold, or when the distance is the max
    if(itr_d_cut == thresholds.end() || distance==1 ) break;
    
    //or continue on...
    merge(i,j);
  }

  //report the final clustering.
  while( itr_d_cut != thresholds.end() ){ 
    output_clusters(itr_d_cut->second);
    itr_d_cut++;
  }

  //  clean_up();
}


void Cluster::output_clusters(string threshold){
      //do the counts
      vector< vector<int> > counts;
      for(int s=0; s<Transcript::samples; s++)
	counts.push_back(get_counts(s));

      ofstream countsFile;
      countsFile.open((file_prefix+string(file_counts+threshold+file_ext)).c_str(),ios_base::app);
      int id=0;
      for(int g=0; g < counts.at(0).size(); g++){
	if(groups.at(g).size()==0) continue; //if the cluster has no transcript (no longer in use) ignore
	countsFile << cluster_id_prefix << get_id() << cluster_id_joiner << id ;
	for(int s=0; s < Transcript::samples; s++)
	  countsFile << "\t"<< counts.at(s).at(g) ;
	countsFile << endl;
	id++;
      }
      countsFile.close();

      //do the clusters
      ofstream clusterFile;
      clusterFile.open((file_prefix+string(file_clusters+threshold+file_ext)).c_str(),ios_base::app);
      id=0; //reset
      for(int g=0; g < groups.size(); g++){
	for(int t=0; t < groups.at(g).size(); t++){
	  clusterFile << get_tran(groups.at(g).at(t))->get_name() 
		      << "\t" << cluster_id_prefix << get_id() << cluster_id_joiner << id << endl;
	}
	if(groups.at(g).size()!=0)
	  id++;
      }
      clusterFile.close();
}


void Cluster::initialise_matrix(){

  //set the position variable for the transcripts
  //loop over transcripts:
  for(int i=0; i<n_trans(); i++)
    get_tran(i)->pos(i);

  //set dimension of distance matrix
  dist.set_size(n_trans());

  read_groups.resize(Transcript::samples);
  for(int s=0; s < Transcript::samples ; s++)
    read_groups.at(s).resize(n_trans());
  read_group_sizes.resize(n_trans());
  for(int t=0; t < n_trans() ; t++)
    read_group_sizes.at(t).resize(Transcript::samples);

  vector< Transcript *>::iterator t1;
  vector< Transcript *>::iterator t2;
  vector< pair<int, int>> non_zero_distance_index; 
  for(int r=0; r<n_reads(); r++){
    Read * read = get_read(r);
    int alignments = read->alignments();
    int sample = read->get_sample();
    for(t1=read->align_begin(); t1!=read->align_end(); t1++){
      int i=(*t1)->pos();
      for(t2=read->align_begin(); t2!=t1; t2++){ //flag the elements in the distance matrix
	int j=(*t2)->pos();                      //which need to be calculated properly.
	if(j<i)
	  non_zero_distance_index.push_back(make_pair(i,j));
	else 
	  non_zero_distance_index.push_back(make_pair(j,i));
      }
      read_groups.at(sample).at(i).push_back(r); 
      read_group_sizes.at(i).at(sample)+=(read->get_weight());  
    }
  }
  
  groups.resize(n_trans());
  for(int n=0; n < n_trans(); n++)
    groups.at(n).push_back(n);

  //now set the distances
  for(int p=0; p< non_zero_distance_index.size(); p++){
    int i=non_zero_distance_index.at(p).first;
    int j=non_zero_distance_index.at(p).second;
    dist.set(i,j,get_dist(i,j)*UCHAR_MAX);
  }

};

void Cluster::print_alignments(){
  
  if(n_trans()<10000) return ;

  for(int r=0; r<n_reads(); r++){
    Read * read = get_read(r);
    int sample = read->get_sample();
    vector<int> als;
    vector< Transcript *>::iterator t1;
    cout << "cluster="<<get_id()<<"\t"<<"sample="<<sample ;
    for(t1=read->align_begin(); t1!=read->align_end(); t1++){
      als.push_back((*t1)->pos());
    }
    sort(als.begin(),als.end());
    cout << "\tthere are "<< als.size() << " reads:";
    for(vector<int>::iterator it=als.begin(); it!=als.end(); it++)
      std::cout << "\t" << *it ;
    std::cout << endl ;
  }
}
