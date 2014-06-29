// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** 
 ** A Cluster object is responsible for holding the information 
 ** associated with a "superCluster" group of transcripts. It also
 ** contains the methods to cluster a "superCluster" into the final
 ** smaller cluster groups. ie. it contains the code for hierarcical
 ** clustering and expression testing. Cluster objects are first created
 ** in the MakeClusters class, then the Cluster::cluster method must be
 ** called to actually perform the clustering. See the description
 ** below for an outline of the Cluster::cluster algorithm.
 **/

#ifndef ONECLUSTER_H
#define ONECLUSTER_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <Read.h>

using namespace std;

typedef vector < vector < int > > group;
typedef vector < vector < vector < int > > > read_group;

class Cluster {

  private:
     //object which store the data
     vector< Transcript * > cluster_;
     vector< Read * > read_;
     read_group read_groups; //3D vector. Dimensions are: sample, transcript, read
     group groups;
     float ** dist;
     int id_;
     vector<int> sample_groups;

     //private methods called from "cluster"
     
     /** get_dis is responsible for returning the distance between
      ** two clusters (the clusters have positions in the distance 
      ** matrix of i and j). This first calculates the distance based on the
      ** proportion of shared reads between the clusters. It then performs
      ** a log likelihood ratio test looking for a difference in the ratio
      ** of expressions cluster i / cluster j between the different experimental
      ** groups. If this is found, then the distance is increased to the maximum. 
      ** The return values will be a distance between 0 and 1. NOTE that the
      ** value returned is actually 1-distance described in our paper. So
      ** 1 = clusters share all reads, 0 = clusters share no reads.
      **/
     float get_dist(int i, int j);

     /** find_next_pair will search through the distance matrix and find
      ** the next smalled distance value. i.e. the next two pairs of 
      ** clusters to be merged together (returned as max_i and max_j)
      **/
     float find_next_pair(int n, int & max_i, int & max_j);

     /** the merge method takes two clusters (at position i and j in
      ** the distance matrx) and will merge them together. This involves
      ** merging their list of reads together, recalculating distance
      ** based on the new allocation of reads and updating the distance matrix.
      **/
     void merge(int i, int j, int n);

     /** Before starting to cluster all transcripts are allocated to
      ** either own cluster, and the distance matrix is set-up */
     void initialise_matrix();

     void clean_up();

     /** get_counts will work out how many reads have been allocated to
      ** each cluster given the current configuration of clusters. */
     vector<int> get_counts(int s);

     /** output_cluster will write to the output "clusters.txt" and
      ** "counts.txt" files with the current configuration of clusters. */
     void output_clusters(string threshold);

 public: //basic getter/setter functions
     void add_tran(Transcript * trans); 
     void add_read(Read * read);
     Transcript * get_tran(int i);
     Read * get_read(int i);
     int n_trans();
     int n_reads();
     int get_id(){return id_;};
     void set_id(int id){id_=id;};
     void set_sample_groups(vector<int> & sgroups){ sample_groups=sgroups; };

     /** This method does all the work! A basic outline is:
      ** - initialise_matrix() is called to set up the distance matrix
      ** - Until we have merged all the transcripts into a single cluster:
      **     - Find the closest clusters to merge with find_next_pair()
      **     - Are we over the distance threshold now? If yes then
      **       call get_counts() followed by output_clusters().
      **     - merge the two closest cluster together and update the
      **       distance matrix with merge(). merge() will call get_dist() 
      **       to recalculate the distances.
      ** - When finished, free the memory we've been using with clean_up()
      **/
     void cluster(map < double, string > & thresholds);

     static double D_cut;
     static string file_prefix;
     const static string file_counts;
     const static string file_clusters;
     const static string file_ext;
     const static string cluster_id_prefix;
     const static string cluster_id_joiner;
     const static string cluster_id_prefix_no_reads;

     void print_alignments();

};

#endif



