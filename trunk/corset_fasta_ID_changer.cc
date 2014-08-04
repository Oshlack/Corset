// Copyright 2014 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** corset_fasta_ID_changer is a simple program to change the contig IDs
 ** within a fasta file. The cluster to which the contig belong will prefix
 ** the contig ID. This allows for each identification of the sequences within
 ** a cluster. e.g. using grep. Thanks to Marco Salvemini for suggesting this
 ** program
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 4 August 2014
 **/ 

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <map>
#include <stdlib.h>

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cout << endl;
  cout << "corset_fasta_ID_changer modifies a fasta file by prefixing each contig ID with its associated cluster ID" << endl;
  cout << endl;
  cout << "Usage: corset_fasta_ID_changer <cluster file> <fasta file>  >  <out file>" << endl;
  cout << "         <cluster file> - Needs to be a table with two columns: contig ID and cluster ID." << endl;
  cout << "                          If you ran corset, this should be cluster.txt" << endl; 
  cout << "         <fasta file>   - This is the file for which you want to alter the contig IDs" << endl; 
  cout << "         <out file>     - The fasta file where the output should be written." << endl;
  cout << endl;
  cout << "Please see https://code.google.com/p/corset-project/ for more information"<< endl;
  cout << endl;
}

// the real stuff starts here.
int main(int argc, char **argv){

  map<string,string> cluster_contig_map;
  if(argc!=3){
    print_usage();
    exit(1);
  }

  ifstream file;
  //Open the cluster file and get the mapping between contigs and clusters
  file.open(argv[1]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  string line;
  string contig;
  string cluster;
  while ( getline (file,line) ){
    istringstream line_stream(line);
    line_stream >> contig;
    line_stream >> cluster;
    cluster_contig_map[contig]=cluster;
  }
  file.close();

  //Now open the fasta file and replace the IDs
  file.open(argv[2]);
  if(!(file.good())){
    cout << "Unable to open file " << argv[1] << endl;
    exit(1);
  }
  while ( getline (file,line) ){
    int start=line.find(">")+1;
    if(start==1){ //if this is the ID line...
      int end=line.find_first_of("\t\n ")-1;
      string id=line.substr(start,end);
      string cluster_id=cluster_contig_map[id]; //look up the cluster ID
      if(cluster_id=="") //case when ID is missing
	cluster_id="UnknownCluster";
      string new_id=cluster_id+"--"+id; //prefix ID
      line.replace(start,end,new_id);
    }
    cout << line << endl; //output
  }
  file.close();

  return(0);
}




