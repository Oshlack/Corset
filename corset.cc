// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** This file contains the main function for the corset program. It
 ** controls the i/o including command line options.
 **
 ** The program runs in four main stages:
 ** 1 - All the alignments in the bam files are read. See the function
 **     read_bam_file below and add_alignment() in ReadList for more
 **     detail. 
 ** 2 - Filter out any transcripts which have fewer than a certain number
 **     of reads aligning. By deafult this is 10.
 ** 3 - The read-transcript associations are parsed and the transcripts
 **     undergo an initial clustering, see MakeClusters::makeSuperClusters(),
 **     in which transcripts are grouped together if they share a single read
 **     with another transcript in the same group.
 ** 4 - For each of the clusters above, hierarchical clustering is performed,
 **     see Cluster::cluster(). The group is split into multiple smaller
 **     clusters and counts are output for each of these smaller groups.
 **
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 20 January 2016
 **/ 

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <map>
//#include <pthread.h>

#include <MakeClusters.h>
#include <Read.h>
#include <Transcript.h>

#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <signal.h>

#include <sam.h>
#include <bam.h>

//#include <gperftools/profiler.h>

#define MAX_BAM_LINE_LENGTH 10000

using namespace std;

const string corset_extension=".corset-reads";

// a function to parse a bam file. It required samtools
// to read it. The alignments are stored in a ReadList object.
ReadList * read_bam_file(string all_file_names, TranscriptList * trans, int sample){

   ReadList * rList = new ReadList(trans);
   string filename;
   stringstream ss(all_file_names);

   //read each file in the comma separated list
   while(getline(ss,filename,',')){  
     cout << "Reading bam file : "<< filename << endl;
     if(!(filename.find(corset_extension)==string::npos) ){
       cout << "The input files look like corset read files (filenames with "<< corset_extension;
       cout << "). Perhaps you meant to run with the -i corset option?" << endl;
       exit(1);
     }
     samfile_t *in = 0 ;
     if ((in = samopen(filename.c_str(), "br", NULL)) == 0 | in->header == 0) {
       cerr << "fail to open "<< filename << " for reading." << endl;
       exit(1);
     }
     bam1_t *b = bam_init1();
     int r;
     int i=0;

     //read the header first..... only really need this when we want to
     //output all the transcripts including those with no reads
     for(int tid=0; tid < in->header->n_targets; tid++)
       trans->insert(in->header->target_name[tid]);

     while ((r = samread(in, b)) >= 0) { // read one alignment from `in'
       //add the alignment into our data object
       string read_name=string(bam1_qname(b));
       int tid=b->core.tid;
       if( tid != -1 )  //unmapped reads have tid=-1
	 rList->add_alignment(read_name, string(in->header->target_name[tid]),sample); 
       if(i % 200000 == 0) //output some info about how it's proceeding.
	 cout << float(i)/float(1000000) << " million alignments read" <<endl; 
       //       if(i>5000000) break;
       i++;
     }
     bam_destroy1(b);
     samclose(in);
     cout << "Done reading "<< filename << endl;
   }
   rList->compactify_reads(trans); //reduce the memory usage
   return rList; //return the alignment list
}

ReadList * read_corset_file(string all_file_names, TranscriptList * trans, int sample){
   ReadList * rList = new ReadList(trans);
   string filename;
   stringstream ss(all_file_names);

   //read each file in the comma separated list
   while(getline(ss,filename,',')){  
     cout << "Reading corset file : "<< filename << endl;
     //check that the filename has the correct file extension
     //i.e. that we can't accidently pass a bam file.
     if(filename.find(corset_extension)==string::npos ){
       cout << "The input files don't have the extension, "<< corset_extension;
       cout << ". Please check them." << endl;
       exit(1);
     }
     ifstream file(filename);
     string line;
     while(getline(file, line)){
	 istringstream istream(line);
	 int weight;
	 vector<string> transNames;
	 string name;
	 istream >> weight;
	 while(istream >> name)
	   transNames.push_back(name);
	 rList->add_alignment(transNames,sample,weight);
     }
   }
   return rList;
}

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cout << endl;
  cout << "Corset Version "<< VERSION << endl;
  cout << endl;
  cout << "Corset clusters contigs and counts reads from de novo assembled transcriptomes." << endl;
  cout << endl;
  cout << "Usage: corset [options] <input bam files>" << endl;
  cout << endl;
  cout << "Input bam files:" << endl;
  cout << "\t The input files should be multi-mapped bam files. They can be single, paired-end or mixed" << endl;
  cout << "\t and do not need to be indexed. A space separated list should be given." << endl;
  cout << "\t e.g. corset sample1.bam sample2.bam sample3.bam" << endl;
  cout << "\t or just: corset sample*.bam" << endl;
  cout << endl;
  cout << "\t If you want to combine the results from different transcriptomes. i.e. the same reads have " << endl;
  cout << "\t been mapped twice or more, you can used a comma separated list like below:" << endl;
  cout << "\t corset sample1_Trinity.bam,sample1_Oases.bam sample2_Trinity.bam,sample2_Oases.bam ..." << endl;
  cout << endl;
  cout << "Options are:" << endl;
  cout << endl;
  cout << "\t -d <double list> A comma separated list of distance thresholds. The range must be" << endl;
  cout << "\t                  between 0 and 1. e.g -d 0.4,0.5. If more than one distance threshold" << endl;
  cout << "\t                  is supplied, the output filenames will be of the form:" << endl;
  cout << "\t                  counts-<threshold>.txt and clusters-<threshold>.txt " << endl;
  cout << "\t                  Default: 0.3"  << endl;
  cout << endl;
  cout << "\t -D <double>      The value used for thresholding the log likelihood ratio. The default " << endl;
  cout << "\t                  value will depend on the number of degrees of freedom (which is the " << endl;
  cout << "\t                  number of groups -1). By default D = 17.5 + 2.5 * ndf, which corresponds " << endl;
  cout << "\t                  approximately to a p-value threshold of 10^-5, when there are fewer than" << endl;
  cout << "\t                  10 groups." << endl;
  cout << endl;
  cout << "\t -m <int>         Filter out any transcripts with fewer than this many reads aligning." << endl;
  cout << "\t                  Default: 10" << endl;
  cout << endl;
  cout << "\t -g <list>        Specifies the grouping. i.e. which samples belong to which experimental" << endl;
  cout << "\t                  groups. The parameter must be a comma separated list (no spaces), with the " << endl;
  cout << "\t                  groupings given in the same order as the bam filename. For example:" << endl;
  cout << "\t                  -g Group1,Group1,Group2,Group2 etc. If this option is not used, each sample" << endl;
  cout << "\t                  is treated as an independent experimental group." << endl;
  cout << endl;
  cout << "\t -p <string>      Prefix for the output filenames. The output files will be of the form" <<endl;
  cout << "\t                  <prefix>-counts.txt and <prefix>-clusters.txt. Default filenames are:" << endl;
  cout << "\t                  counts.txt and clusters.txt" << endl;
  cout << endl;
  cout << "\t -f <true/false>  Specifies whether the output files should be overwritten if they already exist." << endl;
  cout << "\t                  Default: false" << endl;
  cout << endl;
  cout << "\t -n <string list> Specifies the sample names to be used in the header of the output count file." << endl;
  cout << "\t                  This should be a comma separated list without spaces." << endl;
  cout << "\t                  e.g. -n Group1-ReplicateA,Group1-ReplicateB,Group2-ReplicateA etc." << endl;
  cout << "\t                  Default: the input filenames will be used." << endl;
  cout << endl;
  cout << "\t -r <true/true-stop/false> " << endl;  
  cout << "\t                  Output a file summarising the read alignments. This may be used if you" << endl;
  cout << "\t                  would like to read the bam files and run the clustering in seperate runs" << endl;
  cout << "\t                  of corset. e.g. to read input bam files in parallel. The output will be the" << endl;
  cout << "\t                  bam filename appended with "<< corset_extension <<"." << endl;
  cout << "\t                  Default: false" << endl;
  cout << endl;
  cout << "\t -i <bam/corset>  The input file type. Use -i corset, if you previously ran" << endl ;
  cout << "\t                  corset with the -r option and would like to restart using those" << endl;
  cout << "\t                  read summary files. Running with -i corset will switch off the -r option." << endl;
  cout << "\t                  Default: bam" << endl;
  cout << endl;
  cout << "Citation: Nadia M. Davidson and Alicia Oshlack, Corset: enabling differential gene expression " << endl;
  cout << "          analysis for de novo assembled transcriptomes, Genome Biology 2014, 15:410" << endl;
  cout << endl;
  cout << "Please see https://github.com/Oshlack/Corset/wiki for more information"<< endl;
  cout << endl;
}

// the real stuff starts here.
int main(int argc, char **argv){

  //default parameters:
  string distance_string("0.3");
  bool force=false;
  string sample_names;
  int c;
  int params=1;
  vector<int> groups;
  bool output_reads=false;
  bool stop_after_read=false;
  string input_type="bam";
  //function pointer to the method to read the bam or corset input files
  ReadList * (*read_input)(string, TranscriptList *, int) = read_bam_file;

  cout << endl;
  cout << "Running Corset Version "<< VERSION << endl;

  //parse the command line options
  while((c =  getopt(argc, argv, "f:p:d:n:g:D:m:r:i:")) != EOF){
    switch(c){
    case 'f': { //f=force output to be overwritten?
      std::string value(optarg); 
      transform(value.begin(), value.end(), value.begin(), ::tolower);
      if( value.compare("true")==0 | value.compare("t")==0 | value.compare("1")==0 ){
	force=true;
        cout << "Setting output files to be overridden"<<endl;
      }
      else if (value.compare("false")==0 | value.compare("f")==0 | value.compare("0")==0 )
	force=false;
      else {
	cerr << "Unknown argument passed with -f. Please specify true or false." << endl;
	print_usage();
	exit(1);
      }
      params+=2;
      break; }
    case 'p': //p=the prefix for the output files
      cout << "Setting output filename prefix to " << optarg << endl;
      Cluster::file_prefix=string(optarg)+"-";
      params+=2;
      break;
    case 'd': //d=the distance threshold for the clustering
      cout << "Setting distance threshold to: " << optarg << endl;
      distance_string = optarg;
      params+=2;
      break; 
    case 'n':{ //names to give the samples in the output files
      cout << "Setting sample names to:" << optarg << endl;
      sample_names=string(optarg);
      replace( sample_names.begin(), sample_names.end(), ',', '\t');
      params+=2;
      break;
    }
    case 'g':{ //the sample groups. Used for the likelihood ratio test
      cout << "Setting sample groups:" << optarg;
      stringstream ss(optarg);
      string s;
      vector < string > names;
      while (getline(ss, s, ','))
	names.push_back(s);
      int ngroups=0;
      for(int s1=0; s1 < names.size(); s1++){
	//find the first occurance
	int s2=0;
	for(; s2 < names.size(); s2++){
	  if(names.at(s1).compare(names.at(s2))==0)
	    break;
	}
	if(groups.size()<=s2){
	  groups.push_back(ngroups);
	  ngroups++;
	}
	else
	  groups.push_back(groups.at(s2));
      }
      Transcript::groups=ngroups;
      cout << ", " << Transcript::groups << " groups in total" << endl;
      params+=2;
      break;
    }
    case 'D':{
      cout << "Setting likelihood threshold at "<<optarg<<endl;
      Cluster::D_cut=atof(optarg);
      params+=2;
      break;
    }
    case 'm':{
      cout << "Setting minimum counts to "<<optarg<<endl;
      Transcript::min_counts=atoi(optarg);
      params+=2;
      break;
    }
    case 'r':{ //r=whether to output read information
      std::string value(optarg);
      transform(value.begin(), value.end(), value.begin(), ::tolower);
      if( value.compare("true-stop")==0){
	output_reads=true;
	stop_after_read=true;
	cout << "Setting read alignments to be output to file." << endl;
	cout << "Corset will exit after reading bam files." << endl;
      }
      else if( value.compare("true")==0 | value.compare("t")==0 ){
        output_reads=true;
        cout << "Setting read alignments to be output to file."<<endl;
      }
      else if (value.compare("false")!=0 & value.compare("f")!=0 ){
        cerr << "Unknown argument passed with -r. Please specify true or false." << endl;
        print_usage();
        exit(1);
      }
      params+=2;
      break; 
    }
    case 'i':{
      std::string value(optarg);
      transform(value.begin(), value.end(), value.begin(), ::tolower);
      if( value.compare("corset")==0 ){ 
	read_input=read_corset_file;
	output_reads=false; //override the -r option
	stop_after_read=false;
      }
      else if(value.compare("bam")!=0 ){
	cerr << "Unknown input type, "<< value<<", passed with -i. Please check options." << endl;
	print_usage();
	exit(1);
      }
      params+=2;
      break;
    }
    case '?': //other unknown options
      cerr << "Unknown option.. stopping" << endl;
      print_usage();
      exit(1);
      break;
    }
  }
  
  //convert the distances provided as a string into doubles
  map<float,string> distance_thresholds;
  char * cstr = new char[distance_string.length()+1];
  strcpy(cstr, distance_string.c_str());
  char * single_string = strtok(cstr,",");
  while (single_string != NULL){ 
    double value = atof(single_string);
    //check that the value is valid
    if(value<0|value>1){
      cerr << "The distance provided is invalid. It must be between 0 and 1 (inclusive)." << endl;
      print_usage();
      exit(1);
    }
    distance_thresholds[ value ]=string("-")+single_string;
    single_string = strtok(NULL,",");
  } 
  if(distance_thresholds.size()==1)
    distance_thresholds.begin()->second=""; //no suffix if only one distance cut-off

  //the remaining command line arguments are names of the bam files
  int smpls=argc-params;
  if(smpls==0){
    cerr << "No bam files specified" << endl;
    print_usage();
    exit(1);
  }
  if(sample_names==""){ //if the -n option (sample names) have not been 
    //specified, then use the filenames.
    for(int s = params; s < argc ; s++){
      sample_names+=string(argv[s]);
      if(s < (argc-1)) sample_names+='\t';
    }
  }
  //if the number of sample provided for the header of the counts file
  //is inconsistent with the number of samples print an error and exit
  int n_sample_names=count(sample_names.begin(), sample_names.end(),'\t')+1;
  if(n_sample_names!=smpls){
    cerr << "The number of sample names passed (via the -n option), " << n_sample_names 
	 << ", does not match the number of samples, " << smpls
	 << ". Please check how many values you have passed." << endl;
    exit(1);
  }

  if(groups.size()==0){ // if no -g option (no groups information) is given 
    //treat each file as a separate group
    for(int s=0; s<smpls; s++)
      groups.push_back(s);
    Transcript::groups=smpls;
  }

  //before we go any further. Check that there are enough groups and
  //sample names for all the files
  if(groups.size()!=smpls){
    cerr << "The number of experimental groups passed (via the -g option) "
	 << "does not match the number of input bam files. Please check how "
	 << "many values you have passed." << endl;
    exit(1);
  }

  //if D was not set, then set the default value
  if(Cluster::D_cut==0)
    Cluster::D_cut=17.5+2.5*(groups.size()-1);      
  
  //test whether the output files already exist
  for (map<float,string>::iterator it=distance_thresholds.begin(); 
       stop_after_read==false & it!=distance_thresholds.end(); ++it){
    string type[]={Cluster::file_counts,Cluster::file_clusters};
    for(int t=0; t < 2; t++){
      string filename(Cluster::file_prefix+type[t]+it->second+Cluster::file_ext);
      ifstream file(filename.c_str());
      if(file.good()){
	if(force & remove( filename.c_str() )!=0){
	  cerr << "Could not replace the file, "<< filename << endl;
	  exit(1);
	}
	else if(!force){
	  cerr << "File already exists, " << filename << ". Use \"-f true\" option to overwrite " << endl;
	  exit(1);
	}
      }
      file.close();
      ofstream ofile(filename.c_str());
      //output the header for the counts files.
      if(type[t].compare(Cluster::file_counts)==0){
	ofile << '\t' + sample_names << endl;
      }
      ofile.close();
    }
  }

  //set the number of samples
  Transcript::samples=smpls;
  //set-up out data structures ready to receive the alignment information
  TranscriptList * tList = new TranscriptList;
  vector<ReadList*> rList; //one ReadList per sample  

  for(int bam_file=0; bam_file < smpls; bam_file++){
    rList.push_back(read_input(string(argv[params+bam_file]),tList,bam_file));
    //This is where we output the read alignment summary file for future runs of corset
    if(output_reads){
      string bam_filename = string(argv[params+bam_file]) ;
      string base_filename = bam_filename.substr(bam_filename.find_last_of("/\\") + 1) ;
      rList[bam_file]->write(base_filename+corset_extension);
    }
  } 

  cout << "Done reading all files. "<< endl;
  
  if(stop_after_read==true) exit(0);

  //output clusters with no reads in this special case
  StringSet<Transcript>::iterator it;
  int n=0;
  for(it=tList->begin(); it!=tList->end(); it++){
    if(!(*it).second->reached_min_counts()){
      if(Transcript::min_counts==0){ 
	n++;
	//output clusters
	for (map<float,string>::iterator dis=distance_thresholds.begin(); dis!=distance_thresholds.end(); ++dis){
	  ofstream clusterFile;
	  string filename(Cluster::file_prefix+Cluster::file_clusters+dis->second+Cluster::file_ext);
	  clusterFile.open(filename.c_str(),ios_base::app);
	  clusterFile << (*it).second->get_name() << "\t" << Cluster::cluster_id_prefix_no_reads << n << endl;
	  clusterFile.close();
	}
      }
      (*it).second->remove();
    }
  };

  delete tList; //don't need this anymore, so why not free some space  

  cout << "Start to cluster the reads" << endl;

  //Now the reads are parsed and the clustering and counting is performed.
  //  ProfilerStart("gprof.prof");
  MakeClusters cList(rList, distance_thresholds,groups); 
  //  ProfilerStop();

  cout << "Finished" << endl;

  return(0);
}




