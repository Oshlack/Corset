// Copyright 2013 Nadia Davidson for Murdoch Childrens Research
// Institute Australia. This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

#ifndef MAKECLUSTERS_H
#define MAKECLUSTERS_H

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <map>
#include <Cluster.h>
#include <Read.h>
#include <Transcript.h>

using namespace std;

class MakeClusters {
   private:
     vector<Cluster*> clusterList;
     Cluster * current_cluster;

   public:
     void setCurrentCluster(Transcript * trans){ current_cluster=getMapElement(trans)->second; };

 private:
     map< Transcript *, Cluster * > transMap;
     pair< Transcript * const, Cluster * > * getMapElement(Transcript * trans);
     void checkAgainstCurrentCluster(Transcript * trans);
     void makeSuperClusters(vector<ReadList*> & readLists);
     void processSuperClusters(map<double,string> & distance_thresholds, vector<int> & groups);
     
 public:
     MakeClusters(vector<ReadList*> & readLists, map<double,string> & distance_thresholds, vector<int> & groups);
     
};

#endif


