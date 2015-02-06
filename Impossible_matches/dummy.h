//
//  dummy.h
//  Impossible_matches
//
//  Created by Guillaume Haeringer on 4/28/14.
//  Copyright (c) 2014 Guillaume Haeringer. All rights reserved.
//

#ifndef __Impossible_matches__dummy__
#define __Impossible_matches__dummy__


#include <iostream>
#include <fstream>

#include <string>
#include <vector>

using namespace std;


class DUMMY {
private:
    int numGone; // whenever you mark someone as gone, increment this. When numGone=number of agents, stop TTC.

public:
    
    DUMMY();

    // Main Steps
    
    void loadData(string dataFile); // loads data from dataFile
    void analyzeData(); // loads data from dataFile
    void postAnalysis();
    void printResults(string outputFile);


    // Main functions
    
    int RankedOnce();
    int RankedMutiple();
    void SimpleBlocks();    // calls in a loop RankedOnce() and RankedMultiple(): to eliminate easy cases of impossible matches (see web appendix)
    void buildPossible_K();
    void AdmissibleGraph(); // to get a comprehensive matching
    void CheckComprehensiveness();
    void ComprehensiveMatching(); // construct a comprehensive matching
    void Simplify_i0(int i0, int s0); // construct rankings from original_rankings
    void CheckImpossible(int i0, int s0); // main function to check a candidate is impossible
    void buildGamma(vector<int> &J_0, int s0, int i0); // construct the set Gamma_cand
    void buildJ_0(int i0, int s0);  // construct the set J_0
    void GraphOfJ();    // construct the edge set used to compute a maximum matching (rankings restricted to Gamma_cand, truncated at K)

    
    // Maximum matching functions
    int dfs(int a);
    int dfsExp(int a);
    int bipMatch();
   
    // Variables, vectors, ...
    
    vector<vector<int> > original_rankings; // The data for each year
    vector<vector<int> > rankings;  // The data tailored for a pair candidate-position
    vector<vector<int> > impossibles; // Recording who is an impossible match. Same dimensions are original_rankings
    vector<vector<int> > edge;  // Graph to compute a matching
    vector<int> Candidates; // set of candidates that are in rankings
    vector<int> Positions;
    vector<int> matching_candidate; // matching_candidate[i]=j <=> the i-th candidate in Candidates[] is matched to the j-th position in rankings[][]
    vector<int> matching_department;  // matching_department[j]=i <=> the j-th position in rankings[][] is matched to the i-th candidate in Candidates[]
    vector<vector<int> > OpePostes; // the data from the file
    int TempRanking[20]; // used to create original_rankings from OpePostes
    vector<int> K;  // set of students used for the truncation (see definition of a block in the paper)
    vector<int> possible_K; // sets of students who can be in K
    vector<int> SelectedIndex; // used to select a set K in possible_K
    vector<int> visited; // used to compute the maximum matching
    vector<int> J_0;    // Set J_0 (see paper)
    vector<int> Gamma_cand; // The set that is (or not) a block ($\mathbf{J}$ in the proof, see paper)
    vector<int> Gamma_dep; // Acceptable positions for candidates in Gamma_cand
    unsigned long sizeK;
    signed long int LoopNumber;


    // Additional functions
    
    int Index(const vector<int> &S, int a); // Get the index of "a" in a vector "S" (spans all coordinate starting from 0)
    int IndexInRankings(const vector<int> &S, int a);   // Get the index of "candidate" a in a ranking "S" (spans all coordinate starting from 1)
    void Truncate(vector<int> &S, int a);
    bool presentInVector(const vector<int> &S, int z);  // Return 1 if "z" is in "S", 0 otherwise
    signed long int coeff(int n, int k); // Compute # of combinations for a set of size k in a set of size n.
    void next_subset(vector<int> &SelectedIndex, int n, int k); // construct the next subset K.

    
};



#endif /* defined(__Impossible_matches__dummy__) */
