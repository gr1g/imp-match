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
    
public:
    
    DUMMY();
    
    ofstream output;
    ofstream theMatchings;
    ofstream details;
    
    // Main Steps
    void loadData(string dataFile); // loads data from dataFile
    void stepOne(); // main algorithm
    void stepTwo(); // refinement with stars
    void printResults(string outputFile);
    
    // Main functions
    
    int rankedOnce();
    int rankedMutiple();
    
    int newrankedOnce(vector<vector<int> > &rankingsOfPositions);
    int newrankedMutiple(vector<vector<int> > &rankingsOfPositions);
    
    // petitsblocks = calls in a loop rankedOnce() and RankedMultiple(): to eliminate easy cases of impossible matches (see web appendix)
    void petitsblocks(vector<vector<int> > &rankingsOfPositions);

    // construct rankings from profile
    vector<vector<int> > simplify_i0(vector<vector<int> > profile, int i0, int s0);
    
    // construct a comprehensive matching
    void comprehensiveMatching(vector<vector<int> > rankings, vector<int> candidates, vector<int> &matching_candidate, vector<int> &matching_department);
    void checkComprehensiveness(vector<vector<int> > rankings, vector<int> candidates, string step);
    // admissibleGraph = main function to check a candidate is impossible
    vector<vector<int> > admissibleGraph(vector<vector<int> > rankings,vector<int> candidates, vector<int> matchCandidates);
    
    // main function to check a candidate is impossible
    void checkImpossible(vector<vector<int> > profile, int i0, int s0, vector<vector<int> > &impossible);
    
    vector<int> buildPossible_K(vector<int> Gamma_cand, vector<vector<int> > rankings, vector<int> J_0);
    
    // construct the set Gamma_cand
    void buildGamma(vector<vector<int> > rankings, vector<int> candidates, int i0, vector<int> matching_candidate, vector<int> matching_department, vector<int> &J_0, vector<int> &Gamma_dep, vector<int> &Gamma_cand);
    
    // construct the edge set used to compute a maximum matching (rankings restricted to Gamma_cand, truncated at K)
    vector<vector<int> > graphOfJ(vector<vector<int> > rankings, vector<int> candidates, vector<int> Gamma_cand, vector<int> Gamma_dep, vector<int> K);
    
    // construct the edge set used to compute a maximum matching (rankings restricted to Gamma_cand, truncated at K)
    void matchingStars(vector<vector<int> > profile, int dim, vector<int> &theLoops);
    void getStars(vector<vector<int> > profile, vector<int> &theStars, vector<vector<int> > &theirChoices);
    vector<vector<int> > implementStarsChoices(vector<int> chosenPositions, vector<vector<int> > choiceSet, vector<vector<int> > profile);
    
    void nextChoices(vector<int> &theSet, vector<int> &capacity);
    
    int marketCleared(vector<vector<int> > &profile);
    void copyRankings(vector<vector<int> > &source, vector<vector<int> > &copy);
    int loopK(vector<int> K, vector<vector<int> > rankings, vector<int> candidates, vector<int> Gamma_cand, vector<int> Gamma_dep, vector<int> J_0, int s0, int i0);
    void notImpossibles(vector<int> candidate, vector<vector<int> > rankings, vector<int> matching_candidate, vector<vector<int> > &impossibles);
    
    bool thereAreStars(vector<vector<int> > profile);
    
    // Maximum matching functions
    int dfs(vector<vector<int> > rankings, int a, vector<int> candidates, vector<vector<int> > edge, vector<int> &matching_candidate, vector<int> &matching_department, vector<int> &visited);
    int dfsExp(vector<vector<int> > rankings, int a, vector<int> candidates, vector<vector<int> > edge, vector<int> &matching_candidate, vector<int> &matching_department);
    int bipMatch(vector<vector<int> > rankings, vector<int> candidates, vector<vector<int> > edge, vector<int> &matching_candidate, vector<int> &matching_department);
    
    
    // Additional functions
    
    int index(const vector<int> S, int a); // Get the index of "a" in a vector "S" (spans all coordinate starting from 0)
    int indexInRankings(const vector<int> S, int a);   // Get the index of "candidate" a in a ranking "S" (spans all coordinate starting from 1)
    void Truncate(vector<int> &S, int a);
    bool presentInVector(const vector<int> S, int z);  // Return 1 if "z" is in "S", 0 otherwise
    signed long int coeff(int n, int k); // Compute # of combinations for a set of size k in a set of size n.
    void nextSubset(vector<int> &selectedIndex, int n, int k); // construct the next subset K.
    
    // Variables, vectors, ...
    string theFunction; // for debugging purposes
    int maxIter;
    vector<vector<int> > originalRankings; // The data for each year
    vector<vector<int> > originalImpossible; // Recording who is an impossible match. Same dimensions are originalRankings
    vector<vector<int> > hired;
    vector<int> stepOnePredicted;
    
    
    int numberPositionsCleared;
    int totalNumberCombinations;
    int maxDim;
    int firstStep;
    int preferencesSimulated;
    int year;
    
    string fileNameResults;
    string fileNameMatchings;
    int candidatesDone;
    
};

#endif /* defined(__Impossible_matches__dummy__) */

// OLD
//	vector<vector<int> > rankings;  // The data tailored for a pair candidate-position
//	vector<vector<int> > edge;  // Graph to compute a matching
//	vector<int> Candidates; // set of candidates that are in rankings
//	vector<int> Positions;
//	vector<int> matching_candidate; // matching_candidate[i]=j <=> the i-th candidate in Candidates[] is matched to the j-th position in rankings[][]
//	vector<int> matching_department;  // matching_department[j]=i <=> the j-th position in rankings[][] is matched to the i-th candidate in Candidates[]

//	vector<int> K;  // set of students used for the truncation (see definition of a block in the paper)
//	vector<int> possible_K; // sets of students who can be in K
//	vector<int> selectedIndex; // used to select a set K in possible_K
//	vector<int> visited; // used to compute the maximum matching
//	vector<int> J_0;    // Set J_0 (see paper)
//	vector<int> Gamma_cand; // The set that is (or not) a block ($\mathbf{J}$ in the proof, see paper)
//	vector<int> Gamma_dep; // Acceptable positions for candidates in Gamma_cand