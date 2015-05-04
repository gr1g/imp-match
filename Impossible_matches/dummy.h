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
	
public:
	
	DUMMY();
	
	ofstream output;
	
	// Main Steps
	
	void loadData(string dataFile); // loads data from dataFile
	void analyzeData(); // loads data from dataFile
	void postAnalysis();
	void printResults(string outputFile);
	
	// Main functions
	
	int rankedOnce();
	int rankedMutiple();
	
	int newrankedOnce(vector<vector<int> > &rankingsOfPositions);
	int newrankedMutiple(vector<vector<int> > &rankingsOfPositions);
	
	void simpleBlocks();    // calls in a loop rankedOnce() and RankedMultiple(): to eliminate easy cases of impossible matches (see web appendix)
	void petitsblocks(vector<vector<int> > &rankingsOfPositions);
	//void buildPossible_K(vector<vector<int> > &profile);
	
	void buildPossible_K(vector<int> Gamma_cand, vector<vector<int> > rankings, vector<int> J_0, vector<int> &possible_K);
	vector<vector<int> > admissibleGraph(vector<vector<int> > &rankings,vector<int> &candidates, vector<int> &matchCandidates); // to get a comprehensive matching
	void checkComprehensiveness(vector<vector<int> > &rankings, vector<int> &candidates, string step);
	void comprehensiveMatching(vector<vector<int> > &rankings, vector<int> &candidates, vector<int> &matching_candidate, vector<int> &matching_department); // construct a comprehensive matching
	vector<vector<int> > simplify_i0(vector<vector<int> > &profile, int i0, int s0); // construct rankings from profile
	void checkImpossible(vector<vector<int> > profile, int i0, int s0); // main function to check a candidate is impossible
	void buildGamma(vector<vector<int> > rankings, vector<int> candidates, int i0, vector<int> matching_candidate, vector<int> matching_department, vector<int> &J_0, vector<int> &Gamma_dep, vector<int> &Gamma_cand); // construct the set Gamma_cand
	vector<vector<int> > graphOfJ(vector<vector<int> > rankings, vector<int> candidates, vector<int> Gamma_cand, vector<int> Gamma_dep, vector<int> K);    // construct the edge set used to compute a maximum matching (rankings restricted to Gamma_cand, truncated at K)
	
	void matchingStars(vector<vector<int> > profile, int dim);   // assign candidates ranked only 1st (and at least twice) to a position.
	void getStars(vector<int> &theStars, vector<vector<int> > &theirChoices);
	void implementStarsChoices(vector<int> &chosenPositions, vector<vector<int> > &choiceSet);
	
	void nextChoices(vector<int> &theSet, vector<int> &capacity);
	
	int marketCleared(vector<vector<int> > &profile);
	
	void copyRankings(vector<vector<int> > &source, vector<vector<int> > &copy);
	
	int loopK(vector<int> K, vector<vector<int> > rankings, vector<int> candidates, vector<int> Gamma_cand, vector<int> Gamma_dep, vector<int> J_0, int s0, int i0);
	
	void notImpossibles(vector<int> candidate, vector<vector<int> > rankings, vector<int> matching_candidate, vector<vector<int> > &impossibles);
	
	
	// Maximum matching functions
	/*
	int dfs(vector<vector<int> > &rankings, int a, vector<int> &candidates, vector<int> & matching_candidate,vector<int> & matching_department, vector<vector<int> > &edge);
	int dfsExp(vector<vector<int> > &rankings, int a, vector<int> &candidates, vector<int> & matching_candidate,vector<int> & matching_department, vector<vector<int> > &edge);
	int bipMatch(vector<vector<int> > &rankings, vector<int> &candidates, vector<vector<int> > &edge, vector<int> & matching_candidate, vector<int> &matching_department);
	*/
	string theFunction;
	
	// Variables, vectors, ...
	
	int dfs(vector<vector<int> > rankings, int a, vector<int> candidates, vector<vector<int> > edge, vector<int> &matching_candidate, vector<int> &matching_department, vector<int> &visited);
	int dfsExp(vector<vector<int> > rankings, int a, vector<int> candidates, vector<vector<int> > edge, vector<int> &matching_candidate, vector<int> &matching_department);
	int bipMatch(vector<vector<int> > rankings, vector<int> candidates, vector<vector<int> > edge, vector<int> &matching_candidate, vector<int> &matching_department);

	
	
	vector<vector<int> > originalRankings; // The data for each year
//	vector<vector<int> > rankings;  // The data tailored for a pair candidate-position
	vector<vector<int> > impossibles; // Recording who is an impossible match. Same dimensions are originalRankings
//	vector<vector<int> > edge;  // Graph to compute a matching
//	vector<int> Candidates; // set of candidates that are in rankings
//	vector<int> Positions;
//	vector<int> matching_candidate; // matching_candidate[i]=j <=> the i-th candidate in Candidates[] is matched to the j-th position in rankings[][]
//	vector<int> matching_department;  // matching_department[j]=i <=> the j-th position in rankings[][] is matched to the i-th candidate in Candidates[]

//	vector<int> K;  // set of students used for the truncation (see definition of a block in the paper)
//	vector<int> possible_K; // sets of students who can be in K
//	vector<int> selectedIndex; // used to select a set K in possible_K
//	vector<int> visited; // used to compute the maximum matching
	//vector<int> J_0;    // Set J_0 (see paper)
	//vector<int> Gamma_cand; // The set that is (or not) a block ($\mathbf{J}$ in the proof, see paper)
	//vector<int> Gamma_dep; // Acceptable positions for candidates in Gamma_cand
	
	vector<vector<int> > hired;
	
	
	// Additional functions
	
	int index(const vector<int> &S, int a); // Get the index of "a" in a vector "S" (spans all coordinate starting from 0)
	int indexInRankings(const vector<int> &S, int a);   // Get the index of "candidate" a in a ranking "S" (spans all coordinate starting from 1)
	void Truncate(vector<int> &S, int a);
	bool presentInVector(const vector<int> &S, int z);  // Return 1 if "z" is in "S", 0 otherwise
	signed long int coeff(int n, int k); // Compute # of combinations for a set of size k in a set of size n.
	void next_subset(vector<int> &selectedIndex, int n, int k); // construct the next subset K.
	
	
};



#endif /* defined(__Impossible_matches__dummy__) */
