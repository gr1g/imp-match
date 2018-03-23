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
    ofstream details;
    ofstream preferences;
    ofstream randomClearing;
    ofstream finalRestults;
    
    // Main Steps
    void loadData(string dataFile); // loads data from dataFile
    void stepOne(int year); // main algorithm
    void stepTwoReal(); // refinement with stars
    vector<vector<int> > readRankings();
    void predictions(vector<vector<int> > &rankings, string step, int numberRealChoices, int numberPseudoStars);
    
    vector<int> clearedBatch;
    int howManyCleared;
    int totalNumberCandidates;
    int theStep;
    
    // Main functions
    
    void statisticsHired();
    
    
    //////////////////////
    //
    //  stepOne
    //
    
    // singletonBlocks = calls in a loop rankedOnce() and RankedMultiple(): to eliminate easy cases of impossible matches (see web appendix)
    void singletonBlocks(vector<vector<int> > &rankingsOfPositions);
    int newrankedOnce(vector<vector<int> > &rankingsOfPositions);
    
    void newCheckImpossible(const vector<vector<int> > &profile, const int &i0, const int &s0, vector<vector<int> > &impossible);
    
    void matchProcess(const vector<vector<int> > &thisrankings, const vector<int> & matchingDep, bool &isAnImpossible, vector<int> &myMatchingDep, const vector<int> & candToBeMatched);

    void newrankedMutiple(vector<vector<int> > &rankingsOfPositions);

	
	void randomProfile(vector<int> virtualCandidatesQueue, vector<vector<int> > &originalRankings);
	void reshuffleRankings(const vector<vector<int> > &positionsLinked);
    
    //
    //  stepOne
    //
    //////////////////////
    
    //////////////////////
    //
    //  stepTwoReal
    //
    void matchingStarsRealChoices(vector<vector<int> > &profile, int dim, int &numberRealChoices, int &numberPseudoStars, int &numberStars);
    void implementStarsRealChoices(vector<vector<int> > &profile, int &numberRealChoices, int &numberStars);
    void implementPseudoStarsRealChoices(vector<vector<int> > &profile, int &numberRealChoices, int &numberPseudoStars);
    bool anotherStarsLoop(vector<vector<int> > profile);
    bool findAndCleanImpossible(vector<vector<int> > &profile);
    void getChoices(const vector<vector<int> > & thisrankings, const vector<int> & matchingDep, int &toBeMatched, vector<int> &choiceSet, const vector<int> & candToBeMatched);
    void getStars(vector<vector<int> > profile, vector<int> &theStars, vector<vector<int> > &theirChoices);
    void findCandidate(const vector<vector<int> > &thisrankings, const vector<int> &matchingDep, int &toBeMatched);
    //
    //  stepTwoReal
    //
    //////////////////////
    
    // Additional functions
    vector<int> marketCleared(vector<vector<int> > &profile);
    int index(const vector<int> &S, int a); // Get the index of "a" in a vector "S" (spans all coordinate starting from 0)
    int indexInRankings(const vector<int> &S, int a);   // Get the index of "candidate" a in a ranking "S" (spans all coordinate starting from 1)
    void Truncate(vector<int> &S, int a); // Given a ranking S[], delete a and all entries below
    bool presentInVector(const vector<int> &S, int z);  // Return 1 if "z" is in "S", 0 otherwise
    
    // Variables, vectors, ...
    string theFunction; // for debugging purposes
    int maxIter;
    vector<vector<int> > originalRankings; // The data for each year
    vector<vector<int> > originalImpossible; // Recording who is an impossible match. Same dimensions are originalRankings
    vector<vector<int> > hired;
    vector<vector<int> > linkedPositions(vector<vector<int> > rankings);
    
    int numberPositionsCleared;
    int maxDim;
    int year;
    
    string fileNameResults;
    string fileRankingsStepOne;
    int candidatesDone;
    bool firstStep;
    string computeTime(time_t start, time_t end);
    

    void printRankings(vector<vector<int> > rankings);
    
};

#endif /* defined(__Impossible_matches__dummy__) */
