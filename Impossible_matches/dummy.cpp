//
//  dummy.cpp
//  Impossible_matches
//
//  Created by Guillaume Haeringer on 4/28/14.
//  Copyright (c) 2014 Guillaume Haeringer. All rights reserved.
//

#include "dummy.h"

#include <iostream>
#include <iomanip>
#include <sstream>      // std::istringstream
#include <string>       // std::string
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <random>

bool impossible=true;
bool decided = false;

DUMMY::DUMMY()
{
    maxIter = 10;
}



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Some useful functions                                            //
//                                                                      //
//                                                                      //


// Give the index number of "a" in a vector "S"
int DUMMY::index(const vector<int> &S, int a){
	auto it_found = (std::find(S.begin(), S.end(), a));
	if (it_found==S.end()) {
		return -1;
	}
	size_t pose = std::distance(S.begin(), it_found);
	return (int)(pose);
}


// Give the index number of "a" in a vector "S" except first entry
// rankings[j][0] is the ID # of position j, could be the same as that of a candidate
int DUMMY::indexInRankings(const vector<int> &S, int a){
    int  ans=-1;
    if (S.size()>1) {
        for (int i = 1 ; i < S.size(); i++) {
            if (S[i]==a) {
                ans=i;
            }
        }
    }
    return(ans);
}


// Check if "z" is present in a vector "S"
bool DUMMY::presentInVector(const vector<int> &S, int z) {
    bool found = (std::find(S.begin(), S.end(), z) != S.end());
    return found;
}

// Given the ranking of a position (vector S), truncate at candidate "a"
// (i.e., "a" is deleted from the ranking and all the candidates below as well
// "a" is the ID (or NAME) of the candidate, it is NOT its rank number in the ranking of the department.
void DUMMY::Truncate(vector<int> &S, int a){
    if (indexInRankings(S, a)>-1) {
        S.erase(S.begin()+indexInRankings(S,a),S.end());
    } else {
        cout << "Error, cannot truncate, " << a << " not present in vector of size " << S.size() << "\n";
        for (int i = 0 ; i < S.size() ; i++) {
            cout << S[i] << "\t";
        }
        cout << "\n";
        cout << "index = " <<  indexInRankings(S, a) << "\n";
        cout << " theFunction = " << theFunction << "\n";
        exit(0);
    }
}

string DUMMY::computeTime(time_t start, time_t end) {
    string duration = "";
    string second, minute, hour;
    long hours,minutes,seconds;
    second = "second";
    minute = "minute";
    hour = "hour";
    
    hours = (end-start)/3600;
    minutes = (end-start)/60-60*hours;
    seconds = end-start-60*minutes-3600*hours;
    if (minutes > 1 ) {
        minute = "minutes";
    }
    if (seconds > 1 ) {
        second = "seconds";
    }
    if (hours > 1 ) {
        hour = "hours";
    }
    duration = "\nTime spent: ";
    if (hours < 1) {
        if (minutes < 1 ) {
            duration = duration + to_string(seconds) + " " + second + "\n";
        } else {
            duration = duration + to_string(minutes) + " " + minute + ", " + to_string(seconds) + " " + second + "\n";
        }
    } else {
        duration = duration + to_string(hours) + " " + hour + ", " + to_string(minutes) + " " + minute + ", " + to_string(seconds) + " " + second + "\n";
    }
    return duration;
}


//                                                                      //
//     Some useful functions                                            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Truncations at candidates ranked only once                       //
//                                                                      //
//                                                                      //

//
//  Find a candidate only ranked once. Then truncate just below that candidate
//  Then reiterate until there's no such candidate anymore
//
int DUMMY::newrankedOnce(vector<vector<int> > &rankingsOfPositions) {
    theFunction="newRankedOnce";
    
    int DidTruncate;
    int i,j,ii,jj;
    int found,found2,count;
    
    // DidTruncate will serve for the return of the function.
    DidTruncate=0;
    found=1;
    
    size_t numberPositions = rankingsOfPositions.size();
    
    while (found>0) {
        found2=0;
        for (j = 0 ; j < numberPositions ; j++) {
            size_t sizeRanking = rankingsOfPositions[j].size();
            for (i = 1 ; i < sizeRanking ; i++) {
                count=0;
                // count = number of times the candidate is ranked
                for (jj = 0 ; jj < numberPositions ; jj++) {
                    size_t sizeOtherRanking = rankingsOfPositions[jj].size();
                    for (ii = 1 ; ii < sizeOtherRanking ; ii++) {
                        if (rankingsOfPositions[jj][ii]==rankingsOfPositions[j][i]) {
                            count=count+1;
                        }
                    }
                }
                if (count==1) {
                    //  the candidate is ranked only once
                    //  Now check is not the last ranked (otherwise there's nothing to truncate
                    if (i < sizeRanking-1) {
                        //  She/he's not the last ranked
                        //  Truncate at the candidate ranked just after him = originalRankings[j][i+1]
                        Truncate(rankingsOfPositions[j], rankingsOfPositions[j][i+1]);
                        DidTruncate=1;
                        found2=1;
                        break;
                    }
                }
            }
        }
        if (found2==0) {
            found=0;
        }
    }
    return DidTruncate;
}

//                                                                      //
//     Truncations at candidates ranked only once                       //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Truncations at candidates ranked multiple times                  //
//                                                                      //
//                                                                      //


//
//  If k candidates are ranked only k times and for the same k departments
//  Then any candidate ranked below those k candidates for some position is an
//  impossible match. These can can be deleting without affecting the set
//  of impossible matches for the whole profile (see web Appendix)
//


void DUMMY::newrankedMutiple(vector<vector<int> > &rankingsOfPositions) {
    theFunction="newrankedMultiple";
    
    int i,j,ii;
    int found,maxSize;
    int DidTruncate;
    
    vector<vector<int> > numberRanked;
    // for each candidate numberRanked[][] stores the position that ranked him/her
    // numberRanked[i].size() gives the number of times candidate with index i is ranked
    
    vector<int> candidateSet;
    numberRanked.clear();
    numberRanked.resize(0);
    candidateSet.clear();
    candidateSet.resize(0);
    
    // Construct a set of candidates
    
    size_t numberPositions = rankingsOfPositions.size();
    
    for (j = 0 ; j < numberPositions ; j++) {
        size_t sizeRanking = rankingsOfPositions[j].size();
        for (i = 1 ; i < sizeRanking ; i++) {
            found=0;
            int theGuy = rankingsOfPositions[j][i];
            if (!presentInVector(candidateSet, theGuy)) {
                candidateSet.push_back(theGuy);
            }
        }
    }
    
    // For each candidate, list the positions for which he/she is ranked.
    //
    // Note that if both candidates i1 and i2 are ranked by positions j1 and j2, then j1 and j2 will apear in the same order in the candidates' list of positions.
    //
    // candidateSet[i] = the name of the i-th candidate in candidateSet[]
    // numberRanked[i] = a vector listing all the position for which the i-th candidate in candidateSet[] is ranked
    // numberRanked[i][h] = k means that for the i-th candidate, the h-th time he/she is ranked is by the k-th position (k-th in originalRankings).
    //
    
    
    size_t sizeCandidateSet = candidateSet.size();
    numberRanked.resize(sizeCandidateSet);
    for (i = 0 ; i < sizeCandidateSet ; i++) {
        int theGuy;
        for (j = 0 ; j < numberPositions ; j++) {
            if (indexInRankings(rankingsOfPositions[j], theGuy)>0) {
                numberRanked[i].push_back(j);
            }
        }
    }
    
    // Get the highest number of times a candidate is ranked
    maxSize=0;
    for (i = 0 ; i < sizeCandidateSet; i++) {
        if (numberRanked[i].size()>maxSize) {
            maxSize=(int)numberRanked[i].size();
        }
    }
    
    
    // DidTruncate will serve for the return of the function.
    DidTruncate=0;
    
    for (int k = 2 ; k < maxSize ; k++) {
        // for each possible number of times a candidate is ranked multiple times.
        bool foundAGroup;
        vector<int> theCandidates; // a set where we store the candidates ranked for the same positions
        for (i = 0 ; i < sizeCandidateSet-1 ; i++) {
            // Don't need to look at all candidate: we will compare the positions ranked
            // for the i-th candidate with those of candidates with a higher index.
            theCandidates.clear();
            theCandidates.resize(0);
            foundAGroup = false;
            if (numberRanked[i].size()==k) {
                // found a candidate ranked k times
                // add him/her to the set theCandidates[]
                theCandidates.push_back(candidateSet[i]);
                for (ii = i+1 ; ii < sizeCandidateSet ; ii++) {
                    // Look at candidates with a higher index
                    // add them to theCandidates[] each time they are ranked for the same positions as the i-th candidates
                    if (numberRanked[i]==numberRanked[ii]) {
                        // numberRanked[i] and numberRanked[ii] are VECTORS (the list of departments that rank those candidates).
                        theCandidates.push_back(candidateSet[ii]);
                        // Candidates NumerRanked[i] and NumerRanked[ii] ranked for the same positions
                    }
                    if (theCandidates.size()==k) {
                        // as soon as we have found k such candidates, stop the loop
                        // we found a group of k candidates ranked k times for the same k positions.
                        foundAGroup = true;
                        break;
                    }
                }
            }
            if (foundAGroup) {
                for (int h = 0 ; h < k ; h++) {
                    // for each of the k positions, find the rank of the lowest ranked candidate in those k candidates
                    int theLowestRank=-1;
                    int thePosition=numberRanked[i][h]; // the index of the position in originalRankings
                    int theGuy;
                    size_t theCandidatesSize = theCandidates.size();
                    for (ii = 0 ; ii < theCandidatesSize ; ii++) {
                        theGuy=theCandidates[ii]; // the name of the candidate -- just to simplify the code, the compiler will get rid of it.
                        // check if theGuy's rank is below the value of theLowestRank. If yes, then update the value of theLowestRank
                        int theIndex = indexInRankings(rankingsOfPositions[thePosition], theGuy);
                        if (theIndex>theLowestRank) {
                            theLowestRank=theIndex;
                        }
                    }
                    if (theLowestRank<rankingsOfPositions[thePosition].size()-1) {
                        for (int ii = theLowestRank+1 ; ii < rankingsOfPositions[thePosition].size() ; ii++) {
                            originalImpossible[thePosition][ii]=1;
                        }
                    }
                }
            }
        }
    }
}



//                                                                      //
//     Truncations at candidates ranked multiple times                  //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     General loop to eliminate impossible matches from simple blocks  //
//                                                                      //
//                                                                      //

//
//  This script combines the subroutines rankedOnce() and rankedMutiple()
//

void DUMMY::singletonBlocks(vector<vector<int> > &rankingsOfPositions) {
    theFunction="singletonBlocks";
    bool NotDoneYet = true;
    int TempValue=0;
    //
    //  Keep running as long as rankedOnce() performed a truncation.
    //  Stop when there's nothing left to truncate.
    //
    while (NotDoneYet) {
        TempValue=newrankedOnce(rankingsOfPositions);
        if (TempValue==0) {
            NotDoneYet = false;
        }
    }
}

//                                                                      //
//     General loop to eliminate impossible matches from simple blocks  //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Load and prepare the data                                        //
//                                                                      //
//                                                                      //


void DUMMY::loadData(string dataFile) {
    theFunction="loadData";
    
    time_t trand;    // time variables to count the time spent running the program
    
    (void) time(&trand);
    
    //
    //  Reset all variables.
    //
    originalRankings.clear();
    originalRankings.resize(0);
    originalImpossible.clear();
    originalImpossible.resize(0);
    vector<int> Candidates;
    vector<int> Positions;
    Candidates.clear();
    Candidates.resize(0);
    Positions.clear();
    Positions.resize(0);
    vector<vector<int> > OpePostes; // the data from the file
    OpePostes.clear();
    OpePostes.resize(0);
    int TempRanking[30]; // used to create originalRankings from OpePostes
    for (int i = 0 ; i < 30 ; i++) {
        TempRanking[i]=-1;
    }
    ifstream data_read;
    
    //
    //  Load the file and store it in OpePostes[][]
    //
    //  The file has 4 columns. Each row corresponds to a candidate:
    //  1st column: the candidate id
    //  2nd column: 1/0 Whether the candidate is assigned the position (by the Ministry)
    //  3rd column: rank of the candidate for that position
    //  4th column: id of the position
    
    vector<int> rowVector(4);
    int row = 0;
    data_read.open(dataFile.c_str(), ios::in);
    if (data_read.is_open()) {
        while (data_read.good()) {
            OpePostes.push_back(rowVector);
            for (int col = 0 ; col < 4 ; col++) {
                data_read >> OpePostes[row][col];
            }
            row++;
        }
    } else {
        cout << "Unable to opend file\n";
    }
    cout << "There are " << OpePostes.size() << " lines \n";
    cout << "Data read done\n";
    
    for (int j = 0 ; j < OpePostes.size()-1; j++) {
        //  for each position id (in OpePostes[][3]
        //  check if already in Position[]. If not add it.
        if (!presentInVector(Positions, OpePostes[j][3])) {
            Positions.push_back(OpePostes[j][3]);
        }
    }
    cout << Positions.size() << " positions, ";
    if (year>1998) {
        output.open (fileNameResults.c_str(), std::ofstream::out | std::ofstream::app);
        output << Positions.size()-1 << " positions, ";
        output.close();
    }
    hired.clear();
    hired.resize(Positions.size());
    for (int j = 0 ; j < hired.size() ; j++) {
        hired[j].resize(0);
        hired[j].push_back(Positions[j]);
        hired[j].push_back(-1);
    }
    for (int j = 0 ; j < Positions.size() ; j++) {
        for (int jj = 0 ; jj < OpePostes.size() ; jj++) {
            if (OpePostes[jj][3]==Positions[j]) {
                if (OpePostes[jj][1]==1) {
                    hired[j][1]=OpePostes[jj][0];
                }
            }
        }
    }
    //
    //  Construct now the originalRankings[][]
    //
    
    originalRankings.resize((int)(Positions.size()));
    for (int h = 0 ; h < Positions.size() ; h++) {
        // For each position, put in original_ranking[][0] the id of the position
        originalRankings[h].push_back(Positions[h]);
        for (int i = 0 ; i < 30 ; i++) {
            TempRanking[i]=-1;
        }
        // Look for all candidates ranked by Position[row]
        for (int j = 0 ; j < OpePostes.size()-1 ; j++) {
            // if there's a candidate ranked, his/her rank is OpePoste[][2]
            // construct a temporary ranking TempRanking by putting that candidate (with ID OpePostes[][0]
            // to the OpePostes[][2]-th position
            if (OpePostes[j][3]==Positions[h]) {
                TempRanking[OpePostes[j][2]]=OpePostes[j][0];
            }
        }
        // Now that we have all the candidates ranked by the position
        // construct the original_ranking[][] entry for that position
        // by adding the candidates in the order of their rank for that position.
        for (int i = 1 ; i < 30 ; i++) {
            if (TempRanking[i]!=-1) {
                originalRankings[h].push_back(TempRanking[i]);
            }
        }
    }
    
    vector<int> distributionSizeRankings;
    distributionSizeRankings.clear();
    distributionSizeRankings.resize(20);
    for (int j = 0 ; j < distributionSizeRankings.size() ; j++) {
        distributionSizeRankings[j]=0;
    }
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        distributionSizeRankings[originalRankings[j].size()-1]=distributionSizeRankings[originalRankings[j].size()-1]+1;
    }
    
    //
    //  Counting number of candidates and distribution of how many times ranked.
    //
    
    int numberViolations=0;
    vector <int> nonAcceptable;
    nonAcceptable.clear();
    nonAcceptable.resize(0);
    if (year>1998) {
        bool nothired=false;
        for (int j = 0 ; j < originalRankings.size() ; j++) {
            if (hired[j][1]!=-1) {
                if (indexInRankings(originalRankings[j], hired[j][1])>1) {
                    for (int i = 1 ; i < indexInRankings(originalRankings[j], hired[j][1]); i++) {
                        // Check if they are hired somewhere
                        nothired=true;
                        for (int jj = 0 ; jj < hired.size() ; jj++) {
                            if (hired[jj][1]==originalRankings[j][i]) {
                                nothired=false;
                            }
                        }
                        if (nothired) {
                            if (!presentInVector(nonAcceptable, originalRankings[j][i])) {
                                nonAcceptable.push_back(originalRankings[j][i]);
                            }
                            originalRankings[j].erase(originalRankings[j].begin()+i);
                            i = i-1;
                            numberViolations++;
                        }
                    }
                }
            } else {
                cout << "Position " << hired[j][0] << " not hiring anybody\n";
                Truncate(originalRankings[j], originalRankings[j][1]);
            }
        }
        cout << "There are " << numberViolations << " violations of acceptable hypothesis\n";
        finalRestults.open("acceptability_violation.txt", std::ofstream::out | std::ofstream::app);
        finalRestults << year << "\t"  << numberViolations <<  "\t" << nonAcceptable.size() << "\n";
        finalRestults.close();
    }

    
    vector<int> candSet;
    candSet.clear();
    candSet.resize(0);
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        for (int i = 1 ; i < originalRankings[j].size() ; i++) {
            if (!presentInVector(candSet, originalRankings[j][i])) {
                candSet.push_back(originalRankings[j][i]);
            }
        }
    }
    totalNumberCandidates=(int)(candSet.size());
    if (year>1998) {
        cout << candSet.size() << " candidates\n";
        output.open (fileNameResults.c_str(), std::ofstream::out | std::ofstream::app);
        output << candSet.size() << " candidates\n";
        output.close();
    }
    finalRestults.open("final-results.txt", std::ofstream::out | std::ofstream::app);
    finalRestults << "\n" << year << "\t" << originalRankings.size() << "\t" << candSet.size() << "\t";
    finalRestults.close();
    
    vector<int> rankDistribution;
    int maxRanked;
    maxRanked=0;
    for (int i = 0 ; i < candSet.size() ; i++) {
        int counter=0;
        for (int j = 0 ; j < originalRankings.size() ; j++) {
            if(indexInRankings(originalRankings[j], candSet[i])>0) {
                counter++;
            }
        }
        maxRanked=max(counter,maxRanked);
    }
    rankDistribution.clear();
    rankDistribution.resize(maxRanked+1);
    for (int i = 0 ; i < rankDistribution.size() ; i++) {
        rankDistribution[i]=0;
    }
    
    int counter;
    for (int i = 0 ; i < candSet.size() ; i++) {
        counter=0;
        for (int j = 0 ; j < originalRankings.size() ; j++) {
            if(indexInRankings(originalRankings[j], candSet[i])>0) {
                counter++;
            }
        }
        rankDistribution[counter]=rankDistribution[counter]+1;
    }
    if (year>1998) {
        finalRestults.open("ranking-distributions.txt", std::ofstream::out | std::ofstream::app);
        cout << "\nDistribution of number of times being ranked\n";
        finalRestults << year << "\t";
        for (int i = 1 ; i < rankDistribution.size() ; i++) {
            cout <<  i  << "\t";
        }
        cout << "\n";
        for (int i = 1 ; i < rankDistribution.size() ; i++) {
            cout << (float)(rankDistribution[i])/(float)(candSet.size())*100 << "%\t";
            finalRestults << (float)(rankDistribution[i])/(float)(candSet.size())*100 << "%\t";
        }
        cout << "\n";
        finalRestults << "\n";
        for (int i = 1 ; i < rankDistribution.size() ; i++) {
            cout <<  rankDistribution[i] << "\t";
        }
        cout << "\n\n";
        output.open (fileNameResults.c_str(), std::ofstream::out | std::ofstream::app);
        output << "\nDistribution of number of times being ranked\n";
        for (int i = 1 ; i < rankDistribution.size() ; i++) {
            output <<  i << "\t" ;
        }
        output << "\n";
        for (int i = 1 ; i < rankDistribution.size() ; i++) {
            output << (float)(rankDistribution[i])/(float)(candSet.size())*100 << "%\t";
        }
        output << "\n";
        for (int i = 1 ; i < rankDistribution.size() ; i++) {
            output << rankDistribution[i] << "\t";
        }
        output << "\n\n";
        output.close();
        finalRestults.close();
        statisticsHired();
    }
    candSet.clear();
    rankDistribution.clear();
}


//                                                                      //
//     Load and prepare the data                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Construct random profile                                         //
//                                                                      //
//                                                                      //

//	For each year's parameters (# of positions, length of rankings, etc.)
//	Make long list (virtualCandidatesQueue) by concatenating rankings
//	Random shuffle virtualCandidatesQueue
//	Fill dept's rankings using the order in virtualCandidatesQueue ("de-concatenate")
//	Then check no candidate appears twice (or more) in the same ranking


//vector<vector<int> > DUMMY::randomProfile(vector<int> virtualCandidatesQueue) {
void DUMMY::randomProfile(vector<int> virtualCandidatesQueue, vector<vector<int> > &originalRankings) {
    time_t trand;    // time variables to count the time spent running the program
    (void) time(&trand);
    std::srand (unsigned (trand));
    std::vector<int> myvector;
    random_device rd;
    mt19937 g(rd());
    bool OK_Profile = false;
	size_t sizeRankings = originalRankings.size();

    while (!OK_Profile) {
        shuffle(virtualCandidatesQueue.begin(), virtualCandidatesQueue.end(), g);
        int counterCand = 0 ;
        int counterDep = 0;
        for (int j = 0 ; j < sizeRankings ; j++) {
            size_t sizeOfRanking = originalRankings[j].size();
            for (int i = 1 ; i < sizeOfRanking ; i++) {
                originalRankings[counterDep][i]=virtualCandidatesQueue[counterCand];
                counterCand++;
            }
            counterDep++;
            bool hasDuplicate = false;
            for (int i = 1 ; i < sizeOfRanking-1 ; i++) {
                int theGuy = originalRankings[j][i];
                for (int ii = i+1 ; ii < sizeOfRanking ; ii++) {
                    if (theGuy==originalRankings[j][ii]) {
                        hasDuplicate=true;
                    }
                }
            }
            if (hasDuplicate) {
                break;
            } else {
                if (counterCand==sizeOfRanking-1) {
                    OK_Profile = true;
                }
            }
        }
    }
    cout << "profile constructed\t";
}

//                                                                      //
//     Construct random profile                                         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Reshuffle original rankings                                      //
//                                                                      //
//                                                                      //

void DUMMY::reshuffleRankings(const vector<vector<int> > &positionsLinked) {
    time_t trand;    // time variables to count the time spent running the program
    (void) time(&trand);
    std::srand (unsigned (trand));
    random_device rd;
    mt19937 g(rd());
	size_t sizeRankings = originalRankings.size();
	size_t numberLinkages = positionsLinked.size();
	vector<int> positionDone;
	positionDone.clear();
	positionDone.resize(sizeRankings,0);
	for (int j = 0 ; j < sizeRankings ; j++) {
        if (positionDone[j]==0) {
            shuffle(originalRankings[j].begin()+1, originalRankings[j].end(), g);
            for (int h = 0 ; h < numberLinkages-1 ; h++) {
                if (presentInVector(positionsLinked[h], j)) {
                    for (int hh = 0 ; hh < positionsLinked[h].size() ; hh++) {
                        if (positionsLinked[h][hh]!=j) {
							size_t sizeOfRanking = originalRankings[j].size();
                            for (int i = 1 ; i < sizeOfRanking ; i++) {
                                originalRankings[positionsLinked[h][hh]][i]=originalRankings[j][i];
                            }
                            positionDone[positionsLinked[h][hh]]=1;
                        }
                    }
                }
            }
        }
    }
}

//                                                                      //
//     Reshuffle original rankings                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Step One: find impossible matches                                //
//                                                                      //
//                                                                      //

//
//  There's a little bit more data preparation, then searching
//  for the impossible matches
//


void DUMMY::stepOne(int year_index) {
    theFunction="stepOne";
    
    theStep=1;
    
    ofstream fileRankings;
    ifstream dataRankings;
    
    
    // first simplify the ranking profile by eliminating
    // the impossible matches because of simple blocks.
    
    //
    //  Construct set of candidates
    //
    
    vector<int> lesCandidats;
    lesCandidats.clear();
    lesCandidats.resize(0);
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        for (int i = 1 ; i < originalRankings[j].size() ; i++) {
            if (!presentInVector(lesCandidats, originalRankings[j][i])) {
                lesCandidats.push_back(originalRankings[j][i]);
            }
        }
    }
    
    if (year_index!=-1) {
        cout << "There are " << lesCandidats.size() << " candidats\n";
        cout << "There are " << originalRankings.size() << " positions\n";
        cout << "Ratio = " << (float)lesCandidats.size()/(float)originalRankings.size() << " candidates per position\n";
        output.open (fileNameResults.c_str(), std::ofstream::out | std::ofstream::app);
        output << "There are " << lesCandidats.size() << " candidats\n";
        output << "There are " << originalRankings.size() << " positions\n";
        output << "Ratio = " << (float)lesCandidats.size()/(float)originalRankings.size() << " candidates per position\n";
        output.close();
    }
    
    bool skipCheckImpossible = false;
    if (year>1998) {
        cout << "check if need to run Step 1\n";
    }
    dataRankings.open(fileRankingsStepOne.c_str(), ios::in | ios::out);
    if(!dataRankings) { // file couldn't be opened
        skipCheckImpossible = false;
    } else {
        skipCheckImpossible = true;
        originalRankings = readRankings();
    }
    dataRankings.clear();
    dataRankings.close();
    
    singletonBlocks(originalRankings);
    
    
    
    
    //
    //  Declare impossibles (-1 = don't know, 0 = not impossible, 1 = impossible
    //
    originalImpossible.clear();
    originalImpossible.resize(originalRankings.size());
    int candidatesToGo = 0 ;
    for (int j = 0 ; j < originalRankings.size(); j++) {
        originalImpossible[j].push_back(originalRankings[j][0]);
        originalImpossible[j].push_back(0);
        for (int i = 2; i < originalRankings[j].size(); i++) {
            candidatesToGo++;
            originalImpossible[j].push_back(-1);
        }
    }

    newrankedMutiple(originalRankings);

    
    vector<vector<int> > linked;
    linked=linkedPositions(originalRankings);
    
    if (!skipCheckImpossible) {
        if (year_index!=-1) {
            cout << "\nComputing stepOne, stepOne ranking file not found\n";
            cout << "Checking Impossible for step One: " << candidatesToGo << " candidates to look at\n";
        }
        candidatesDone=1;
        //
        //	Find the impossible matches.
        //
        for (int j = 0 ; j < originalRankings.size() ; j++) {
            for (int i = 2 ; i < originalRankings[j].size() ; i++) {
                candidatesDone++;
                if (originalImpossible[j][i]==-1) {
					newCheckImpossible(originalRankings, originalRankings[j][i], j, originalImpossible);
                    // If the position is linked, no need to check for the other linked positions.
                    for (int jj = 0 ; jj < linked.size() ; jj++) {
                        if (presentInVector(linked[jj], j)) {
                            for (int h = 0 ; h < linked[jj].size() ; h++) {
                                originalImpossible[linked[jj][h]][i]=originalImpossible[j][i];
                            }
                        }
                    }
                }
            }
        }
        if (year_index!=-1) {
            cout << "\n";
            output << "\n";
        }
        // Eliminates impossible matches from the rankings
        int numberImpossibles = 0;
        for (int j = 0 ; j < originalRankings.size() ; j++) {
            for (int i = 2 ; i < originalRankings[j].size() ; i++) {
                if (originalImpossible[j][i]==1) {
                    numberImpossibles++;
                    originalRankings[j].erase(originalRankings[j].begin()+i);
                    originalImpossible[j].erase(originalImpossible[j].begin()+i);
                    i=i-1;
                }
            }
        }
        if (year>1998) {
            cout << "\nFound " << numberImpossibles << " additional impossible matches\n";
            output << "\nFound " << numberImpossibles << " impossible matches\n";
            fileRankings.open(fileRankingsStepOne.c_str(), std::ofstream::out);
            for (int j = 0 ; j < originalRankings.size() ; j++) {
                for (int i = 0 ; i < originalRankings[j].size() ; i++) {
                    fileRankings << originalRankings[j][i] << "\t";
                }
                fileRankings << "\n";
            }
            fileRankings.close();
        } else {
            
        }
    }
    int hiredCorrect = 0 ;
    int notHiredCorrect = 0;
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        // Look first at positions potentially hiring someone (rank at least one candidate)
        if (originalRankings[j].size()>1) {
            int count=0;
            for (int jj = 0 ; jj < originalRankings.size() ; jj++) {
                if (originalRankings[jj].size()>1) {
                    for (int i = 1 ; i < originalRankings[jj].size() ; i++) {
                        if (originalRankings[jj][i]==originalRankings[j][1]) {
                            // candidate ranked first at position j is ranked.
                            count++;
                        }
                    }
                }
            }
            if (count==1) {
                // ranked only once.
                // Now check whether the prediction is correct.
                for (int jj = 0 ; jj < hired.size() ; jj++) {
                    if (hired[jj][0]==originalRankings[j][0]) {
                        if (hired[jj][1]==originalRankings[j][1]) {
                            hiredCorrect++;
                        } else {
                            if (year>1998) {
                                cout << "Concours " << originalRankings[j][0] << " predicted = " << originalRankings[j][1] << " --- hired = " << hired[jj][1] << "\n";
                            }
                        }
                    }
                }
            }
        } else {
            // position does not rank anybody
            // Check that it is indeed the case
            for (int jj = 0 ; jj < hired.size() ; jj++) {
                if (hired[jj][0]==originalRankings[j][0]) {
                    if (hired[jj][1]==-1) {
                        notHiredCorrect++;
                    } else {
                        if (year > 1998) {
                            cout << "Concours " << originalRankings[j][0] << " predicted no hire but hired " << hired[jj][1] << "\n";
                        }
                    }
                }
            }
        }
    }
    int numberRealChoices = 0;
    int numberPseudoStars = 0;
    predictions(originalRankings, "After complex impossible matches", numberRealChoices, numberPseudoStars);
    
}
//                                                                      //
//     Step One: find impossible matches                                //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//          Find positions with identical rankings                      //
//                                                                      //
//                                                                      //

//  If 2 rankings are identical they are declared as "linked".
//  Note: for the French market linked positions have necessarily the same ranking
//  but 2 positions with the same rankings are not necessarily linked (i.e., recruiting
//  committees could have chosen different rankings
//  for computational purposes, linked position <=> identical rankings

vector<vector<int> > DUMMY::linkedPositions(vector<vector<int> > rankings) {
    vector<vector<int> > pseudoRankings;
    vector<vector<int> > linked;
    linked.clear();
    pseudoRankings.clear();
    pseudoRankings=rankings;
    int sizeOfLinked = 0;
    linked.resize(1);
    for (int j = 0 ; j < pseudoRankings.size() ; j++) {
        pseudoRankings[j].erase(pseudoRankings[j].begin()+0);
    }
    for (int j = 0 ; j < pseudoRankings.size()-1 ; j++) {
        vector<int > theLinked;
        theLinked.clear();
        theLinked.resize(0);
        bool found = false;
        for (int jj = j+1 ; jj < pseudoRankings.size() ; jj++) {
            if (pseudoRankings[j]==pseudoRankings[jj]) {
                //cout << "positions " << rankings[j][0] << " and " << rankings[jj][0] << " are linked\n";
                found=true;
                if (!presentInVector(theLinked, j)) {
                    theLinked.push_back(j);
                }
                if (!presentInVector(theLinked, jj)) {
                    theLinked.push_back(jj);
                }
                
            }
        }
        if (found) {
            bool already=false;
            for (int jj = 0 ; jj < linked.size() ; jj++) {
                for (int jjj = 0 ; jjj < theLinked.size() ; jjj++) {
                    if (presentInVector(linked[jj], theLinked[jjj])) {
                        already=true;
                    }
                }
            }
            if (!already) {
                for (int jj = 0 ; jj < theLinked.size() ; jj++) {
                    linked[sizeOfLinked].push_back(theLinked[jj]);
                }
                sizeOfLinked=sizeOfLinked+1;
                linked.resize(sizeOfLinked+1);
            }
        }
    }
    return linked;
}

//                                                                      //
//          Find positions with identical rankings                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Get the stars and their choice list                              //
//                                                                      //
//                                                                      //

void DUMMY::getStars(vector<vector<int> > profile, vector<int> &theStars, vector<vector<int> > &theirChoices) {
    theFunction="getStars";
    theStars.clear();
    theStars.resize(0);
    theirChoices.clear();
    theirChoices.resize(0);

    int numberTimesTop,numberTimesRanked;
    
    vector<vector<int> > linked;
    linked.clear();
    linked=linkedPositions(profile);
    linked.resize(linked.size()-1);

    // Make the list of candidates (marketStars[]) ranked several times and always 1st.
    for (int j = 0 ; j < profile.size() ; j++) {
        if (profile[j].size()>1) {
            numberTimesTop = 0;
            numberTimesRanked = 0;
            for (int jj = 0 ; jj < profile.size() ; jj++) {
                if (profile[jj].size()>1) {
                    if (profile[j][1]==profile[jj][1]) {
                        numberTimesTop=numberTimesTop+1;
                    }
                    for (int i = 1 ; i < profile[jj].size(); i++) {
                        if (profile[j][1]==profile[jj][i]) {
                            numberTimesRanked=numberTimesRanked+1;
                        }
                    }
                }
            }
            if (numberTimesTop==numberTimesRanked) {
                if (numberTimesTop>1) {
                    if (!presentInVector(theStars, profile[j][1])) {
                        theStars.push_back(profile[j][1]);
                    }
                }
            }
        }
    }
    
    //  For each candidate in marketStars[], make the list of
    //  the positions for which they are ranked 1st
    theirChoices.resize(theStars.size());
    for (int i = 0 ; i < theStars.size() ; i++) {
        vector<int> usedRows;
        usedRows.clear();
        usedRows.resize(0);
        for (int j = 0 ; j < profile.size() ; j++) {
            if (profile[j].size()>1) {
                if (profile[j][1]==theStars[i]) {
                    int linkRow=-1;
                    for (int h = 0 ; h < linked.size() ; h++) {
                        if (presentInVector(linked[h], j)) {
                            linkRow=h;
                        }
                    }
                    if (linkRow==-1) {
                        theirChoices[i].push_back(j);
                    } else {
                        if (!presentInVector(usedRows, linkRow)) {
                            // Haven't added any positions in linked[linkRow] yet.
                            // position #j is added to the choices
                            theirChoices[i].push_back(j);
                            // we won't consider anymore any position in linked[linkRow]
                            usedRows.push_back(linkRow);
                        }
                    }
                }
            }
        }
    }
    
}
//                                                                      //
//     Get the stars and their choice list                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////
//                                                                      //
//      Use market stars's observed choices to find more                //
//      impossible matches                                              //
//                                                                      //
//                                                                      //
void DUMMY::stepTwoReal() {
    theFunction="stepTwoReal";
    
    theStep=2;
    
    vector<int> theLoops;
    vector<string> thePeople;
    theLoops.clear();
    theLoops.resize(0);
    float numberPositions;
    numberPositions = (float)(originalRankings.size());
    
    maxDim=0;
    
    cout << "Starting second step with stars (Real choices)\n\n";
    
    int iter = 1;
    numberPositionsCleared=0;
    
    
    int numberRealChoices = 0;
    int numberPseudoStars = 0;
    int numberStars = 0;
    
    matchingStarsRealChoices(originalRankings, iter, numberRealChoices, numberPseudoStars, numberStars);
    
    finalRestults.open("marketStars.txt", std::ofstream::out | std::ofstream::app);
    finalRestults << year << "\t"  << numberRealChoices << "\t" << numberStars << "\t" << numberPseudoStars << "\n";
    finalRestults.close();

  
    predictions(originalRankings, "After Step Two using Real Choices", numberRealChoices, numberPseudoStars);
    
    finalRestults.open("final-results.txt", std::ofstream::out | std::ofstream::app);
    //finalRestults << "\n";
    finalRestults.close();
}
//                                                                      //
//      Use market stars's observed choices to find more                //
//      impossible matches                                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//          Assigning marketStars to a position                         //
//                                                                      //
//                                                                      //

// This is the "findAll" that needs to be recursive
// Remember: get profile and dim as inputs

void DUMMY::matchingStarsRealChoices(vector<vector<int> > &thisRankings, int dim, int &numberRealChoices, int &numberPseudoStars, int &numberStars) {
    //    void DUMMY::matchingStarsSimulated(vector<vector<int> > profile, int dim, vector<int> &theLoops, int oldDim, vector<string> &thePeople) {
    //theFunction="matchingStars";
    bool anotherLoop = false;
    
    if (dim>maxDim) {
        maxDim=dim;
    }
    
    //  Get the stars and implement their choices
    implementStarsRealChoices(thisRankings, numberRealChoices, numberStars);
    //  Findthere and eliminate impossible matches (if any)
    
    
    anotherLoop = findAndCleanImpossible(thisRankings);
    
    
    if (anotherStarsLoop(thisRankings)) {
        // There are new stars showing up in the updated rankings, do another loop
        matchingStarsRealChoices(thisRankings, dim+1, numberRealChoices, numberPseudoStars, numberStars);
    } else {
        int previousNumberRealChoices = numberRealChoices;
        implementPseudoStarsRealChoices(thisRankings, numberRealChoices, numberPseudoStars);
        if (previousNumberRealChoices!=numberRealChoices) {
            // Some pseudo stars have been found
            // First clean the impossible
            anotherLoop = findAndCleanImpossible(thisRankings);
            // Then do another loop
            matchingStarsRealChoices(thisRankings, dim+1, numberRealChoices, numberPseudoStars, numberStars);
        } else {
            cout << "\n\nFinal dim = " << dim << "\n";
            vector<int> newCleared;
            newCleared.clear();
            newCleared = marketCleared(thisRankings);
            cout << newCleared[0] + newCleared[1] << " positions cleared\n";
            cout << "Total number of observed choices implemented = " << numberRealChoices << "\n";
            numberPositionsCleared=numberPositionsCleared+newCleared[0] + newCleared[1];
            cout << "\n\nEND OF REAL CHOICES LOOPS\n";
        }
    }
}



//                                                                      //
//          Assigning marketStars to a position                         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//          Check if candidate-position is an impossible match          //
//                                                                      //
//                                                                      //


void DUMMY::newCheckImpossible(const vector<vector<int> > &profile, const int &i0, const int &s0, vector<vector<int> > &impossible) {
    
    size_t sizeProfile = profile.size();
    
    vector<int> matchingDep;
    matchingDep.clear();
    matchingDep.resize(0);
    vector<int> myMatchingDep;
    myMatchingDep.clear();
    myMatchingDep.resize(0);
    
    matchingDep.resize(sizeProfile,-1);
    myMatchingDep.resize(sizeProfile,-1);
    
    matchingDep[s0]=i0;
    bool isAnImpossible=true;
    decided = false;
    vector<int> candToBeMatched;
    candToBeMatched.clear();
    int indexOfi0 = indexInRankings(profile[s0], i0);
    for (int i = 1 ; i < indexOfi0 ; i++) {
        candToBeMatched.push_back(profile[s0][i]);
    }
    
    matchProcess(profile, matchingDep, isAnImpossible, myMatchingDep, candToBeMatched);

	
    if (isAnImpossible) {
        impossible[s0][indexInRankings(profile[s0], i0)]=1;
        for (int j = 0 ; j < sizeProfile ; j++) {
            if (myMatchingDep[j]!=-1) {
                if (indexInRankings(profile[j], myMatchingDep[j])>-1) {
                    impossible[j][indexInRankings(profile[j], myMatchingDep[j])]=0;
                } else {
                    cout << "Error, candidate " << myMatchingDep[j] << " matched to position " << profile[j][0] << " but not in ranking of position " << profile[j][0] << "\n";
                    cout << "Matching:\n";
                    for (int h = 0 ; h < sizeProfile ; h++) {
                        if (myMatchingDep[h]!=-1) {
                            cout << profile[h][0] << "\tmatched to\t" << myMatchingDep[h]<< "\n";
                        }
                    }
                    cout << "\n";
                    cout << "i0 = " << i0 << "\n";
                    printRankings(profile);
                }
            }
        }
    }
    
}

//                                                                      //
//          Check if candidate-position is an impossible match          //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Count # positions cleared                                        //
//                                                                      //
//                                                                      //

vector<int> DUMMY::marketCleared(vector<vector<int> > &profile) {
    theFunction="marketCleared";
    vector<int> numberClearedPositions;
    numberClearedPositions.clear();
    numberClearedPositions.resize(2);
    numberClearedPositions[0]=0;
    numberClearedPositions[1]=0;
    int count;
    //int counter;
    // counter=0;
    size_t sizeProfile = profile.size();
    
    for (int j = 0 ; j < sizeProfile ; j++) {
        if (profile[j].size()==2) {
            // There are just 2 entries in originalRankings[j]:
            // originalRankings[j][0] is the position ID and originalRankings[j][1] is the candidate ranked 1st.
            count=0;
            // Now check that originalRankings[j][1] is ranked only once.
            for (int jj = 0 ; jj < sizeProfile ; jj++) {
                if (profile[j].size()>1) {
                    size_t sizeRanking = profile[jj].size();
                    for (int i = 1 ; i < sizeRanking ; i++) {
                        if (profile[jj][i]==profile[j][1]) {
                            count=count+1;
                        }
                    }
                }
            }
            if (count==1) {
                // originalRankings[j][1] only ranked once, so he/she gets that position
                // => one position solved.
                numberClearedPositions[0]=numberClearedPositions[0]+1;
                //counter=counter+1;
            }
        }
        if (profile[j].size()==1) {
            numberClearedPositions[1]=numberClearedPositions[1]+1;
            //counter++;
        }
    }
    return numberClearedPositions;
}


//                                                                      //
//     Count # positions cleared                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Matching pseudo stars with observed choices                      //
//                                                                      //
//                                                                      //


//
//  A pseudo star is a candidate ranked alway weakly above to the hired candidate
//  but not necessarily always ranked first. Stability constraint implies that
//  the position that hires him/her is the most preferred.
//


void DUMMY::implementPseudoStarsRealChoices(vector<vector<int> > &profile, int &numberRealChoices, int &numberPseudoStars) {
    vector<int> thePseudoStars;
    thePseudoStars.clear();
    thePseudoStars.resize(0);
    
    vector<vector<int> > candidates;
    candidates.clear();
    candidates.resize(2);
    // Construct the set of releventcandidates: ranked at least twice.
    // candidates[0] = vector of candidates
    // candidates[1] = for each candidate the number of times he's ranked
    size_t sizeProfile = profile.size();
    
    for (int j = 0 ; j < sizeProfile ; j++) {
        if (profile[j].size()>1) {
            size_t sizeTheProfile = profile[j].size();
            for (int i = 1 ; i < sizeTheProfile ; i++) {
                int count=0;
                for (int jj = 0 ; jj < sizeProfile ; jj++) {
                    if (profile[jj].size()>1) {
                        size_t sizeTheOtherProfile = profile[jj].size();
                        for (int ii = 1 ; ii < sizeTheOtherProfile ; ii++) {
                            if (profile[jj][ii]==profile[j][i]) {
                                count++;
                            }
                        }
                    }
                }
                if (count>1) {
                    if (!presentInVector(candidates[0], profile[j][i])) {
                        candidates[0].push_back(profile[j][i]);
                        candidates[1].push_back(count);
                    }
                }
            }
        }
    }
    //
    //  For each candidate, count the number of times he is
    //  a potential blocker: ranked above the hired candidate
    //
    size_t sizeHired = hired.size();
    for (int i = 0 ; i < candidates[0].size() ; i++) {
        bool trueStar=false;
        int numberTimesBlocker=0;
        for (int j = 0 ; j < sizeProfile ; j++) {
            if (profile[j].size()>1) {
                if (indexInRankings(profile[j], candidates[0][i])>-1) {
                    // candidate present in the ranking.
                    // Now look for the hired candidate;
                    int hiredCandidate=-1;
                    for (int jj = 0 ; jj < sizeHired ; jj++) {
                        if (hired[jj][0]==profile[j][0]) {
                            hiredCandidate=hired[jj][1];
                            int indexHired = indexInRankings(profile[j], hiredCandidate); // ranked of the hired candidate
                            int indexTheGuy = indexInRankings(profile[j], candidates[0][i]); // rank of the candidate
                            if (indexHired>=indexTheGuy) {
                                // he's ranked weakly above (weakly because he could be the hired one).
                                numberTimesBlocker++;
                                if (indexTheGuy==1) {
                                    trueStar=true;
                                }
                            }
                        }
                    }
                }
            }
        }
        if (numberTimesBlocker==candidates[1][i]) {
            // he can block any position he's not assigned to
            // Add it to the list of pseudostars
            if (!presentInVector(thePseudoStars, candidates[0][i])) {
                thePseudoStars.push_back(candidates[0][i]);
            }
            if (!trueStar) {
                cout << "========================================================================Not a true star!\n";
            }
        }
    }
    numberPseudoStars = numberPseudoStars +     (int)(thePseudoStars.size());
    //numberRealChoices = numberRealChoices+(int)(thePseudoStars.size());
    
    //
    //  Now for each pseudostar, "match" him to the position that hired him
    //  match = erase all candidates ranked below him for the position that hired him
    //          and erase him to all other rankings where he may show up.
    //
    
    size_t currentNumberPseudoStars = thePseudoStars.size();
    //size_t numberHired = hired.size();
    //size_t numberProfile = profile.size();
    for (int i = 0 ; i < currentNumberPseudoStars ; i++) {
        for (int j = 0 ; j < sizeHired ; j++) {
            if (hired[j][1]==thePseudoStars[i]) {
                for (int jj = 0 ; jj < sizeProfile ; jj++) {
                    if (profile[jj][0]==hired[j][0]) {
                        int indexPseudoStar = indexInRankings(profile[jj], thePseudoStars[i]);
                        if (indexPseudoStar==1) {
                            numberRealChoices++;
                        }
                        if (profile[jj].size()>indexPseudoStar+1) {
                            // the ranking that "hires" him is profile[jj]
                            // truncate below
                            Truncate(profile[jj], profile[jj][indexInRankings(profile[jj], thePseudoStars[i])+1]);
                        }
                    } else {
                        if (profile[jj].size()>1) {
                            if (indexInRankings(profile[jj], thePseudoStars[i])!=-1) {
                                // We deduce that he/she prefers the position hired[j][0] (which hired him/her)
                                // to the position profile[jj][0]
                                preferences <<  thePseudoStars[i] << "\t" << hired[j][0] << "\t" << profile[jj][0] << "\n";
                                // He's not hired by profile[jj], so just delete him from the ranking.
                                profile[jj].erase(profile[jj].begin()+indexInRankings(profile[jj], thePseudoStars[i]));
                            }
                        }
                    }
                }
            }
        }
    }
 }


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Implement stars' choices                                         //
//                                                                      //
//                                                                      //

//  1. Find the stars (ranked more than once, and always first)
//  2. Get the position that hires them: it's the most preferred
//  3. Store in file the deduced partial preferences (hiring position preferred to all other ones)
//  4. Eliminate stars from the rankings of the positions not hiring them
//  5. Eliminate all candidates ranked 2nd and below from the position hiring a star


void DUMMY::implementStarsRealChoices(vector<vector<int> > &profile, int &numberRealChoices, int &numberStars) {
    // Implement choices
    //theFunction="implementStarsRealChoices";
    vector<int> theStars;
    theStars.clear();
    theStars.resize(0);
    vector<vector<int> > theirChoices;
    theirChoices.clear();
    theirChoices.resize(0);
    //int numberTimesTop,numberTimesRanked;
    bool foundStars;
    foundStars=false;
    
    
    
    
    getStars(profile, theStars, theirChoices);
    
    numberRealChoices = numberRealChoices+(int)(theStars.size());
    numberStars = numberStars + (int)(theStars.size());
    
    vector<vector<int> > linked;
    linked.clear();
    linked=linkedPositions(profile);
    linked.resize(linked.size()-1);

    
    if(theStars.size()>0) {
        for (int i = 0 ; i < theStars.size() ; i++) {
            // take a star, theStars[i]
            for (int j = 0 ; j < hired.size() ; j++) {
                // search the position hiring him
                if (hired[j][1]==theStars[i]) {
                    // position[j] hires him
                    for (int jj = 0 ; jj < profile.size() ; jj++) {
                        // look which position it is (in case positions' indices in profile[] and hired[] do not match
                        if (profile[jj][0]==hired[j][0]) {
                            int lengthRanking=-1;
                            lengthRanking = (int)(profile[jj].size());
                            int theIndex = -1;
                            int currentStart = -1;
                            currentStart = theStars[i];
                            theIndex = indexInRankings(profile[jj], currentStart);
                            if (lengthRanking>theIndex+1) {
                                int theTruncated = -1;
                                theTruncated = profile[jj][theIndex+1];
                                Truncate(profile[jj], theTruncated);
                            }
                        } else {
                            if (profile[jj].size()>1) {
                                // get theStars[i] rank in profile[jj][];
                                if (indexInRankings(profile[jj], theStars[i])>-1) {
                                    preferences << profile[jj][1] << "\t" << hired[j][0] << "\t" << profile[jj][0] << "\n";
                                    profile[jj].erase(profile[jj].begin()+indexInRankings(profile[jj], theStars[i]));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

//                                                                      //
//     Implement stars' choices                                         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Check if need for another run of implement real choices          //
//                                                                      //
//                                                                      //

bool DUMMY::anotherStarsLoop(vector<vector<int> > thisRankings) {
    
    int numberTimesRankedFirst,numberTimesRanked;
    numberTimesRankedFirst = 0;
    numberTimesRanked = 0;
    bool indeed;
    indeed=false;

    // Make the list of candidates (marketStars[]) ranked several times and always 1st.
    vector<int> lesStars;
    lesStars.clear();
    lesStars.resize(0);
    
    for (int j = 0 ; j < thisRankings.size() ; j++) {
        if (thisRankings[j].size()>1) {
            numberTimesRankedFirst = 0;
            numberTimesRanked = 0;
            for (int jj = 0 ; jj < thisRankings.size() ; jj++) {
                if (thisRankings[jj].size()>1) {
                    if (thisRankings[j][1]==thisRankings[jj][1]) {
                        numberTimesRankedFirst=numberTimesRankedFirst+1;
                    }
                    for (int i = 1 ; i < thisRankings[jj].size() ; i++) {
                        if (thisRankings[j][1]==thisRankings[jj][i]) {
                            numberTimesRanked=numberTimesRanked+1;
                        }
                    }
                }
            }
            if (numberTimesRankedFirst==numberTimesRanked) {
                if (numberTimesRanked>1) {
                    indeed=true;
                    if (!presentInVector(lesStars, thisRankings[j][1])) {
                        lesStars.push_back(thisRankings[j][1]);
                    }
                }
            }
        }
    }
    return indeed;
}

//                                                                      //
//     Check if need for another run of implement real choices          //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Load rankings obtained from stepOne                              //
//                                                                      //
//                                                                      //

vector<vector<int> > DUMMY::readRankings() {
    
    ifstream data;
    string line;
    vector<string> list;
    list.resize(0);
    vector<vector<string> > sRankings;
    vector<vector<int> > rankings;
    rankings.clear();
    rankings.resize(0);
    
    
    
    data.open(fileRankingsStepOne.c_str(), ios::in | ios::out);
    //data.open(fileRankingsStepOne, ios::in | ios::out);
    if(!data) { // file couldn't be opened
        cout << "Error, file " << fileRankingsStepOne << " not found\n";
        exit(1);
    } else {
        cout << "Skipping stepOne, stepOne ranking file found" << endl;
        size_t lines = 0;
        ifstream in(fileRankingsStepOne);
        for (string s; getline(in,s); ) {
            ++lines;
        }
        sRankings.resize(lines);
        rankings.resize(lines);
        
        while (std::getline(data, line)) {
            list.push_back(line);
        }
        for (int i = 0 ; i < list.size() ; i++) {
            string token;
            string mystring = list[i];
            while (token!=mystring) {
                token=mystring.substr(0,mystring.find_first_of("	"));
                mystring = mystring.substr(mystring.find_first_of("	") + 1);
                sRankings[i].push_back(token);
            }
            
        }
        for (int j = 0 ; j < sRankings.size() ; j++) {
            if (sRankings[j].size()>2) {
                for (int i = 0 ; i < sRankings[j].size()-1 ; i++) {
                    int candidate;
                    candidate = std::stoi(sRankings[j][i]);
                    rankings[j].push_back(candidate);
                }
            } else {
                int candidate;
                candidate = std::stoi(sRankings[j][0]);
                rankings[j].push_back(candidate);
            }
        }
    }
    data.clear();
    data.close();
    return rankings;
}

//                                                                      //
//     Load rankings obtained from stepOne                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//    Find an unmatched candidate                                       //
//                                                                      //
//                                                                      //

//  Given a matching, find a candidate who is not matched but ranked for
//  a position above a candidate matched (to that position)

void DUMMY::findCandidate(const vector<vector<int>> &thisrankings, const vector<int> &matchingDep, int &toBeMatched) {
    // get a candidate to be matched
    size_t numberRankings = thisrankings.size();
    for (int j = 0 ; j < numberRankings ; j++) {
        bool stop = false;
        if (matchingDep[j]!=-1) {
            int indexTheMatchedCandidate = indexInRankings(thisrankings[j], matchingDep[j]);
            for (int i = 1 ; i < indexTheMatchedCandidate ; i++) {
                // Now see if he's matched.
                if (!presentInVector(matchingDep, thisrankings[j][i])) {
                    toBeMatched=thisrankings[j][i];
                    stop = true;
                    break;
                }
            }
        }
        if (stop) {
            break;
        }
    }
}
//                                                                      //
//    Find an unmatched candidate                                       //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//    For an unmatched candidate get all available positions            //
//    (that rank that candidate)                                        //
//                                                                      //

void DUMMY::getChoices(const vector<vector<int>> & thisrankings, const vector<int> & matchingDep, int &toBeMatched, vector<int> &choiceSet, const vector<int> & candToBeMatched) {
    vector<int> possiblePositions;
    possiblePositions.clear();
    vector<int> possibleRanks;
    possibleRanks.clear();
    choiceSet.clear();
    
    
    
    vector<int> totalChoices;
    totalChoices.clear();
    totalChoices.resize(candToBeMatched.size(), 0);
    /*
     for (int i = 0 ; i < candToBeMatched.size() ; i++) {
     totalChoices[i]=0;
     }
     */
    size_t numbRankings = thisrankings.size();
    vector<int> requestedPositions;
    requestedPositions.clear();
    requestedPositions.resize(numbRankings, 0);
    for (int j = 0 ; j < numbRankings ; j++) {
        if (matchingDep[j]==-1) {
            size_t lengthRanking = thisrankings[j].size();
            for (int i = 1 ; i < lengthRanking ; i++) {
                int theGuy = thisrankings[j][i];
                int myIndex = index(candToBeMatched, theGuy);
                if (myIndex>-1) {
                    if (originalImpossible[j][i]!=1) {
                        totalChoices[myIndex]++;
                    }
                    requestedPositions[j]=1;
                }
            }
        }
    }
    
    int totalRequests = 0;
    for (int j = 0 ; j < numbRankings ; j++) {
        totalRequests = totalRequests + requestedPositions[j];
    }
    
    
    bool weStop = false;
    if (totalRequests<candToBeMatched.size()) {
        weStop = true;
    }
    size_t sizeChoices = totalChoices.size();
    if (!weStop) {
        for (int i = 0 ; i <  sizeChoices ; i++) {
            if (totalChoices[i]==0) {
                weStop=true;
                /*
                 
                 int total = 1;
                 for (int ii = 0 ; ii < totalChoices.size() ; ii++) {
                 if (totalChoices[ii]!=0) {
                 total=total*totalChoices[ii];
                 }
                 }
                 //cout << "Don't need to go further (with " << totalChoices.size() << " candidates remaining";
                 //cout << " and " << total << " combinations avoided!)\n";
                 */
                break;
            }
        }
    }
    toBeMatched=-1;
    int min=(int)(numbRankings);
    for (int i = 0 ; i < sizeChoices ; i++) {
        int numberChoices = totalChoices[i];
        if (numberChoices>0) {
            if (numberChoices<min) {
                min=totalChoices[i];
                toBeMatched=candToBeMatched[i];
            }
        }
    }
    if (!weStop) {
        //cout << "Don't stop for candidate " << matchingDep[0] << " for position " << thisrankings[0][0] << "\n";
        for (int j = 0 ; j < numbRankings ; j++) {
            if (matchingDep[j]==-1) {
                int myindex = indexInRankings(thisrankings[j], toBeMatched);
                if (myindex > -1) {
                    if (originalImpossible[j][myindex]!=1) {
                        possiblePositions.push_back(j);
                        //possibleRanks.push_back(indexInRankings(thisrankings[j], toBeMatched));
                        possibleRanks.push_back(myindex);
                        //choiceSet.push_back(j);
                    }
                }
            }
        }
        size_t sizePossibleRanks = possibleRanks.size();
        for (int hh = 0 ; hh < sizePossibleRanks ; hh++) {
            for (int h = 0 ; h < sizePossibleRanks-1 ; h++) {
                int pos,rank;
                if (possibleRanks[h+1]<possibleRanks[h]) {
                    pos=possiblePositions[h];
                    rank=possibleRanks[h];
                    possiblePositions[h]=possiblePositions[h+1];
                    possibleRanks[h]=possibleRanks[h+1];
                    possiblePositions[h+1]=pos;
                    possibleRanks[h+1]=rank;
                }
            }
        }
        choiceSet=possiblePositions;
    }
}


//                                                                      //
//    For an unmatched candidate get all available positions            //
//    (that rank that candidate)                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//    Match process                                                     //
//                                                                      //
//                                                                      //

//  Given a matching, find a candidate to be matched (not matched but ranked
//  above some matched candidate. Match the candidate to one of the possible (unmatched)
//  position ranking him (explore all possible choices), and then repeat.
//  Stop if there's nobody to be matched (toBeMatched = -1) of when all possibilities
//  have been explored (the boolean isAnImpossible has been initialized to "true" before
//  calling matchProcess.
//


void DUMMY::matchProcess(const vector<vector<int>> &thisrankings, const vector<int> & matchingDep, bool &isAnImpossible, vector<int> &myMatchingDep, const vector<int> & candToBeMatched) {

    vector<int> choiceSet;
    choiceSet.clear();
    choiceSet.resize(0);
    int toBeMatched=-1;

    // get a candidate
    findCandidate(thisrankings, matchingDep, toBeMatched);
    
    if (toBeMatched==-1) {
        isAnImpossible = false; // If nobody to be matched we're done. Not an impossible
        decided=true;   // that will stop the loop below: if (choiceSet.size()!=0)
        myMatchingDep = matchingDep;    // we'll use that matching to avoid computing further impossible matches
    }
    
    getChoices(thisrankings, matchingDep, toBeMatched, choiceSet, candToBeMatched);
	
    // get the choices of the candidate
    
    if (choiceSet.size()!=0) {
        // The candidate can be matched somewhere.
        vector<int> thisMatchingDep;
        //vector<int> thisMatchingCand;
        //thisMatchingCand=matchingCand;
        thisMatchingDep=matchingDep;
		vector<int> thiscandToBeMatched;
        thiscandToBeMatched.clear();
        thiscandToBeMatched.resize(0);
		bool ilestla=presentInVector(candToBeMatched, toBeMatched);
        // implement choices: for (choice = 0 ; choice.size()) {do...}
        for (int h = 0 ; h < choiceSet.size() ; h++) {
            if (!decided) {
                // toBeMatched no longer to be matched
				thiscandToBeMatched=candToBeMatched;
                if (ilestla) {
                    thiscandToBeMatched.erase(thiscandToBeMatched.begin()+index(thiscandToBeMatched, toBeMatched));
                }
                // Now add all the candidates not matched but ranked above toBeMatched
				int myend = indexInRankings(thisrankings[choiceSet[h]],toBeMatched);
                for (int ii = 1 ; ii <  myend; ii++) {
                    int aCandidate = thisrankings[choiceSet[h]][ii];
                    if (!presentInVector(matchingDep, aCandidate)) {
                        if (!presentInVector(thiscandToBeMatched, aCandidate)) {
                            thiscandToBeMatched.push_back(aCandidate);
                        }
                    }
                }
                if (h>0) {
                    thisMatchingDep[choiceSet[h-1]]=-1;
                }
                thisMatchingDep[choiceSet[h]]=toBeMatched;
                matchProcess(thisrankings, thisMatchingDep, isAnImpossible, myMatchingDep, thiscandToBeMatched);
            }
        }
    }
}


//                                                                      //
//    Match process                                                     //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//    Predictions                                                       //
//                                                                      //
//                                                                      //

//  Match of a position is predicted if either:
//  * only rank 1 candidate and that candidate is ranked only once
//  * ranks nobody
//  Then check whether the prediction corresponds to actual hires from the data.


void DUMMY::predictions(vector<vector<int> > &rankings, string step, int numberRealChoices, int numberPseudoStars) {
    
    vector<int> numberClearedPositions;
    numberClearedPositions.clear();
    numberClearedPositions.resize(2);
    numberClearedPositions[0]=0;
    numberClearedPositions[1]=0;
    int hiredCorrect = 0 ;
    int notHiredCorrect = 0;
    
    
    if (year>1998) {
        for (int j = 0 ; j < rankings.size() ; j++) {
            // Look first at positions potentially hiring someone (rank at least one candidate)
            if (rankings[j].size()>1) {
                int count=0;
                for (int jj = 0 ; jj < rankings.size() ; jj++) {
                    if (rankings[jj].size()>1) {
                        for (int i = 1 ; i < rankings[jj].size() ; i++) {
                            if (rankings[jj][i]==rankings[j][1]) {
                                // candidate ranked first at position j (=rankings[j][1]) is ranked.
                                count++;
                            }
                        }
                    }
                }
                if (count==1) {
                    if (rankings[j].size()>2) {
                        Truncate(rankings[j], rankings[j][2]);
                    }
                    // rankings[j][1] is ranked only once (and first by j-th position!)
                    // Now check whether the prediction is correct.
                    for (int jj = 0 ; jj < hired.size() ; jj++) {
                        if (hired[jj][0]==rankings[j][0]) {
                            if (hired[jj][1]==rankings[j][1]) {
                                hiredCorrect++;
                            } else {
                                cout << "year = " << year << "\n";
                                cout << "Not correct: ";
                                cout << "Prediction is " << rankings[j][1] << " but it should be " << hired[jj][1];
                                cout << " for position " << rankings[j][0] << "\n";
                                cout << "And " << rankings[j][1] << " is in fact hired by ";
                                bool foundHisMatch= false;
                                for (int h = 0 ; h < hired.size() ; h++) {
                                    if (hired[h][1]==rankings[j][1]) {
                                        foundHisMatch=true;
                                        cout << hired[h][0] << " (=" << rankings[h][0] << "? from rankings)\n";
                                    }
                                }
                                if (!foundHisMatch) {
                                    cout << "nobody\n";
                                }
                            }
                        }
                    }
                }
            } else {
                // position does not rank anybody
                // Check that it is indeed the case
                for (int jj = 0 ; jj < hired.size() ; jj++) {
                    if (hired[jj][0]==rankings[j][0]) {
                        if (hired[jj][1]==-1) {
                            notHiredCorrect++;
                        }
                    }
                }
            }
        }

    }
    
    numberClearedPositions = marketCleared(rankings);
    
    //
    //  Counting nunber of predicted candidates
    //
    vector<int> remainingCandidates;
    remainingCandidates.clear();
    remainingCandidates.resize(0);
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        for (int i = 1 ; i < originalRankings[j].size() ; i++) {
            if (!presentInVector(remainingCandidates, originalRankings[j][i])) {
                remainingCandidates.push_back(originalRankings[j][i]);
            }
        }
    }
    int numberPredictedCandidates;
    numberPredictedCandidates = totalNumberCandidates - (int)(remainingCandidates.size()) + (int)(numberClearedPositions[0]);
    
    
    
    
    if (year>1998) {
        output.open(fileNameResults.c_str(), ofstream::out | ofstream::app);
        cout << "\n" << step << "\n";
        cout << "There are " <<  numberClearedPositions[0] << " hiring positions cleared, " << hiredCorrect << " are correct\n";
        cout << "We found  " <<  numberClearedPositions[1] << " not hiring positions cleared, " << notHiredCorrect << " are correct\n\n";
        output << "\n" << step << "\n";
        output << "There are " <<  numberClearedPositions[0] << " hiring positions cleared, " << hiredCorrect << " are correct\n";
        output << "We found  " <<  numberClearedPositions[1] << " not hiring positions cleared, " << notHiredCorrect << " are correct\n\n";
        finalRestults.open("final-results.txt", std::ofstream::out | std::ofstream::app);
        if (theStep==1) {
            finalRestults << numberPredictedCandidates << "\t";
            finalRestults << numberClearedPositions[0] + numberClearedPositions[1] << "\t";
        } else {
            finalRestults << numberPredictedCandidates-numberRealChoices << "\t";
            finalRestults << numberClearedPositions[0] + numberClearedPositions[1] - numberRealChoices<< "\t";
        }
        finalRestults.close();
        if (numberClearedPositions[0]==hiredCorrect && numberClearedPositions[1]==notHiredCorrect) {
            cout << "******** 100% Correct prediction ********\n\n";
            output << "******** 100% Correct prediction ********\n\n";
        }
        output.close();
    } else {
        if (step.compare("Before complex impossible matches")!=0) {
            cout << "\t" << numberClearedPositions[0] << "\t" << numberClearedPositions[1] << "\n";
            howManyCleared = numberClearedPositions[0];
        }
    }
    
}

//                                                                      //
//    Predictions                                                       //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//    Find an Eliminate Impossible Matched                              //
//                                                                      //
//                                                                      //

bool DUMMY::findAndCleanImpossible(vector<vector<int>> &thisRankings) {
    bool foundImpossible=false;
    vector<vector<int> > thisImpossible;
    thisImpossible.clear();
    thisImpossible.resize(0);
    thisImpossible.clear();
    thisImpossible.resize(thisRankings.size());
    for (int j = 0 ; j < thisRankings.size(); j++) {
        thisImpossible[j].push_back(thisRankings[j][0]);
        if (thisRankings[j].size()>1) {
            thisImpossible[j].push_back(0);
            if (thisRankings[j].size()>2) {
                for (int i = 2; i < thisRankings[j].size(); i++) {
                    thisImpossible[j].push_back(-1);
                }
            }
        }
    }
    int numberCases=0;
    // Just for curiosity: how many "checkImpossible" we have to run
    for (int j = 0 ; j < thisRankings.size() ; j++) {
        if (thisRankings[j].size()>2) {
            numberCases=numberCases+(int)(thisRankings[j].size()-2);
        }
    }
    //cout << "Number of candidates to look at = " << numberCases <<"\n";
    
    int currentCase=1;
    // Find the impossible matches
    for (int j = 0 ; j < thisRankings.size() ; j++) {
        if (thisRankings[j].size()>2) {
            for (int i = 2 ; i < thisRankings[j].size() ; i++) {
                //cout << currentCase << "/" << numberCases << "\t";
                currentCase++;
                if (thisImpossible[j][i] == -1) {
                    newCheckImpossible(thisRankings, thisRankings[j][i], j, thisImpossible);
                }
            }
        }
    }
    // Eliminates the impossible matches and sets next profile
    for (int j = 0 ; j < thisRankings.size() ; j++) {
        if (thisRankings[j].size()>2) {
            for (int i = 2 ; i < thisRankings[j].size() ; i++) {
                if (thisImpossible[j][i]==1) {
                    foundImpossible = true;
                    thisRankings[j].erase(thisRankings[j].begin()+i);
                    thisImpossible[j].erase(thisImpossible[j].begin()+i);
                    i=i-1;
                }
            }
        }
    }
    return foundImpossible;
}

//                                                                      //
//    Find an Eliminate Impossible Matched                              //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//    Statistics about hired candidates                                 //
//                                                                      //
//                                                                      //

void DUMMY::statisticsHired() {
    vector<int> rankedHired;
    rankedHired.clear();
    rankedHired.resize(15,0);
    vector<int> ranksOfHired;
    ranksOfHired.clear();
    ranksOfHired.resize(15,0);
    
    int totalHires = 0;
    for (int j = 0 ; j < hired.size() ; j++) {
        if (hired[j][1]!=-1) {
            totalHires++;
            // Look for the position in original rankings;
            for (int jj = 0 ; jj < originalRankings.size() ; jj++) {
                if (originalRankings[jj][0]==hired[j][0]) {
                    for (int i = 1 ; i < originalRankings[jj].size() ; i++) {
                        if (originalRankings[jj][i]==hired[j][1]) {
                            rankedHired[i]++;
                        }
                    }
                }
            }
        }
    }
    
    finalRestults.open("Statistics-hired.txt", std::ofstream::out | std::ofstream::app);
    for (int i = 1 ; i < rankedHired.size() ; i++) {
        cout << i << "\t";
    }
    cout << "\n";
    for (int i = 1 ; i < rankedHired.size() ; i++) {
        cout << 100*((float)(rankedHired[i]))/((float)(totalHires)) << "%\t";
    }
    cout << "\n";
    finalRestults << year << "\t";
    for (int i = 1 ; i < rankedHired.size() ; i++) {
        finalRestults << 100*((float)(rankedHired[i]))/((float)(totalHires)) << "%\t";
    }
    finalRestults << "\n";
    finalRestults.close();

    for (int j = 0 ; j < hired.size() ; j++) {
        if (hired[j][1]!=-1) {
            int counter = 0;
            // search how many times he's ranked
            for (int jj = 0 ; jj < originalRankings.size() ; jj++) {
                if (indexInRankings(originalRankings[jj], hired[j][1])>-1) {
                    counter++;
                }
            }
            ranksOfHired[counter]=ranksOfHired[counter]+1;
        }
    }
    cout << "\n\nRank of Hired candidates\n";
    finalRestults.open("Statistics-hired-partII.txt", std::ofstream::out | std::ofstream::app);
    finalRestults << year << "\t";
    for (int i = 1 ; i < 15 ; i++) {
        cout << ranksOfHired[i] << "\t";
        finalRestults << ranksOfHired[i] << "\t";
    }
    cout << "\n";
    finalRestults << "\n";
    finalRestults.close();
}

//                                                                      //
//    Statistics about hired candidates                                 //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


void DUMMY::printRankings(vector<vector<int> > rankings) {
    cout << "\n\n";
    for (int j = 0 ; j < rankings.size() ; j++) {
        for (int i = 0 ; i < rankings[j].size() ; i++) {
            cout << rankings[j][i] << "\t";
        }
        cout << "\n";
    }
    cout << "\n";
}




