//
//  dummy.cpp
//  Impossible_matches
//
//  Created by Guillaume Haeringer on 4/28/14.
//  Copyright (c) 2014 Guillaume Haeringer. All rights reserved.
//

#include "dummy.h"

DUMMY::DUMMY()
{
    
}



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Some useful functions                                            //
//                                                                      //
//                                                                      //


// Give the index number of "a" in a vector "S"
int DUMMY::Index(const vector<int> &S, int a){
    int  ans=-1;
    bool isPresent = (std::find(S.begin(), S.end(), a) != S.end());
    if (isPresent) {
        for (int i = 0 ; i < S.size(); i++) {
            if (S[i]==a) {
                ans=i;
            }
        }
    }
    return(ans);
}


// Give the index number of "a" in a vector "S" except first entry
// rankings[j][0] is the ID # of position j, could be the same as that of a candidate

int DUMMY::indexInRankings(const vector<int> &S, int a){
    int  ans=-1;
    bool isPresent = (std::find(S.begin(), S.end(), a) != S.end());
    if (isPresent) {
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
    if (found==1) {
        if (Index(S,z)!=-1) {
            return 1;
        } else {
            return 0;
        }
    } else {
        return 0;
    }
}




// Given the ranking of a position (vector S), truncate at candidate "a"
// (i.e., a is deleted from the ranking and all the candidates below as well
// "a" is the ID (or NAME) of the candidate, not its rank number in the ranking of the department.

void DUMMY::Truncate(vector<int> &S, int a){
    S.erase(S.begin()+indexInRankings(S,a),S.end());
}



//                                                                      //
//     Some useful functions                                            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Maximum Matching sub-Routine                                     //
//                                                                      //
//                                                                      //


// find a maximum matching (with Depth First Search) given an adjency matrix
//  edge[][]
//


int DUMMY::dfs(int a) {
    if(a<0) return(1);
    if(visited[a]) return(0);
    visited[a]=1;
    int i;
    for(i=0;i<rankings.size();i++) if(edge[a][i]) //* see remark
    {
        if(dfs(matching_department[i]))
        {
            matching_candidate[a]=i;matching_department[i]=a;
            return(1);
        }
    }
    return(0);
}

int DUMMY::dfsExp(int a) {
    visited.resize(Candidates.size());
    int i;
    for(i=0;i<Candidates.size();i++) visited[i]=0;
    return dfs(a);
}

int DUMMY::bipMatch()
{
    int i;
    int ans=0;
    for(i=0;i<Candidates.size();i++) {
        if(matching_candidate[i]<0) ans+=dfsExp(i);
    }
    return(ans);
}

//                                                                      //
//     Maximum Matching sub-Routine                                     //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Construct a graph (edge[][]) for comprehensive matching          //
//                                                                      //
//                                                                      //


//  Given a comprehensive matching, make a graph such that for each position there is an edge
//  between the 1st ranked candidate and the position, the 2nd ranked candidate and the position, etc.
//  until we attain the first candidate who is not matched to any position: make also an edge for that candidate
//  but then stop.
//
//  this ensures that if we construct a maximum matching with that graph the resulting matching
//  is necessarily comprehensive



void DUMMY::admissibleGraph() {
    
    int i,j;
    unsigned long stop=1;
    
    // First reinitialize the graph
    edge.clear();
    edge.resize(Candidates.size());
    for (int i = 0 ; i < Candidates.size() ; i++) {
        for (int j = 0 ; j < rankings.size(); j++) {
            edge[i].push_back(0);
        }
    }
    // For each position j and the candidate ranked 1st in j's ranking, say, i (i.e., ranking[j][1]=i there's an edge.
    // (the edge is done using the index of ranking[j][1] in Candidates[].)
    for (j = 0 ; j < rankings.size() ; j++) {
        if (rankings[j][1]>-1) { // so rankings[j][1] = 0, 1, 2, ... => the 1st ranked is a candidate (0, 1, 2, ... is the name of the candidate
            edge[Index(Candidates, rankings[j][1])][j]=1; // edge created between that candidate and that department.
            // edge[i][j] = 1 => edge between the i-th candidate in Candidates and the j-th position.
        }
    }
    
    //
    //  For each position, search the highest ranked candidate that is not matched (to that position or to another position).
    //
    for (j = 0 ; j < rankings.size(); j++) {
        stop=rankings[j].size()-1;
        for (int i = 1 ; i < rankings[j].size(); i++) {
            if (matching_candidate[Index(Candidates, rankings[j][i])]==-1) {
                stop=i;
                break;
            }
        }
        // Once that candidate is identified, make an edge between this candidate and the position,
        // and the same for all better ranked candidate.
        for (i = 1 ; i < stop+1; i++) {
            edge[Index(Candidates, rankings[j][i])][j]=1;
        }
    }
}

//                                                                      //
//     Construct a graph (edge[][]) for comprehensive matching          //
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


int DUMMY::rankedOnce() {
    
    int DidTruncate;
    int i,j,ii,jj;
    int found,found2,count;
    
    // DidTruncate will serve for the return of the function.
    DidTruncate=0;
    found=1;
    
    while (found>0) {
        found2=0;
        for (j = 0 ; j < originalRankings.size() ; j++) {
            for (i = 1 ; i < originalRankings[j].size() ; i++) {
                count=0;
                // count = number of times the candidate is ranked
                for (jj = 0 ; jj < originalRankings.size() ; jj++) {
                    for (ii = 1 ; ii < originalRankings[jj].size() ; ii++) {
                        if (originalRankings[jj][ii]==originalRankings[j][i]) {
                            count=count+1;
                        }
                    }
                }
                if (count==1) {
                    //  the candidate is ranked only once
                    //  Now check is not the last ranked (otherwise there's nothing to truncate
                    if (i < originalRankings[j].size()-1) {
                        //  She/he's not the last ranked
                        //  Truncate at the candidate ranked just after him = originalRankings[j][i+1]
                        Truncate(originalRankings[j], originalRankings[j][i+1]);
                        DidTruncate=1;
                        found2=1;
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


int DUMMY::rankedMutiple(){
    
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
    
    for (j = 0 ; j < originalRankings.size() ; j++) {
        for (i = 1 ; i < originalRankings[j].size() ; i++) {
            found=0;
            for (ii = 0 ; ii < candidateSet.size() ; ii++) {
                if (candidateSet[ii]==originalRankings[j][i]) {
                    found=1;
                }
            }
            if (found==0) {
                candidateSet.push_back(originalRankings[j][i]);
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
    
    numberRanked.resize(candidateSet.size());
    for (i = 0 ; i < candidateSet.size() ; i++) {
        for (j = 0 ; j < originalRankings.size() ; j++) {
            for (ii = 1 ; ii < originalRankings[j].size() ; ii++) {
                if (originalRankings[j][ii]==candidateSet[i]) {
                    numberRanked[i].push_back(j);
                }
            }
        }
    }
    
    // Get the highest number of times a candidate is ranked
    maxSize=0;
    for (i = 0 ; i < numberRanked.size(); i++) {
        if (numberRanked[i].size()>maxSize) {
            maxSize=(int)numberRanked[i].size();
        }
    }
    
    // DidTruncate will serve for the return of the function.
    DidTruncate=0;
    
    for (int k = 2 ; k < maxSize ; k++) {
        // for each possible number of times a candidate is ranked multiple times.
        int foundAGroup;
        vector<int> theCandidates; // a set where we store the candidates ranked for the same positions
        for (i = 0 ; i < candidateSet.size()-1 ; i++) {
            // Don't need to look at all candidate: we will compare the positions ranked
            // for the i-th candidate with those of candidates with a higher index.
            theCandidates.clear();
            theCandidates.resize(0);
            foundAGroup=0;
            if (numberRanked[i].size()==k) {
                // found a candidate ranked k times
                // add him/her to the set theCandidates[]
                theCandidates.push_back(candidateSet[i]);
                for (ii = i+1 ; ii < numberRanked.size() ; ii++) {
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
                        foundAGroup=1;
                        break;
                    }
                }
            }
            if (foundAGroup==1) {
                for (int h = 0 ; h < k ; h++) {
                    // for each of the k positions, find the rank of the lowest ranked candidate in those k candidates
                    int theLowestRank=-1;
                    int thePosition=numberRanked[i][h]; // the index of the position in originalRankings
                    //cout << "The position is " << thePosition << "\n";
                    int theGuy;
                    for (ii = 0 ; ii < theCandidates.size() ; ii++) {
                        theGuy=theCandidates[ii]; // the name of the candidate -- just to simplify the code, the compiler will get rid of it.
                        // check if theGuy's rank is below the value of theLowestRank. If yes, then update the value of theLowestRank
                        //cout << "Candidate " << theGuy << " ranked " << indexInRankings(originalRankings[thePosition], theGuy)<< "\n";
                        if (indexInRankings(originalRankings[thePosition], theGuy)>theLowestRank) {
                            theLowestRank=indexInRankings(originalRankings[thePosition], theGuy);
                        }
                    }
                    if (theLowestRank<originalRankings[thePosition].size()-1) {
                        // if the theLowestRank is not the lowest rank of the position, then there's something to truncate.
                        // if not, there's nothing to truncate (and the DidTruncate is not update to take the value 1.
                        Truncate(originalRankings[thePosition], originalRankings[thePosition][theLowestRank+1]);
                        DidTruncate=1;
                    }
                }
            }
        }
    }
    return DidTruncate;
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

void DUMMY::simpleBlocks() {
    
    int NotDoneYet=1;
    int TempValue=0;
    
    //
    //  Keep running as long as rankedOnce() or RankedMultiple() performed a truncation.
    //  Stop when there's nothing left to truncate.
    //
    
    
    while (NotDoneYet>0) {
        TempValue=rankedOnce();
        TempValue=TempValue+rankedMutiple();
        if (TempValue==0) {
            NotDoneYet=0;
        }
    }
    
    /////////////////////////////////////////////
    //  Count number of positions solved
    int count,counter;
    int i,j,jj;
    counter=0;
    for (j = 0 ; j < originalRankings.size() ; j++) {
        if (originalRankings[j].size()==2) {
            // There are just 2 entries in originalRankings[j]:
            // originalRankings[j][0] is the position ID and originalRankings[j][1] is the candidate ranked 1st.
            count=0;
            // Now check that originalRankings[j][1] is ranked only once.
            for (jj = 0 ; jj < originalRankings.size() ; jj++) {
                for (i = 1 ; i < originalRankings[jj].size() ; i++) {
                    if (originalRankings[jj][i]==originalRankings[j][1]) {
                        count=count+1;
                    }
                }
            }
            if (count==1) {
                // originalRankings[j][1] only ranked once, so he/she gets that position
                // => one position solved.
                counter=counter+1;
            }
        }
    }
    cout <<  counter << " positions solved out of " << originalRankings.size() << " (" << 100*((float)counter/(float)originalRankings.size())  << "%)\n";
    //  Count number of positions solved
    /////////////////////////////////////////////
    
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


void DUMMY::loadData(std::string dataFile) {
    
    rankings.clear();
    impossibles.clear();
    Positions.clear();
    Positions.resize(0);
    OpePostes.clear();
    OpePostes.resize(0);
    Candidates.clear();
    
    originalRankings.clear();
    originalRankings.resize(0);
    
    for (int i = 0 ; i< 20 ; i++) {
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
    
    
    data_read.open(dataFile, ios::in | ios::out ); // opens the file
    if(!data_read) { // file couldn't be opened
        cerr << "Error: data file  could not be opened" << endl;
        exit(1);
    }
    size_t lines = 0;
    ifstream in(dataFile);
    for (string s; getline(in,s); ) {
        ++lines;
    }
    OpePostes.resize(lines);
    for (int i = 0 ; i < OpePostes.size() ; i++) {
        OpePostes[i].resize(0);
    }
    
    
    int row = 0; // Row counter
    cout << "\n";
    while (!data_read.eof()) {
        int data;
        for (int i = 0 ; i < 4 ; i++) {
            data_read >> data;
            cout << data << "\t";
            OpePostes[row].push_back(data);
        }
        cout << "\n";
        row++;
    }
    data_read.clear();
    data_read.close();
    cout << "\n";
    
    
    for (row = 0 ; row < OpePostes.size(); row++) {
        //  for each position id (in OpePostes[][3]
        //  check if already in Position[]. If not add it.
        if (!presentInVector(Positions, OpePostes[row][3])) {
            Positions.push_back(OpePostes[row][3]);
        }
    }
    
    //
    //  Construct now the originalRankings[][]
    //
    
    originalRankings.resize(Positions.size());
    for (int h = 0 ; h < Positions.size() ; h++) {
        // For each position, put in original_ranking[][0] the id of the position
        originalRankings[h].push_back(Positions[h]);
        for (int i = 0 ; i < 20 ; i++) {
            TempRanking[i]=-1;
        }
        // Look for all candidates ranked by Position[row]
        for (int j = 0 ; j < OpePostes.size() ; j++) {
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
        for (int i = 1 ; i < 20 ; i++) {
            if (TempRanking[i]!=-1) {
                originalRankings[h].push_back(TempRanking[i]);
            }
        }
    }
}


//                                                                      //
//     Load and prepare the data                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Analyze the data                                                 //
//                                                                      //
//                                                                      //

//
//  There's a little bit more data preparation, then searching
//  for the impossible matches
//


void DUMMY::analyzeData() {
    
    
    // first simplify the ranking profile by eliminating
    // the impossible matches because of simple blocks.
    simpleBlocks();
    
    //
    //  Declare impossibles (-1 = don't know, 0 = not impossible, 1 = impossible
    //
    impossibles.resize(originalRankings.size());
    for (int i = 0 ; i < originalRankings.size(); i++) {
        impossibles[i].push_back(originalRankings[i][0]);
        impossibles[i].push_back(0);
        for (int j = 2 ; j < originalRankings[i].size(); j++) {
            impossibles[i].push_back(-1);
        }
    }
    
    //
    //  Construct set of candidates
    //
    Candidates.clear();
    Candidates.resize(0);
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        for (int i = 1 ; i < originalRankings[j].size() ; i++) {
            if (!presentInVector(Candidates, originalRankings[j][i])) {
                Candidates.push_back(originalRankings[j][i]);
            }
        }
    }
    
    
    int TotalCases=0;
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        for (int i = 2 ; i < originalRankings[j].size() ; i++) {
            TotalCases=TotalCases+1;
        }
    }
    
    //
    //	Find the impossible matches.
    //
    int CaseNumber=1;
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        for (int i = 2 ; i < originalRankings[j].size() ; i++) {
            if (impossibles[j][i]==-1) {
                checkImpossible(originalRankings[j][i], j);
            }
            CaseNumber=CaseNumber+1;
        }
    }
    
    // Run a post analysis
    postAnalysis();
    
    
}
//                                                                      //
//     Analyze the data                                                 //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Computing how many simulations of preferences we need            //
//                                                                      //
//                                                                      //



void DUMMY::postAnalysis() {
    
    int TopCandidates;
    int NumbCombinations;
    NumbCombinations=1;
    TopCandidates=0;
    int NumberTop,NumberGeneral;
    vector<int> MarketStars;
    MarketStars.clear();
    MarketStars.resize(0);
    
    
    int NumberCases=0;
    int NumberImpossibles=0;
    
    for (int j = 0 ; j < impossibles.size() ; j++) {
        for (int i = 1 ; i < impossibles[j].size() ; i++) {
            NumberCases=NumberCases+1;
            if (impossibles[j][i]==1) {
                NumberImpossibles=NumberImpossibles+1;
            }
        }
    }
    
    cout << "There are " << NumberImpossibles << " impossible pairs out of " << NumberCases << "\n";
    
    
    cout << "\n\nSimulations of (partial) preferences (next step)\n";
    for (int j = 0 ; j < originalRankings.size() ; j++) {
        NumberTop = 0;
        NumberGeneral = 0;
        for (int jj = 0 ; jj < originalRankings.size() ; jj++) {
            if (originalRankings[j][1]==originalRankings[jj][1]) {
                NumberTop=NumberTop+1;
            }
            for (int i = 1 ; i < originalRankings[jj].size() ; i++) {
                if (originalRankings[j][1]==originalRankings[jj][i]) {
                    NumberGeneral=NumberGeneral+1;
                }
            }
        }
        if (NumberTop==NumberGeneral) {
            if (NumberTop>1) {
                bool CandIsPresent = (std::find(MarketStars.begin(), MarketStars.end(), originalRankings[j][1]) != MarketStars.end());
                if (!CandIsPresent) {
                    MarketStars.push_back(originalRankings[j][1]);
                    TopCandidates=TopCandidates+1;
                    NumbCombinations=NumbCombinations*NumberTop;
                    cout << "Candidate " << originalRankings[j][1] << " ranked " << NumberTop << " times 1st\n";
                }
            }
        }
    }
    cout << "\n\nThere are " << TopCandidates << " top candidates (ranked only first & 2 or more times)\n";
    cout << "Need to consider " << NumbCombinations << " possible pref profiles\n";
}

//                                                                      //
//     Computing how many simulations of preferences we need            //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Find a maximum & comprehensive matching                          //
//                                                                      //
//                                                                      //

void DUMMY::comprehensiveMatching(){
    
    int keepGoing = 1;
    int previousNumber = 0;
    int count = 0;

    // Initialize the variables
    matching_candidate.clear();
    matching_candidate.resize(Candidates.size());
    matching_department.clear();
    matching_department.resize(rankings.size());
    
    for (int i = 0 ; i < Candidates.size() ; i++) {
        matching_candidate[i]=-1;
    }
    for (int j = 0 ; j < rankings.size(); j++) {
        matching_department[j]=-1;
    }
    
    //  Loop to construct a maximum matching.
    while (keepGoing>0) {
        admissibleGraph();
        bipMatch();
        count = 0;
        for (int i = 0 ; i < Candidates.size() ; i++) {
            if (matching_candidate[i]!=-1) {
                count=count+1;
            }
        }
        if (count==previousNumber) {
            keepGoing=0;
        } else {
            previousNumber=count;
        }
    }
}

//                                                                      //
//     Find a maximum & comprehensive matching                          //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Given a pair (candidate,position) simplify the rankings          //
//                                                                      //
//                                                                      //


//  Want to check if candidate i0 is an impossible match for position s0
//  We construct from originalRankings[][] a profile rankings[][] specfic for that problem:
//  - the first position int rankings[][] is s0
//  - candidate i0 only appears in s0's ranking.


void DUMMY::Simplify_i0(int i0, int s0) {
    
    Candidates.clear();
    Candidates.resize(0);
    Positions.clear();
    Positions.resize(0);
    rankings.clear();
    rankings.resize(Positions.size());

    
    for (int i = 1 ; i < (indexInRankings(originalRankings[s0], i0)+1) ; i++) {
        // index(i0)+1 to include i0, and nobody after.
        Candidates.push_back(originalRankings[s0][i]);
    }
    Positions.push_back(s0);
    // so now Positions[0] = s0;
    
    // Candidates[] first contains the names of all candidates better ranked than i0 at s0
    // Positions[0] = s0 (more precisely, the index of position s0 in originalRankings[][]).
    //
    // What we do below:
    // for each candidate i in Candidates[], find all candidates better ranked than i at any position where i is ranked,
    // then add those candidates to Candidates[]
    // then add the position to Positions[].
    // Positions[j] = the j-th entry in Positions, the index of the position in originalRankings.
    // so, if the 3rd entry in Positions[] is the 17th position in originalRankings[][], whose name is 84 we have
    // Positions[3]=17, and originalRankings[17][0]=84.
    
    for (int i = 0 ; i < Candidates.size() ; i++) {
        if (Candidates[i]!=i0) {
            for (int j = 0 ; j < originalRankings.size() ; j++) {
                if (indexInRankings(originalRankings[j], Candidates[i])>0) {
                    for (int ii = 1 ; ii < indexInRankings(originalRankings[j], Candidates[i]) ; ii++) {
                        if (!presentInVector(Candidates, originalRankings[j][ii])) {
                            Candidates.push_back(originalRankings[j][ii]);
                        }
                    }
                    if (!presentInVector(Positions, j)) {
                        Positions.push_back(j);
                    }
                }
            }
        }
    }
    
    
    // Construct rankings[][], a copy of originalRankings[][]
    // such that i0 is deleted everywhere except at s0.
    // and such that rankings[][] are only about the candidates in Candidates[] and the positions in Positions[]
    
    for (int j = 0 ; j < Positions.size() ; j++) {
        rankings[j].resize(0);
        rankings[j].push_back(originalRankings[Positions[j]][0]);
        //
        //	The 3rd entry in Positions/rankings is the 17th position in originalRankings[], whose name is 84:
        // rankings[3][0] = 17, originalRankings[17][0] = 84
        //
        //	Next we fill rankings[j][1], rankings[j][2], etc, but only with candidates in Candidates[].
        //	If j != 0 (it's not j0, then we skip i0
        //	If j=0 we don't skip i0.
        for (int i = 1 ; i < originalRankings[Positions[j]].size(); i++) {
            if (originalRankings[Positions[j]][i]!=i0) {
                bool CandIsPresent = (std::find(Candidates.begin(), Candidates.end(), originalRankings[Positions[j]][i]) != Candidates.end());
                if (CandIsPresent) {
                    rankings[j].push_back(originalRankings[Positions[j]][i]);
                }
            } else {
                if (j==0) {
                    bool CandIsPresent = (std::find(Candidates.begin(), Candidates.end(), originalRankings[Positions[j]][i]) != Candidates.end());
                    if (CandIsPresent) {
                        rankings[j].push_back(originalRankings[Positions[j]][i]);
                    }
                }
            }
        }
    }
}

//                                                                      //
//     Given a pair (candidate,position) simplify the rankings          //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Check if a matching is comprehensive                             //
//                                                                      //
//                                                                      //

void DUMMY::checkComprehensiveness(){
    
    for (int j = 0 ; j < rankings.size(); j++) {
        for (int i = 2 ; i < rankings[j].size() ; i++) {
            if (matching_candidate[Index(Candidates,rankings[j][i])]==j) {
                for (int ii = 1 ; ii < i ; ii++) {
                    if (matching_candidate[Index(Candidates,rankings[j][ii])]==-1) {
                        cout << "\n\n\n ERROR NOT COMPREHENSIVE \n\n\n\n";
                        exit (1);
                    }
                }
                
            }
        }
    }
}

//                                                                      //
//     Check if a matching is comprehensive                             //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Check if i0 is an impossible match for s0                        //
//                                                                      //
//                                                                      //


void DUMMY::checkImpossible(int i0, int s0) {
    
    // i0 comes from originalRankings[j][i] = name of a candidate
    // s0 comes from j = index of the deparment in originalRankings[][] (not the name!)
    
    int decided = 0;
    
    Simplify_i0(i0, s0);
    comprehensiveMatching();
    checkComprehensiveness();
    
    //  Gamma_dep and Gamma_cand are sets of positions and candidates that will
    //  be use to check whether i0 is an impossible match
    //  Gamma_cand is the set $\mathbf{J}$ in the characterization of a block (see the paper)
    //  Gamma_dep is the set of acceptable positions for candidates in Gamma_cand
    //
    
    Gamma_dep.clear();
    Gamma_dep.resize(0);
    Gamma_cand.clear();
    Gamma_cand.resize(0);
    J_0.clear();
    J_0.resize(0);
    
    if (matching_candidate[Index(Candidates, i0)]==0) {
        // so i0 is matched to Positions[0]
        // but Positions[0] is s0.
        // so i0 is not an impossible match for s0.
        decided=1;
        impossibles[s0][Index(originalRankings[s0], i0)]=0;
        //cout << "i0 matched, not impossible";
        // recall that s0 is the s0-th department in originalRankings.
        
        
        // Since i0 is matched at a comprehensive in rankings[][]
        // the matching is still comprenhensive in originalRankings[][]
        // so all the students matched to a position are not impossible for that position
        
        for (int i = 0 ; i < Candidates.size() ; i++) {
            if (matching_candidate[i]!=-1) {
                // The candidate is matched to a position
                int theCandidate=Candidates[i];
                int thePosition = matching_candidate[i]; // the candidate is matched to the matching_candidate[i]-th position in rankings[][]
                int indexPosition;
                int indexCandidate;
                // Now find the index of that position in originalRankings[][]
                for (int j = 0 ; j < originalRankings.size() ; j++) {
                    if (originalRankings[j][0]==rankings[thePosition][0]) {
                        indexPosition=j;
                        break;
                    }
                }
                indexCandidate=Index(originalRankings[indexPosition], theCandidate);
                // The candidate is ranked indexCandidate-th at the position in originalRankings
                // It may not be the sama as in rankings because in rankings we deleted i0 from the rankings (except at s0).
                if (originalRankings[indexPosition][indexCandidate]!=theCandidate) {
                    cout << "\n\n\nError!!!!!!\n";
                    exit(1);
                }
                impossibles[indexPosition][indexCandidate]=0;
            }
        }
    }
    if (decided==0) {
        // i0 is not matched to s0 => don't know yet if impossible.
        
        //////////////////
        //  build J_0
        //
        for (int i = 1 ; i < indexInRankings(rankings[0], i0); i++) {
            J_0.push_back(rankings[0][i]);
        }
        if (J_0.size()==0) {
            // shouldn't happen. If J0 is empty then i0 should be matched to s0 at the comprehensive matching.
            cout << "\n\nERROR --- J0 empty\n\n";
            exit(1);
        }
        
        //////////////////
        //  build Gamma
        //
        buildGamma(J_0,i0);
    }
    while (decided==0) {
        // Eliminate s0 from Gamma_dep: we want to reserve it for i0.
        // Only candidates in Gamma_cand will be matched
        // and i0 is NOT in Gamma_cand
        
        Gamma_dep.erase(Gamma_dep.begin()+0);
        
        if (Gamma_dep.size()>Gamma_cand.size()) {
            cout << "Error: Gamma dep > Gamma_cand ";
            exit(1);
        }
        
        // we'll need to truncate at some set K
        // the vector possible_K[] contains all the candidates
        // that we can put in K (candidates ranked 1st for a position are not eligible for K, they are prevalent --- see paper).
        //
        
        buildPossible_K();
        sizeK = (Gamma_cand.size()-Gamma_dep.size());
        // There's no the "+1" in the formula of sizeK like in the paper because s0 eliminated from Gamma_dep[].
        // (In the paper, s0 is part of Gamma_dep[])
        
        
        if (possible_K.size()>0) {
            //cout <<  "combinations: C(" << possible_K.size() << "," << sizeK << ") = ";
            //cout << coeff((int) possible_K.size(), (int) sizeK) << "\t";
        } else {
            // K is necessarily empty, so i0 is an impossible match.
            impossibles[s0][Index(originalRankings[s0], i0)]=1;
            //cout << "impossible\n";
            break;
        }
        if (possible_K.size()<sizeK) {
            // In this case the 2nd condition of the characterization of a block has no bite
            // The set Gamma_cand satisfies the 1st condition, so i0 is an impossible match for s0
            impossibles[s0][Index(originalRankings[s0], i0)]=1;
            //cout << "impossible\n";
            break;
        } else {
            
            LoopNumber=1;
            int stop=0;
            //cout << "going for the loops ";
            while (stop==0) {
                // we stop when we find either a maximum matching at the truncation for some K
                // or we have considered all possible sets K.
                K.clear();
                K.resize(0);
                // construct the set K
                for (int j = 0 ; j < sizeK ; j++) {
                    K.push_back(possible_K[selectedIndex[j]]);
                }
                // Construct the graph such that only admissible candidates in rankings[][]
                // truncated at K have and edge.
                graphOfJ();
                
                
                // reinitialize the matching
                
                matching_candidate.clear();
                matching_candidate.resize(Candidates.size());
                matching_department.clear();
                matching_department.resize(rankings.size());
                
                for (int i = 0 ; i < Candidates.size() ; i++) {
                    matching_candidate[i]=-1;
                }
                for (int j = 0 ; j < rankings.size() ; j++) {
                    matching_department[j]=-1;
                }
                
                // find a maximum matching and count how many students are eligible
                // The function graphOfJ() reconstructed the vector Candidates[], it now contains only
                // admissible candidates
                
                bipMatch();
                int NumberOfMatchedCandidates=0;
                for (int i = 0 ; i < Candidates.size(); i++) {
                    if (matching_candidate[i]!=-1) {
                        NumberOfMatchedCandidates++;
                    }
                }
                int NumberCandidateToBeMatched=0;
                for (int i = 0 ; i < Candidates.size() ; i++) {
                    for (int j = 0 ; j < rankings.size() ; j++) {
                        if (edge[i][j]==1) {
                            NumberCandidateToBeMatched=NumberCandidateToBeMatched+1;
                            break;
                        }
                    }
                }
                
                if (NumberOfMatchedCandidates==Candidates.size()) {
                    // all admissible candidates are matched. So condition 2 of the definition of a block is violated
                    // i0 is not an impossible match
                    impossibles[s0][Index(originalRankings[s0], i0)]=0;
                    //cout << "not impossible";
                    stop=1; // don't need to do another loop (with a different K)
                    decided=1; // The whole thing stop
                    break;
                } else {
                    if ((int) coeff((int) possible_K.size(), (int) sizeK)==(int)LoopNumber) {
                        // we looked at all possible cases: condition 2 is not violated
                        // i0 is an impossible match for s0.
                        impossibles[s0][indexInRankings(originalRankings[s0], i0)]=1;
                        //cout << "impossible";
                        // Gamma_cand is a block. Check if the candidates below i0 are in Gamma_cand. If not, Gamma_cand is also a block for them.
                        for (int ii = indexInRankings(originalRankings[s0], i0)+1 ; ii < originalRankings[s0].size() ; ii++) {
                            if (!presentInVector(Gamma_cand, originalRankings[s0][ii])) {
                                impossibles[s0][ii]=1;
                            }
                        }
                        stop=1;
                        decided=1;
                        break;
                    } else {
                        // need to go for another subset K taken from possible_K
                        next_subset(selectedIndex, (int) possible_K.size(), (int) sizeK);
                        LoopNumber++;
                        //cout << "loop number = " << LoopNumber << "/" <<  coeff((int) possible_K.size(), (int) sizeK) << "\n";
                    }
                }
            }
        }
    }
}

//                                                                      //
//     Check if i0 is an impossible match for s0                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     build the set of candidates eligible for the set K               //
//                                                                      //
//                                                                      //

void DUMMY::buildPossible_K() {
    possible_K.clear();
    possible_K.resize(0);
    //  First make possible_K a copy of Gamma_cand
    for (int i = 0 ; i < Gamma_cand.size(); i++) {
        possible_K.push_back(Gamma_cand[i]);
    }
    // Candidates in J_0 can never be in K.
    for (int i = 0 ; i < J_0.size() ; i++) {
        possible_K.erase(possible_K.begin()+(Index(possible_K, J_0[i])));
    }
    // Candidates ranked first for a position are prevalent (see the paper)
    // don't need them in K.
    for (int j = 0 ; j < rankings.size() ; j++) {
        if (presentInVector(possible_K, rankings[j][1])) {
            possible_K.erase(possible_K.begin()+(Index(possible_K, rankings[j][1])));
        }
    }
    selectedIndex.clear();
    selectedIndex.resize(0);
    
    for (int i = 0 ; i < possible_K.size() ; i++) {
        selectedIndex.push_back(i);
    }
}

//                                                                      //
//     build the set of candidates eligible for the set K               //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Constructing Gamma_cand (and Gamma_dep)                          //
//                                                                      //
//                                                                      //


//
//  One of the most critical part of the code. The bigger Gamma_cand
//  The more likely we'll have a LARGE number of combinations to look at
//  Gamma_cand is initially built like in the paper, but some hacks are added
//  to make is (possibly) smaller.
//

void DUMMY::buildGamma(vector<int> &J_0, int i0) {
    
    // Initialize the beast...
    Gamma_dep.clear();
    Gamma_dep.resize(0);
    Gamma_cand.clear();
    Gamma_cand.resize(0);
    
    ///////////////////////////////////////////////
    //
    //  Construct Gamma_cand as in the paper
    //
    //  (the sets J^h and S^h, h=1,...)
    //
    
    // Gamma_cand first contains all candidates in J_0, i.e., candidates ranked better than i0 at s0.
    for (int i=0 ; i < J_0.size(); i++) {
        Gamma_cand.push_back(J_0[i]);
    }
    // Then add all the candidates who are matched to some department.
    for (int i = 0 ; i < Candidates.size() ; i++) {
        if (matching_candidate[i]!=-1) {
            if (!presentInVector(Gamma_cand, Candidates[i])) {
                Gamma_cand.push_back(Candidates[i]);
            }
        }
    }
    
    // construction of out_dep[], a set of positions to withdraw.
    // out_dep[] is a set of positions' indices.
    // out_dep[j]=h => the j-th position in out_dep is the h-th position in rankings[],
    // whose name is rankings[h][0]=rankings[out_dep[j]][0]
    
    vector<int> out_dep;
    out_dep.clear();
    out_dep.resize(0);
    vector<int> out_cand;
    out_cand.clear();
    out_cand.resize(0);
    
    for (int j = 0 ; j < rankings.size(); j++) {
        if (matching_department[j]==-1) {
            out_dep.push_back(j);
        }
    }
    
    int found =1;
    int found2=0;
    
    while (found>0) {
        // Compute J^h
        // Add to out_cand[] the candidates in Gamma_cand[] who are acceptable for a department in out_dep[]
        for (int i = 0 ; i < Gamma_cand.size() ; i++) {
            for (int j = 0 ; j < out_dep.size() ; j++) {
                if (Index(rankings[out_dep[j]], Gamma_cand[i])>-1) {
                    if (!presentInVector(out_cand, Gamma_cand[i])) {
                        out_cand.push_back(Gamma_cand[i]);
                        //cout << "Eliminating candidate " << Gamma_cand[i] << " because acceptable for " << rankings[out_dep[j]][0] << "\n";
                    }
                }
            }
        }
        // eliminate from Gamma_cand[] all the candidates in out_cand[]
        for (int i = 0 ; i < out_cand.size() ; i++) {
            if (presentInVector(Gamma_cand, out_cand[i])) {
                Gamma_cand.erase(Gamma_cand.begin()+Index(Gamma_cand, out_cand[i]));
            }
            if (presentInVector(J_0, out_cand[i])) {
                J_0.erase(J_0.begin()+(Index(J_0, out_cand[i])));
            }
        }
        
        //////// Compute S^h
        //
        found2=0;
        // search a position that is not yet in out_dep
        for (int j = 0 ; j < rankings.size(); j++) {
            if (!presentInVector(out_dep, j)) {
                // the j-th position in rankings[][] is not in out_dep.
                if (matching_department[j]!=-1) {
                    // that position is matched to the matching_department[j]-th candidate in Candidates.
                    // The name of that candidate is Candidates[matching_department[j]]
                    if (!presentInVector(Gamma_cand, Candidates[matching_department[j]])) {
                        // that candidate is not in Gamma_cand, so we add the position to out_dep[]
                        out_dep.push_back(j);
                        //cout << "adding dep " << rankings[j][0] << " to out_dep (matched to a candidate not in Gamma_cand)\n";
                        found2=1;
                    }
                }
            }
        }
        if (Gamma_cand.size()==0) {
            cout << "Gamma cand size = 0\n";
            exit(1);
        }
        if (found2==0) {
            found=0;
        }
    }
    //
    //	Now that we have the set out_dep, reconstruct Gamma_dep that contains all the department
    //	that are not in out_dep
    //
    Gamma_dep.clear();
    Gamma_dep.resize(0);
    for (int j = 0 ; j < rankings.size(); j++) {
        bool DepIsPresent = (std::find(out_dep.begin(), out_dep.end(), j) != out_dep.end());
        if (!DepIsPresent) {
            Gamma_dep.push_back(j);
        }
    }
    // Check that s0 (= position 0) is not in out_dep
    bool DepIsPresent = (std::find(out_dep.begin(), out_dep.end(), 0) != out_dep.end());
    if (DepIsPresent) {
        cout << "\n\nError\n s0 is in out_dep[]\n\n";
        exit(1);
    }
    // Check that Gamma_dep[0] is s0
    if (Gamma_dep[0]!=0) {
        cout << "\n\nError\n s0 is not Gamma_dep[0]\n\n";
        exit(1);
    }
    //
    //  Construct Gamma_cand as in the paper
    //
    //  (the sets J^h and S^h, h=1,...)
    //
    ///////////////////////////////////////////////
    
    
    //
    //		Gamma_cand only contains matched candidates. We need to add unmatched candidates. 
    //		The while loop below does that: Should only add unmatched candidates that are ranked better than a candidate in Gamma.
    //		But as we add new candidates in Gamma, new potential candidates for Gamma may show up. Loop stops
    //		when we don't add new candidates.
    
    found =1;
    while (found>0) {
        found2=0;
        for (int j = 1; j < Gamma_dep.size() ; j++) {
            for (int i = 1 ; i < rankings[Gamma_dep[j]].size() ; i++) {
                // pick the candidate ranked i-th
                int IndexOfCandidate = Index(Candidates, rankings[Gamma_dep[j]][i]);
                // IndexOfCandidate is the index in Candidates[] of the candidate
                if (matching_candidate[IndexOfCandidate]==-1) {
                    // The candidate is not matched
                    if (rankings[Gamma_dep[j]][i]!=i0) {
                        // if it's no i0, see whether he/she's in Gamma_cand.
                        if (!presentInVector(Gamma_cand, rankings[Gamma_dep[j]][i])) {
                            Gamma_cand.push_back(rankings[Gamma_dep[j]][i]);
                            //cout << "adding candidate " << rankings[Gamma_dep[j]][i] << " not matched but acceptable\n";
                            found2 = 1;
                        }
                    }
                }
            }
        }
        if (found2==0) {
            found=0;
        }
    }
    
    // Check that all candidates in J0 are also acceptable for a position different from s0
    // if there is a candidate for which it is not the case (i.e., ranked only by s0)
    // then candidates below (including i0) are impossible for s0
    // and thus in the reduced profile (originalRankings[s0][]) i0 is not ranked by s0.
    for (int i = 0 ; i < J_0.size() ; i++ ) {
        found = 0;
        for (int j = 0 ; j < Gamma_dep.size() ; j++) {
            if (Index(rankings[Gamma_dep[j]], J_0[i])>-1) {
                found=1;
            }
        }
        if (found==0) {
            cout << "\n\nError --- candidate in J_0 only acceptable for s0\n\n";
            exit(1);
        }
    }
    
    
    
    //  Up to know the set Gamma_cand is the same as the set J in the paper
    //
    //	Now construct the connected component of Gamma_cand that includes J0
    //	(i.e., there exists a path from J_0 to any other candidate in Gamma
    //
    //	Then reconstruct Gamma_dep)
    //
    
    // First copy Gamma_cand to temp_Gamma_cand, and then empty Gamma_cand (that we will reconstruct).
    vector<int>  temp_Gamma_cand;
    temp_Gamma_cand.clear();
    temp_Gamma_cand.resize(0);
    for (int i = 0 ; i < Gamma_cand.size() ; i++) {
        temp_Gamma_cand.push_back(Gamma_cand[i]);
    }
    Gamma_cand.clear();
    Gamma_cand.resize(0);
    
    
    // Do the same with Gamma_dep (copy it to temp_Gamma_dep, etc.).
    vector<int>  temp_Gamma_dep;
    temp_Gamma_dep.clear();
    temp_Gamma_dep.resize(0);
    for (int j = 0 ; j < Gamma_dep.size() ; j++) {
        temp_Gamma_dep.push_back(Gamma_dep[j]);
    }
    Gamma_dep.clear();
    Gamma_dep.resize(0);
    
    // add s0 to Gamma_dep.
    
    
    
    if (presentInVector(Gamma_cand, i0)) {
        cout << "i0 in Gamma 1\n";
        exit(1);
    }
    
    //
    // We start with s0 (want to get the connected component containing s0)
    //
    Gamma_dep.push_back(0);
    
    found=1;
    while (found>0) {
        for (int j = 0 ; j < Gamma_dep.size(); j++) {
            // take a position in Gamma_dep
            found2=0;
            for (int i = 1 ; i < rankings[Gamma_dep[j]].size(); i++) {
                if (rankings[Gamma_dep[j]][i]!=i0) {
                    // Take a candidate ranked by Gamma_dep
                    //cout << "Candidat " << rankings[Gamma_dep[j]][i] << "\n";
                    // Find candidates that are in temp_Gamma_cand that are acceptable for a department in Gamma_dep
                    if (presentInVector(temp_Gamma_cand, rankings[Gamma_dep[j]][i])) {
                        // It's a candidate we want in Gamme_dep
                        // Now check if that candidate is already in Gamma_cand. If not add it to Gamma_cand.
                        if (!presentInVector(Gamma_cand, rankings[Gamma_dep[j]][i])) {
                            //cout << "\tadded to Gamma_cand (" << rankings[Gamma_dep[j]][i] << ")\n";
                            Gamma_cand.push_back(rankings[Gamma_dep[j]][i]);
                            found2=1;
                            // Now add all departments where candidate rankings[Gamma_dep[j]][i] is acceptable
                            int TheCurrentCandidate = rankings[Gamma_dep[j]][i];
                            for (int jj = 0 ; jj < rankings.size() ; jj++) {
                                // Find positions that rank the candidate
                                if (Index(rankings[jj], TheCurrentCandidate)>-1) {
                                    // Check if the position where he's ranked is in Gamma_dep. If not, add the position to Gamma_dep.
                                    if (!presentInVector(Gamma_dep, jj)) {
                                        Gamma_dep.push_back(jj);
                                        //cout << "Adding department " << rankings[jj][0] << "\n";
                                        found2 =1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        if (found2==0) {
            found=0;
        }
    }
    
    //
    //  Gamma_dep contains s0. But we withdraw if from Gamma_dep:
    //  Later we'll try to match candidates in Gamma_cand to a position in Gamma_dep
    //  but without matching any to s0 (that position is reserved for i0).
    //
    
    Gamma_dep.erase(Gamma_dep.begin()+0);
    
    // Check J_0 is not empty
    if (J_0.size()==0) {
        cout << "Error: J_0 empty after computing Gamma\n";
        exit(1);
    }
    
    if (Gamma_dep.size() > Gamma_cand.size()) {
        cout << "error size, Gamma_dep > Gamma_cand\n";
        exit(1);
    }
    if (presentInVector(Gamma_cand, i0)) {
        cout << "i0 in Gamma\n";
        exit(1);
    }
}


//                                                                      //
//     Constructing Gamma_cand (and Gamma_dep)                          //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Computing how many subsets of size k from a set of size n        //
//                                                                      //
//                                                                      //

signed long int DUMMY::coeff(int n, int k) {
    int i,j;
    signed  long int prev[n];
    signed  long int current[n];
    
    for (i = 0 ; i < n ; i++) {
        prev[i]=1;
        current[i]=1;
    }
    for (i = 1 ; i < n ; i++) {
        for (j = 1 ; j < i+1; j++) {
            current[j]=prev[j-1]+prev[j];
        }
        for (j = 0 ; j < n ; j++) {
            prev[j]=current[j];
        }
    }
    if (n==k) {
        return (1);
    }
    return current[k];
}

//                                                                      //
//     Computing how many subsets of size k from a set of size n        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Adjacency matrix for Gamma_can truncated at some set K           //
//                                                                      //
//                                                                      //

//
//  Given Gamma_cand and Gamma_dep and a set K\subseteq Gamma_can\J_0
//  we truncate the rankings at K and re-construct the edge[][].
//  The graph that we obtain will be used to find a maximum matching.
//


void DUMMY::graphOfJ(){
    
    // Given rankings[][] and K[]
    // construct a graph such that it contains only students in Gamma_cand
    // when rankings[][] is truncated at K[];
    
    // Declare temp_rankings
    // We don't want to truncate rankings[][], they might serve for the next subset K.
    // So we just construct a temporary ranking.
    
    
    
    edge.clear();
    edge.resize(Candidates.size());
    for (int i = 0 ; i < Candidates.size() ; i++) {
        for (int j = 0 ; j < rankings.size(); j++) {
            edge[i].push_back(0);
        }
    }
    
    
    for (int j = 0 ; j < rankings.size() ; j++) {
        if (presentInVector(Gamma_dep, j)) {
            for (int i = 1 ; i < rankings[j].size() ; i++) {
                if (presentInVector(K, rankings[j][i])) {
                    break;
                } else {
                    if (presentInVector(Gamma_cand, rankings[j][i])) {
                        int theCandidate;
                        theCandidate=Index(Candidates, rankings[j][i]);
                        edge[theCandidate][j]=1;
                    }
                }
                
            }
            
        }
    }
}


//                                                                      //
//     Adjacency matrix for Gamma_can truncated at some set K           //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Enumerating all subsets of size k                                //
//                                                                      //
//                                                                      //

void DUMMY::next_subset(vector<int> &selectedIndex, int n, int k) {
    int i,j;
    int movable[k];
    for (i = 0 ; i < k ; i++) {
        movable[k]=0;
    }
    for (i = 0 ; i < k ; i++) {
        if (selectedIndex[i]==n-k+i) {
            movable[i]=0;
        } else {
            movable[i]=1;
        }
    }
    if (movable[k-1]==1) {
        selectedIndex[k-1]=selectedIndex[k-1]+1;
    } else {
        for (i = 0 ; i < k ; i++) {
            if (movable[k-i-1]==1) {
                selectedIndex[k-i-1]=selectedIndex[k-i-1]+1;
                for (j = k-i; j < k ; j++) {
                    selectedIndex[j]=selectedIndex[j-1]+1;
                }
                i=k;
            }
        }
    }
}

//                                                                      //
//     Enumerating all subsets of size k                                //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////






