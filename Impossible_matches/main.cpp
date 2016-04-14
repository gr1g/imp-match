//
//  main.cpp
//  Impossible_matches
//
//  Created by Guillaume Haeringer on 4/28/14.
//  Copyright (c) 2014 Guillaume Haeringer. All rights reserved.
// Hello Guillaume

#include <iostream>
#include "dummy.h"
#include <iostream>
#include <iomanip>
#include <time.h>

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using namespace std;

int main() {
    int year_index;
    string filename;
    time_t t1,t2,startSimulations,endSimulations;    // time variables to count the time spent running the program
    time_t bigStart,bigEnd;
    filename="none";
    ofstream output;
    vector<vector<int> > dataRankings;
    year_index=0;
    string currentYear;
    
    bool doingReshuffle = true;
    bool doingRandom = false;
    
    
    // What     year_index
    // Random   -1
    // 1999     0
    // 2000		1
    // 2001		2
    // 2002		3
    // 2003		4
    // 2004		5
    // 2005		6
    // 2006		7
    // 2007		8
    // 2008		9
    // 2009		10
    // 2010		11
    // 2011		12
    // 2012		13
    // 2013		14
    DUMMY dummy;
    
    dummy.finalRestults.open("final-results.txt", std::ofstream::out | std::ofstream::app);
    dummy.finalRestults << "Year\tNumber Positions\tNumber Candidate\tCand Predicted\tPositions predicted\tCand Predicted (final)\tPositions predicted (final)";
    dummy.finalRestults.close();
    
    dummy.finalRestults.open("marketStars.txt", std::ofstream::out | std::ofstream::app);
    dummy.finalRestults << "year\tPositions extracted\tNumber Stars\tNumber Pseudo Stars\n";
    dummy.finalRestults.close();

    
    dummy.finalRestults.open("acceptability_violation.txt", std::ofstream::out | std::ofstream::app);
    dummy.finalRestults << "Year\tNumber Violations\tNumber Violators\n";
    dummy.finalRestults.close();
    
    dummy.finalRestults.open("Statistics-hired.txt", std::ofstream::out | std::ofstream::app);
    dummy.finalRestults << "year\t";
    for (int i = 1 ; i < 15 ; i++) {
        dummy.finalRestults << i << "\t";
    }
    dummy.finalRestults << "\n";
    dummy.finalRestults.close();
    
    
    
    (void) time(&bigStart);
    
    dummy.finalRestults.open("ranking-distributions.txt", std::ofstream::out | std::ofstream::app);
    dummy.finalRestults << "year\t";
    for (int i = 1 ; i < 20 ; i++) {
        dummy.finalRestults << i << "\t";
    }
    dummy.finalRestults << "\n";
    dummy.finalRestults.close();
    
    
    for (year_index = 0 ; year_index < 15 ; year_index++) {
        
        string currentYear = to_string(year_index+1999);
        filename = currentYear + "-data.txt";
        dummy.fileNameResults = "Results-" +(currentYear) + ".txt";
        dummy.fileRankingsStepOne = "Rankings-StepOne-" +(currentYear) + ".txt";
        output.open(dummy.fileNameResults.c_str(), ofstream::out | ofstream::app);
        cout << "\n\n";
        cout << "*************************************\n*\n*";
        output << "\n\n";
        output << "*************************************\n*\n*";
        if (year_index==-1) {
            cout << "\tRandom\n*\n*";
            output << "\tRandom\n*\n*";
        } else {
            cout << "\tYear " << year_index+1999 << "\n*\n*";
            output << "\tYear " << year_index+1999 << "\n*\n*";
        }
        cout << "*************************************\n\n";
        output << "*************************************\n\n";
        output.close();
        (void) time(&t1);
        dummy.year = year_index+1999;
        if (year_index==-1) {
            for (int y = 0 ; y < 15 ; y++) {
                (void) time(&startSimulations);
                currentYear = to_string(y+1999);
                filename = currentYear + "-data.txt";
                dataRankings.clear();
                dataRankings.resize(0);
                dummy.loadData(filename.c_str());
                dataRankings=dummy.originalRankings;
                size_t sizeRankings = dummy.originalRankings.size();
                if (doingRandom) {
                    vector<int> virtualCandidatesQueue;
                    virtualCandidatesQueue.clear();
                    virtualCandidatesQueue.resize(0);
                    for (int j = 0 ; j < sizeRankings ; j++) {
                        for (int i = 1 ; i < dummy.originalRankings[j].size() ; i++) {
                            virtualCandidatesQueue.push_back(dummy.originalRankings[j][i]);
                        }
                    }
                    dummy.fileNameResults = "Results-random_profiles-" + currentYear + ".txt";
                    cout << "\n\n\n*************************************\n*";
                    cout << "\t\t" << y+1999 << "\tRandom Profiles\n";
                    cout << "*************************************\n\n";
                    output.open(dummy.fileNameResults.c_str(), ofstream::out | ofstream::app);
                    output.close();
                    int counter = 0;
                    int sizeSubLoop = 10;
                    dummy.clearedBatch.resize(sizeSubLoop,0);
                    for (int h = 0 ; h < 10 ; h++) {
                        cout << "Run " << h+1 << ": ";
                        dummy.originalRankings.clear();
                        dummy.originalRankings=dataRankings;
                        dummy.randomProfile(virtualCandidatesQueue,dummy.originalRankings);
                        dummy.firstStep=true;
                        dummy.stepOne(year_index);
                        dummy.clearedBatch[counter]=dummy.howManyCleared;
                        counter++;
                        if (counter==sizeSubLoop) {
                            output.open(dummy.fileNameResults.c_str(), ofstream::out | ofstream::app);
                            for (int n = 0 ; n < sizeSubLoop ; n++) {
                                output  << dummy.clearedBatch[n] << "\n";
                            }
                            output.close();
                            dummy.clearedBatch.resize(sizeSubLoop,0);
                            counter = 0;
                        }
                    }
                }
                if (doingReshuffle) {
                    vector<vector<int> > positionsLinked;
                    positionsLinked.clear();
                    positionsLinked = dummy.linkedPositions(dummy.originalRankings);
                    dummy.fileNameResults = "Results-reshuffled_profiles-" + currentYear + ".txt";
                    output.open(dummy.fileNameResults.c_str(), ofstream::out | ofstream::app);
                    output.close();
                    cout << "\n\n\n*************************************\n*";
                    cout << "\t" << y+1999 << "\tProfiles reshuffled\n";
                    cout << "*************************************\n\n";
                    int counter = 0;
                    int sizeSubLoop = 5;
                    dummy.clearedBatch.resize(sizeSubLoop,0);
                    for (int h = 0 ; h < 10 ; h++) {
                        cout << "Run " << h+1 << ": ";
                        dummy.originalRankings.clear();
                        dummy.originalRankings=dataRankings;
                        dummy.reshuffleRankings(positionsLinked);
                        dummy.firstStep=true;
                        dummy.stepOne(year_index);
                        dummy.clearedBatch[counter]=dummy.howManyCleared;
                        counter++;
                        if (counter==sizeSubLoop) {
                            output.open(dummy.fileNameResults.c_str(), ofstream::out | ofstream::app);
                            for (int n = 0 ; n < sizeSubLoop ; n++) {
                                output  << dummy.clearedBatch[n] << "\n";
                            }
                            output.close();
                            dummy.clearedBatch.resize(sizeSubLoop,0);
                            counter = 0;
                        }
                    }
                }
                (void) time(&endSimulations);
                cout << "\n\n\n*************************************\n*";
                cout << "*\t\t" << dummy.computeTime(startSimulations, endSimulations) << "\n";
                cout << "\n\n\n*************************************\n*";
            }
        } else {
            dummy.loadData(filename.c_str());
            dummy.firstStep=true;
            dummy.stepOne(year_index);
            dummy.firstStep=false;
            dummy.stepTwoReal();
        }
        (void) time(&t2);
        if (year_index>-1) {
            output.open(dummy.fileNameResults, std::ofstream::out | std::ofstream::app);
            output << dummy.computeTime(t1, t2);
            output.close();
        }
        cout << dummy.computeTime(t1, t2);
    }
    (void) time(&bigEnd);
    
    cout << "Total " << dummy.computeTime(bigStart, bigEnd);
    
    return 0;
}
