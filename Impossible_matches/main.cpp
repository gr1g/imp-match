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


using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using namespace std;

//int main(int argc, const char * argv[])
int main() {
    int year_index;
    string filename;
    time_t t1,t2;    // time variables to count the time spent running the program
    string second, minute, hour;
    long hours,minutes,seconds;
    second = "second";
    minute = "minute";
    hour = "hour";
    filename="none";
    ofstream output;
    
    
    year_index=0;
    
    for (year_index = 1 ; year_index < 12 ; year_index++) {
        if (year_index == 0) {
            filename = "1999-data.txt";
        }
        if (year_index == 1) {
            filename = "2000-data.txt";
        }
        if (year_index == 2) {
            filename = "2001-data.txt";
        }
        if (year_index == 3) {
            filename = "2002-data.txt";
        }
        if (year_index == 4) {
            filename = "2003-data.txt";
        }
        if (year_index == 5) {
            filename = "2004-data.txt";
        }
        if (year_index == 6) {
            filename = "2005-data.txt";
        }
        if (year_index == 7) {
            filename = "2006-data.txt";
        }
        if (year_index == 8) {
            filename = "2007-data.txt";
        }
        if (year_index == 9) {
            filename = "2008-data.txt";
        }
        if (year_index == 10) {
            filename = "2009-data.txt";
        }
        if (year_index == 11) {
            filename = "2010-data.txt";
        }
        
        DUMMY dummy;
        std::string currentYear = std::to_string(year_index+1999);
        dummy.fileNameResults = "Results-" +(currentYear) + ".txt";
        dummy.fileNameMatchings = "Matchings-" +(currentYear) + ".txt";
        
        output.open(dummy.fileNameResults, ofstream::out);
        
//        output.open(fileNameResults, ofstream::out | ofstream::app);
        //output.open("results.txt", ofstream::out | ofstream::app);
        
        cout << "\n\n";
        cout << "*************************************\n*\n*";
        cout << "\tYear " << year_index+1999 << "\n*\n*";
        cout << "*************************************\n\n";
        
        output << "\n\n";
        output << "*************************************\n*\n*";
        output << "\tYear " << year_index+1999 << "\n*\n*";
        output << "*************************************\n\n";
        
        
        output.close();
        
        
        (void) time(&t1);
        
        
        
        
        dummy.year = year_index+1999;
        dummy.loadData(filename);
        dummy.stepOne();
        dummy.stepTwo();
        
        
        (void) time(&t2);
        hours = (t2-t1)/3600;
        minutes = (t2-t1)/60-60*hours;
        seconds = t2-t1-60*minutes-3600*hours;
        if (minutes > 1 ) {
            minute = "minutes";
        }
        if (seconds > 1 ) {
            second = "seconds";
        }
        if (hours > 1 ) {
            hour = "hours";
        }
        
        
        output.open(dummy.fileNameResults, std::ofstream::out | std::ofstream::app);

        cout << "\nTime spent: ";
        if (hours < 1) {
            if (minutes < 1 ) {
                cout << seconds << " " << second << endl;
                output << seconds << " " << second << endl;
            } else {
                cout << minutes << " " << minute << ", " << seconds << " " << second <<endl;
                output << minutes << " " << minute << ", " << seconds << " " << second <<endl;
            }
        } else {
            cout << hours << " " << hour << ", " << minutes << " " << minute << ", " << seconds << " " << second <<endl;
            output << hours << " " << hour << ", " << minutes << " " << minute << ", " << seconds << " " << second <<endl;
        }
        output.close();
    }
    return 0;
}