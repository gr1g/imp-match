//
//  dummy.cpp
//  Impossible_matches
//
//  Created by Guillaume Haeringer on 4/28/14.
//  Copyright (c) 2014 Guillaume Haeringer. All rights reserved.
//

#include "dummy.h"

#define COLS 20

DUMMY::DUMMY()
{
	numGone=0;
}


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
// rankings[j][0] is the ID # of a position, could be the same as that of a candidate

int DUMMY::IndexInRankings(const vector<int> &S, int a){
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
bool DUMMY::CandidateisPresent(const vector<int> &S, int z) {
	bool found = (std::find(S.begin(), S.end(), z) != S.end());
	if (found==1) {
		if (Index(S,z)!=0) {
			return 1;
		} else {
			return 0;
		}
	} else {
		return 0;
	}
}



//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Maximum Matching sub-Routine                                     //
//                                                                      //
//                                                                      //


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
	////cout << "Matching done \n";
	return(ans);
}

//                                                                      //
//     Maximum Matching sub-Routine                                     //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////



void DUMMY::AdmissibleGraph() {
	edge.clear();
	edge.resize(Candidates.size());
	for (int i = 0 ; i < Candidates.size() ; i++) {
		for (int j = 0 ; j < rankings.size(); j++) {
			edge[i].push_back(0);
		}
	}
	int i,j;
	unsigned long stop=1;
	for (j = 0 ; j < rankings.size() ; j++) {
		if (rankings[j][1]>-1) { // so rankings[j][1] = 0, 1, 2, ... => the 1st ranked is a candidate (0, 1, 2, ... is the name of the candidate
			edge[Index(Candidates, rankings[j][1])][j]=1; // edge created between that candidate and that department.
			// edge[i][j] = 1 => edge between the i-th candidate in Candidates and the j-th position.
		}
	}
	for (j = 0 ; j < rankings.size(); j++) {
		stop=rankings[j].size()-1;
		for (int i = 1 ; i < rankings[j].size(); i++) {
			if (matching_candidate[Index(Candidates, rankings[j][i])]==-1) {
				stop=i;
				i=(int)rankings[j].size();
			}
		}
		for (i = 1 ; i < stop+1; i++) {
			edge[Index(Candidates, rankings[j][i])][j]=1;
		}
	}
}


//
//	loadData: load the data and identify the impossible matches.
//
void DUMMY::loadData(std::string dataFile) {
	
	original_rankings.clear();
	impossibles.clear();
	
	ifstream data_read;
	
	data_read.open(dataFile, ios::in | ios::out ); // opens the file
	if(!data_read) { // file couldn't be opened
		cerr << "Error: data file  could not be opened" << endl;
		exit(1);
	}
	size_t lines = 0;
	cout << "there are " << lines << " lines\n";
	ifstream in(dataFile);
	for (string s; getline(in,s); ) {
		cout << s << "\n";
		++lines;
	}
	cout << "there are " << lines << " lines\n";
	original_rankings.resize(lines);
	impossibles.resize(lines);
	for (int i = 0 ; i < original_rankings.size() ; i++) {
		original_rankings[i].resize(0);
		impossibles[i].resize(0);
	}
	int row = 0; // Row counter
	while (!data_read.eof()) {
		int data;
		for (int i = 0 ; i < 20 ; i++) {
			data_read >> data;
			if (data!=-1) {
				original_rankings[row].push_back(data);
			}
		}
		row++;
	}
	data_read.clear();
	data_read.close();
	
	cout << "data read " << original_rankings.size() << "\n";
	
	//
	//  Declare impossibles (-1=don't know, 0=not impossible, 1=impossible
	//
	for (int i = 0 ; i < original_rankings.size(); i++) {
		impossibles[i].push_back(original_rankings[i][0]);
		impossibles[i].push_back(0);
		for (int j = 2 ; j < original_rankings[i].size(); j++) {
			impossibles[i].push_back(-1);
		}
	}
	
	
	//
	//  Construct set of candidates
	//
	Candidates.clear();
	Candidates.resize(0);
	for (int i = 0 ; i < original_rankings.size() ; i++) {
		for (int j = 1 ; j < original_rankings[i].size() ; j++) {
			int found=0;
			for (int h = 0 ; h < Candidates.size(); h++) {
				if (Candidates[h]==original_rankings[i][j]) {
					found=1;
				}
			}
			if (found==0) {
				Candidates.push_back(original_rankings[i][j]);
			}
		}
	}
	
	/*
	 //
	 //  Construct set of positions
	 //
	 Positions.clear();
	 Positions.resize(0);
	 for (int i = 0 ; i < original_rankings.size() ; i++) {
		int found=0;
		for (int h = 0 ; h < Positions.size(); h++) {
	 if (Positions[h]==original_rankings[i][0]) {
	 found=1;
	 }
		}
		if (found==0) {
	 Positions.push_back(original_rankings[i][0]);
		}
		
	 }
	 */
	
	
	for (int j = 0 ; j < original_rankings.size() ; j++) {
		cout << original_rankings[j][0] << "\t";
		for (int i = 1 ; i < original_rankings[j].size() ; i++) {
			if (original_rankings[j][i]!=-1) {
				cout << original_rankings[j][i] << " ";
			}
		}
		cout << "\n";
	}
	
	//
	//	Find the impossible matches.
	//
//	for (int j = 15 ; j < original_rankings.size() ; j++) {
//		for (int i = 10 ; i < original_rankings[j].size() ; i++) {
	for (int j = 1 ; j < 16 ; j++) {
		for (int i = 3 ; i < original_rankings[j].size() ; i++) {
			if (impossibles[j][i]==-1) {
				cout << "\n********\nChecking for i0 = " << original_rankings[j][i] << ", s0 = " << original_rankings[j][0] << "\n";
				for (int ii = 0 ; ii < original_rankings[j].size() ; ii++) {
					cout << original_rankings[j][ii] << "(" << ii <<") ";
				}
				cout << "\n";
				cout << "j = " << j << ", i = " << i << "\n";
				CheckImpossible(original_rankings[j][i], j);
				cout << "\n";
			} else {
				// To be done
				cout << "Already decided\n";
			}
		}
	}
}



void DUMMY::ComprehensiveMatching(){
	matching_candidate.clear();
	matching_candidate.resize(Candidates.size());
	matching_department.clear();
	matching_department.resize(Positions.size());
	for (int i = 0 ; i < Candidates.size() ; i++) {
		matching_candidate[i]=-1;
	}
	for (int j = 0 ; j < rankings.size(); j++) {
		matching_department[j]=-1;
	}
	for (int i = 0 ; i < Candidates.size(); i++) {
		AdmissibleGraph();
		bipMatch();
	}
}


void DUMMY::Simplify_i0(int i0, int s0) {
	// i0 =
	
	Candidates.clear();
	Candidates.resize(0);
	Positions.clear();
	Positions.resize(0);
	
	for (int i = 1 ; i < (Index(original_rankings[s0], i0)+1) ; i++) {
		// index(i0)+1 to include i0, and nobody after.
		Candidates.push_back(original_rankings[s0][i]);
	}
	Positions.push_back(s0);
	// so Positions[0] = s0;
	
	// Candidates[] first contains the names of all candidates better ranked than i0 at s0
	// Positions[0] = s0 (more precisely, the index of position s0).
	//
	// What we do below:
	// for each candidate i in Candidates[], find all candidates better ranked than i at any position where i is ranked,
	//	then add those candidates to Candidates[]
	// then add the position to Positions[].
	// Positions[j] = the j-th entry in Positions, the index of the position in original_rankings.
	// so, if the 3rd entry in Positions[] is the 17th position in original_rankings[][], whose name is 84 we have
	// Positions[3]=17, and original_rankings[17][0]=84.
	
	for (int i = 0 ; i < Candidates.size() ; i++) {
		if (Candidates[i]!=i0) {
			for (int j = 0 ; j < original_rankings.size() ; j++) {
				if (CandidateisPresent(original_rankings[j], Candidates[i])) {
					for (int ii = 1 ; ii < Index(original_rankings[j], Candidates[i]) ; ii++) {
						bool CandIsPresent = (std::find(Candidates.begin(), Candidates.end(), original_rankings[j][ii]) != Candidates.end());
						if (!CandIsPresent) {
							Candidates.push_back(original_rankings[j][ii]);
						}
					}
					bool DepIsPresent = (std::find(Positions.begin(), Positions.end(), j) != Positions.end());
					if (!DepIsPresent) {
						Positions.push_back(j);
					}
				}
			}
		}
	}
	
	
	rankings.clear();
	rankings.resize(Positions.size());
	
	
	// Construct rankings[][], a copy of original_rankings[][]
	// such that i0 is deleted everywhere except at s0.
	// and such that rankings[][] are only about the candidates in Candidates[] and the positions in Positions[]
	
	for (int j = 0 ; j < Positions.size() ; j++) {
		rankings[j].resize(0);
		rankings[j].push_back(original_rankings[Positions[j]][0]);
		//
		//	The 3rd entry in Positions is the 17th position in original_rankings[], whose name is 84:
		// rankings[3][0] = original_rankings[17][0] = 84
		//
		//	Next we fill rankings[j][1], rankings[j][2], etc, but only with candidates in Candidates[].
		//	If j != 0 (it's not j0, then we skip i0
		//	If j=0 we don't skip i0.
		for (int i = 1 ; i < original_rankings[Positions[j]].size(); i++) {
			if (original_rankings[Positions[j]][i]!=i0) {
				bool CandIsPresent = (std::find(Candidates.begin(), Candidates.end(), original_rankings[Positions[j]][i]) != Candidates.end());
				if (CandIsPresent) {
					rankings[j].push_back(original_rankings[Positions[j]][i]);
				}
			} else {
				if (j==0) {
					bool CandIsPresent = (std::find(Candidates.begin(), Candidates.end(), original_rankings[Positions[j]][i]) != Candidates.end());
					if (CandIsPresent) {
						rankings[j].push_back(original_rankings[Positions[j]][i]);
					}
				}
			}
		}
	}
	
	cout << "\nSimplified ranking:\n";
	for (int j = 0 ; j < rankings.size() ; j++) {
		for (int i = 0 ; i < rankings[j].size() ; i++) {
			cout << rankings[j][i] << "\t";
		}
		cout << "\n";
	}
	
	
}



void DUMMY::CheckComprehensiveness(){
	
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

void DUMMY::CheckImpossible(int i0, int s0) {
	
	// i0 comes from original_rankings[j][i] = name of a candidate
	// s0 comes from j = index of the deparment in original_rankings[][] (not the name!)
	
	int decided=0;
	
	Simplify_i0(i0, s0);
	cout << "simplify done\n";
	ComprehensiveMatching();
	cout << "comprehensive matching done\n";
	CheckComprehensiveness();
	cout << "Check comprehensiveness done\n";
	
	cout << "there are " << rankings.size() << "\n";
	for (int j = 0 ; j < rankings.size() ; j++) {
		cout << rankings[j][0] << ", ";
	}
	cout << "\n";
	cout << "Matching: ";
	for (int i = 0 ; i < Candidates.size() ; i++) {
		if (matching_candidate[i]!=-1) {
			cout << Candidates[i] << " ";
//			cout << " (" << rankings[matching_candidate[i]][0] << "), ";
//		} else {
//			cout << " (" << matching_candidate[i] << "), ";
		}
	}
	cout << "\n";

	
	Gamma_dep.clear();
	Gamma_dep.resize(0);
	Gamma_cand.clear();
	Gamma_cand.resize(0);
	
	if (matching_candidate[Index(Candidates, i0)]==0) {
		// so i0 is matched to Positions[0]
		// but Positions[0] is s0.
		// so i0 is not an impossible match for s0.
		decided=1;
		impossibles[s0][Index(original_rankings[s0], i0)]=0;
		// recall that s0 is the s0-th department in original_rankings.
	}
	if (decided==0) {
		// i0 is not matched to s0 => don't know yet if impossible.
		buildJ_0(i0, 0);
		cout << "build J0 done\n";
		cout << "J_0 = ";
		for (int i = 0 ; i < J_0.size() ; i++) {
			cout << J_0[i] << ", ";
		}
		cout << "\n\n";
		if (J_0.size()==0) {
			// shouldn't happen. If J0 is empty then i0 should be matched to s0 at the comprehensive matching.
			cout << "\n\nERROR --- J0 empty\n\n";
			exit(1);
		}
		buildGamma(J_0,0);
		cout << "build Gamma done\n";
	}
	while (decided==0) {
		
		// Eliminate s0 from Gamma_dep: we want to reserve it for i0.
		// Only candidates in Gamma_cand will be matched
		// and i0 is NOT in Gamma_cand
		Gamma_dep.erase(Gamma_dep.begin()+0);
		
		if (Gamma_dep.size()>Gamma_cand.size()) {
			cout << "Gamma dep > Gamma_cand ";
			exit(1);
			for (int i= 0; i < Candidates.size() ; i++) {
				for (int j = 0 ; j < rankings.size() ; j++) {
					edge[i][j]=0;
				}
			}
			for (int j = 0 ; j < Positions.size(); j++) {
				bool DepIsPresent = (std::find(Gamma_dep.begin(), Gamma_dep.end(), j) != Gamma_dep.end());
				if (DepIsPresent) {
					for (int i = 1 ; i < rankings[j].size(); j++) {
						bool CandIsPresent = (std::find(Gamma_cand.begin(), Gamma_cand.end(), rankings[j][i]) != Gamma_cand.end());
						if (!CandIsPresent) {
							rankings[j].erase(rankings[j].begin()+(i));
						}
					}
				} else {
					rankings[j].erase(rankings[j].begin()+1, rankings[j].end());
				}
			}
			for (int i = 0 ; i < Candidates.size() ; i++) {
				matching_candidate[i]=-1;
			}
			for (int j = 0 ; j < Positions.size() ; j++) {
				matching_department[j]=-1;
			}
			ComprehensiveMatching();
			buildGamma(J_0,s0);
			bipMatch();
			int count=0;
			for (int i = 0 ; i < Candidates.size(); i++) {
				if (matching_candidate[i]!=-1) {
					count++;
				}
			}
			if (count==Gamma_cand.size()) {
				//cout << "Not dummy : could be matched\n";
				decided=1;
				//cout << "Decided = " << decided << "\n";
			}
		}
		
		buildPossible_K();
		sizeK = (Gamma_cand.size()-Gamma_dep.size());

		if (possible_K.size()>0) {
			cout <<  " --- combinations: C(" << possible_K.size() << "," << sizeK << ") = ";
			cout << coeff((int) possible_K.size(), (int) sizeK) << " ";
		} else {
			// K is necessarily empty, so i0 is an impossible match.
			impossibles[s0][Index(original_rankings[s0], i0)]=1;
			break;
		}
		if (possible_K.size()<sizeK) {
			//impossibles[s0][Index(original_rankings[s0], i0)]=1;
			cout << "\n\nError --- K too small\n";
			cout << "i0 = " << i0 << "\n";
			cout << "s0 = " << rankings[0][0] << "\n";
			cout << "\nRankings\n";
			for (int j = 0 ; j < rankings.size() ; j++) {
				for (int i = 0 ; i < rankings[j].size() ; i++) {
					cout << rankings[j][i] << ", ";
				}
				cout << "\n";
			}
			cout << "\n";
			cout << "|K| = " << sizeK << "\n";
			cout << "Possible K = ";
			for (int i = 0 ; i < possible_K.size() ; i++) {
				cout << possible_K[i] << ", ";
			}
			cout << "\n\n";
			cout << "Gamma_cand = ";
			for (int i = 0 ; i < Gamma_cand.size() ; i++) {
				cout << Gamma_cand[i] << ", ";
			}
			cout << "\n\n";
			cout << "|Gamma_dep| = " << Gamma_dep.size() << " = " << rankings.size() << " #depts\n";
			cout << "|Gamma_cand| = " << Gamma_cand.size() << "\n";
			cout << "|Possible_K| = " << possible_K.size() << "\n";
			exit(1);
		} else {
			signed long int counter=1;
			// Now going for the loops
			
			int stop=0;
			while (stop==0) {
				K.clear();
				K.resize(0);
				for (int j = 0 ; j < sizeK ; j++) {
					K.push_back(possible_K[SelectedIndex[j]]);
				}
				GraphOfJ();
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
				bipMatch();
				int count=0;
				for (int i = 0 ; i < Candidates.size(); i++) {
					if (matching_candidate[i]!=-1) {
						count++;
					}
				}
				if (count==Candidates.size()) {
					//cout << "NOT IMPOSSIBLE MATCH!!!!!!!!!!!!! i0 = " << i0 << ", s0 " << s0 << "(" << original_rankings[s0][0] << ")\n";
					impossibles[s0][Index(original_rankings[s0], i0)]=0;
					decided=1;
					stop=1;
					break;
				} else {
					if ((int) coeff((int) possible_K.size(), (int) sizeK)==(int)counter) {
						impossibles[s0][IndexInRankings(original_rankings[s0], i0)]=1;
						for (int ii = IndexInRankings(original_rankings[s0], i0)+1 ; ii < original_rankings[s0].size() ; ii++) {
						 bool candidatepresent = (std::find(Gamma_cand.begin(), Gamma_cand.end(), original_rankings[s0][ii]) != Gamma_cand.end());
						 if (!candidatepresent) {
							 impossibles[s0][ii]=1;
							 cout << "\nNext candidate " << original_rankings[s0][ii] << " also impossible\n";
						 }
						}
						stop=1;
						decided=1;
						break;
					} else {
						next_subset(SelectedIndex, (int) possible_K.size(), (int) sizeK);
						counter++;
					}
				}
			}
		}
	}
	//return 100;
}

void DUMMY::rematch(int i0, int s0){
	
	for (int j = 0 ; j < Positions.size() ; j++) {
		for (int i = 0 ; i < Candidates.size() ; i++) {
			edge[i][j]=0;
		}
	}
	
	for (int i = 0 ; i < Gamma_cand.size() ; i++) {
		for (int j = 0 ; j < Positions.size(); j++) {
			if (CandidateisPresent(rankings[j], Gamma_cand[i])) {
				edge[Index(Candidates, Gamma_cand[i])][j]=1;
			}
		}
	}
	
	edge[Index(Candidates, i0)][s0]=1;
	for (int i = 0 ; i < Candidates.size(); i++) {
		matching_candidate[i]=-1;
	}
	for (int j = 0 ; j < Positions.size(); j++) {
		matching_department[j]=-1;
	}
	
	bipMatch();
	
}

void DUMMY::buildPossible_K() {
	possible_K.clear();
	possible_K.resize(0);
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
		bool CandIsPresent = (std::find(possible_K.begin(), possible_K.end(), rankings[j][1]) != possible_K.end());
		if (CandIsPresent) {
			possible_K.erase(possible_K.begin()+(Index(possible_K, rankings[j][1])));
		}
	}
	SelectedIndex.clear();
	SelectedIndex.resize(0);
	
	for (int i = 0 ; i < possible_K.size() ; i++) {
		SelectedIndex.push_back(i);
	}
	// There's no the "+1" like in the paper because s0 eliminated from Gamma_dep[].
	// (In the paper, s0 is part of Gamma_dep[])
}


void DUMMY::buildJ_0(int i0, int s0){
	J_0.clear();
	J_0.resize(0);
	for (int i = 1 ; i < IndexInRankings(rankings[0], i0); i++) {
		J_0.push_back(rankings[0][i]);
	}
}

void DUMMY::buildGamma(vector<int> &J_0, int s0) {
	
	
	Gamma_dep.clear();
	Gamma_dep.resize(0);
	Gamma_cand.clear();
	Gamma_cand.resize(0);
	
	
	// Gamma_cand first contains all candidates in J_0, i.e., candidates ranked better than i0 at s0.
	for (int i=0 ; i < J_0.size(); i++) {
		Gamma_cand.push_back(J_0[i]);
	}
	// Then add all the candidates who are matched to some department.
	for (int i = 0 ; i < Candidates.size() ; i++) {
		if (matching_candidate[i]!=-1) {
			bool CandIsPresent = (std::find(Gamma_cand.begin(), Gamma_cand.end(), Candidates[i]) != Gamma_cand.end());
			if (!CandIsPresent) {
				Gamma_cand.push_back(Candidates[i]);
			}
		}
	}
	
	cout << "\n";
	cout << "\nGamma cand with just J0: ";
	for (int i = 0 ; i < Gamma_cand.size() ; i++) {
		cout << Gamma_cand[i] << ", ";
	}
	cout << "\n";
	vector<int> out_dep;
	out_dep.clear();
	out_dep.resize(0);
	vector<int> out_cand;
	out_cand.clear();
	out_cand.resize(0);
	
	// construction of out_dep[], a set of positions to withdraw.
	// out_dep[] is a set of positions' indices.
	// out_dep[j]=h => the j-th position in out_dep is the h-th position in rankings[],
	// whose name is rankings[h][0]=rankings[out_dep[j]][0]
	
	for (int j = 0 ; j < rankings.size(); j++) {
		if (matching_department[j]==-1) {
			out_dep.push_back(j);
		}
	}
	
	cout << "\nUnmatched depts: ";
	for (int j = 0 ; j < out_dep.size() ; j++) {
		cout << rankings[out_dep[j]][0] << " ";
	}
	cout << "\n";
	
	int found =1;
	int found2=0;
	
	while (found>0) {
		// Compute J^h
		// Add to out_cand[] the candidates in Gamma_cand[] who are acceptable for a department in out_dep[]
		for (int i = 0 ; i < Gamma_cand.size() ; i++) {
			for (int j = 0 ; j < out_dep.size() ; j++) {
				if (CandidateisPresent(rankings[out_dep[j]], Gamma_cand[i])) {
					bool CandIsPresent = (std::find(out_cand.begin(), out_cand.end(), Gamma_cand[i]) != out_cand.end());
					if (!CandIsPresent) {
						out_cand.push_back(Gamma_cand[i]);
						cout << "adding cand " << Gamma_cand[i] << " because acceptable for " << rankings[out_dep[j]][0] << "\n";
					}
				}
			}
		}
		// eliminate from Gamma_cand[] all the candidates in out_cand[]
		for (int i = 0 ; i < out_cand.size() ; i++) {
			bool CandIsPresent = (std::find(Gamma_cand.begin(), Gamma_cand.end(), out_cand[i]) != Gamma_cand.end());
			if (CandIsPresent) {
				Gamma_cand.erase(Gamma_cand.begin()+Index(Gamma_cand, out_cand[i]));
			}
			bool CandIsPresent2 = (std::find(J_0.begin(), J_0.end(), out_cand[i]) != J_0.end());
			if (CandIsPresent2) {
				J_0.erase(J_0.begin()+(Index(J_0, out_cand[i])));
			}
		}
		
		//////// Compute S^h
		//
		found2=0;
		// search a position that is not yet in out_dep
		for (int j = 0 ; j < rankings.size(); j++) {
			bool DepIsPresent = (std::find(out_dep.begin(), out_dep.end(), j) != out_dep.end());
			if (!DepIsPresent) {
				// the j-th position in rankings[][] is not in out_dep.
				if (matching_department[j]!=-1) {
					// that position is matched to the matching_department[j]-th candidate in Candidates.
					// The name of that candidate is Candidates[matching_department[j]]
					bool candidatepresent = (std::find(Gamma_cand.begin(), Gamma_cand.end(), Candidates[matching_department[j]]) != Gamma_cand.end());
					if (!candidatepresent) {
						// that candidate is not in Gamma_cand, so we add the position to out_dep[]
						out_dep.push_back(j);
						cout << "adding dep " << rankings[j][0] << " matched to \n";
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
	
	
	cout << "J^h done\n";
	
	Gamma_dep.clear();
	Gamma_dep.resize(0);
	for (int j = 0 ; j < rankings.size(); j++) {
		bool DepIsPresent = (std::find(out_dep.begin(), out_dep.end(), j) != out_dep.end());
		if (!DepIsPresent) {
			Gamma_dep.push_back(j);
		}
	}
	
	// Check that s0 is not in out_dep
	bool DepIsPresent = (std::find(out_dep.begin(), out_dep.end(), s0) != out_dep.end());
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
	//		Gamma_cand only contains matched candidates. We need to add unmatched candidates. 
	//		Loop below: Should only add unmatched candidates that are ranked better than a candidate in Gamma.
	//		But as we add new candidates in Gamma, new potential candidates for Gamma may show up. Loop stops
	//		when we don't add new candidates.
	
	found =1;
	while (found>0) {
		found2=0;
		for (int j = 1; j < Gamma_dep.size() ; j++) {
			for (int i = 1 ; i < rankings[Gamma_dep[j]].size() ; i++) {
				if (matching_candidate[Index(Candidates, rankings[Gamma_dep[j]][i])]==-1) {
					bool CandIsPresent = (std::find(Gamma_cand.begin(), Gamma_cand.end(), rankings[Gamma_dep[j]][i]) != Gamma_cand.end());
					if (!CandIsPresent) {
						Gamma_cand.push_back(rankings[Gamma_dep[j]][i]);
						found2 = 1;
					}
				}
			}
		}
		if (found2==0) {
			found=0;
		}
	}
	
	cout << "adding unmatched candidates\n";
	
	// Check that all candidates in J0 are also acceptable for a position different from s0
	// if there is a candidate for which it is not the case (i.e., ranked only by s0)
	// then candidates below (including i0) are impossible for s0
	// and thus in the reduced profile (original_rankings[s0][]) i0 is not ranked by s0.
	
	
	for (int i = 0 ; i < J_0.size() ; i++ ) {
		found = 0;
		for (int j = 0 ; j < Gamma_dep.size() ; j++) {
			bool CandIsPresent = (std::find(rankings[Gamma_dep[j]].begin(), rankings[Gamma_dep[j]].end(), J_0[i]) != rankings[Gamma_dep[j]].end());
			if (CandIsPresent) {
				Gamma_cand.push_back(rankings[Gamma_dep[j]][i]);
				found = 1;
			}
		}
		if (found==0) {
			cout << "\n\nError --- candidate in J_0 only acceptable for s0\n\n";
			exit(1);
		}
	}
	
	cout << "Check all J0 acceptable\n";
	
	//
	//	Construct the connected component of Gamma_cand that includes J0
	//	(i.e., there exists a path from J_0 to any other candidate in Gamma
	//
	//	Then reconstruct Gamma_dep
	
	
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
	Gamma_dep.push_back(0);
	
	
	cout << "Ready to construct the connected compenent\n";

	cout << "Gamma dep = ";
	for (int j = 0 ; j < Gamma_dep.size(); j++) {
		cout << rankings[Gamma_dep[j]][0] << " ";
	}
	cout << "\n";
	
	found=1;
	while (found>0) {
		for (int j = 0 ; j < Gamma_dep.size(); j++) {
			found2=0;
			for (int i = 1 ; i < rankings[Gamma_dep[j]].size(); i++) {
				// Find candidates that are in temp_Gamma_cand that are acceptable for a department in Gamma_dep
				bool CandIsPresent = (std::find(temp_Gamma_cand.begin(), temp_Gamma_cand.end(), rankings[Gamma_dep[j]][i]) != Gamma_cand.end());
				cout << "Candidat " << rankings[Gamma_dep[j]][i] << "\n";
				if (CandIsPresent) {
					// Now check if that candidate is already in Gamma_cand. If not add it to Gamma_cand.
					bool CandIsPresent2 = (std::find(Gamma_cand.begin(), Gamma_cand.end(), rankings[Gamma_dep[j]][i]) != Gamma_cand.end());
					if (!CandIsPresent2) {
						cout << "\tadded to Gamma_cand (" << rankings[Gamma_dep[j]][i] << ")\n";
						Gamma_cand.push_back(rankings[Gamma_dep[j]][i]);
						found2=1;
						// Now add all departments where candidate rankings[Gamma_dep[j]][i] is acceptable
						int TheCurrentCandidate = rankings[Gamma_dep[j]][i];
						for (int jj = 0 ; jj < rankings.size() ; jj++) {
							// Find positions that rank the candidate
							bool CandIsPresent3 = (std::find(rankings[jj].begin(), rankings[jj].end(), TheCurrentCandidate) != rankings[jj].end());
							if (CandIsPresent3) {
								// Check if the position where he's ranked is in Gamma_dep. If not, add the position to Gamma_dep.
								bool DepIsPresent = (std::find(Gamma_dep.begin(), Gamma_dep.end(), jj) != Gamma_dep.end());
								if (!DepIsPresent) {
									Gamma_dep.push_back(jj);
									found2 =1;
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
	
	cout << "Connected compenent done\n";

	
	cout << "\n\nGamma_dep size = " << Gamma_dep.size() << "\n";
	cout << "\n\nGamma_cand size = " << Gamma_cand.size() << "\n";
	
	//	TO BE DONE
	
	
	// Check J_0 is not empty
	if (J_0.size()==0) {
		cout << "\n\nError: J_0 empty after Gamma\n";
		exit(1);
	}
	
	if (Gamma_dep.size() > Gamma_cand.size()) {
		cout << "error size, Gamma_dep > Gamma_cand\n";
		cout << "Position " << rankings[s0][0] << "\n";
		cout << "Gamma_cand (" << Gamma_cand.size() << "): ";
		for (int i = 0 ; i < Gamma_cand.size(); i++) {
			cout << Gamma_cand[i] << " ";
		}
		cout << "\n";
		cout << "Gamma_dep (" << Gamma_dep.size() << "): ";
		for (int i = 0 ; i < Gamma_dep.size(); i++) {
			cout << rankings[Gamma_dep[i]][0] << " ";
		}
		cout << "\n";
		cout << "Position " << rankings[s0][0] << "\n";
		cout << "Gamma_cand (" << Gamma_cand.size() << "): ";
		for (int i = 0 ; i < Gamma_cand.size(); i++) {
			cout << Gamma_cand[i] << " ";
		}
		cout << "\n";
		cout << "Gamma_dep (" << Gamma_dep.size() << "): ";
		for (int i = 0 ; i < Gamma_dep.size(); i++) {
			cout << rankings[Gamma_dep[i]][0] << " ";
		}
		cout << "\n";
		
		
		
		cout << "\nEXIT\n";
		
		cout << "Gamma dep > Gamma cand\n";
		cout << "for position " << " " << original_rankings[s0][0] << "\n";
		exit(1);
	}
	
	
	
	
}

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
			//cout << current[j] << "\t";
		}
		//cout << prev[i+1];
		//cout << "\n";
		for (j = 0 ; j < n ; j++) {
			prev[j]=current[j];
		}
	}
	//cout << "C(" << n << "," << k << ") = " << current[k];
	
	if (n==k) {
		return (1);
	}
	
	return current[k];
	
}



void DUMMY::GraphOfJ(){
	
	// Given rankings[][] and K[]
	// construct a graph such that it contains only students in Gamma_cand
	// when rankings[][] is truncated at K[];
	
	// Declare temp_rankings;
	vector<vector<int> > temp_rankings;
	temp_rankings.clear();
	temp_rankings.resize(Gamma_dep.size());
	for (int j = 0 ; j < Gamma_dep.size(); j++) {
		temp_rankings[j].resize(0);
	}
	
	temp_rankings[0].push_back(rankings[Gamma_dep[0]][0]);
	for (int i = 1 ; i < rankings[Gamma_dep[0]].size() ; i++) {
		if (CandidateisPresent(Gamma_cand, rankings[Gamma_dep[0]][i])) {
			temp_rankings[0].push_back(rankings[Gamma_dep[0]][i]);
		}
	}
	
	
	for (int j = 1 ; j < Gamma_dep.size() ; j++) {
		temp_rankings[j].push_back(rankings[Gamma_dep[j]][0]);
		for (int i = 1 ; i < rankings[Gamma_dep[j]].size() ; i++) {
			if (CandidateisPresent(Gamma_cand, rankings[Gamma_dep[j]][i])) {
				temp_rankings[j].push_back(rankings[Gamma_dep[j]][i]);
			}
		}
	}
	
	
	// truncating at K
	for (int j = 0 ; j < temp_rankings.size() ; j++) {
		for (int i = 1 ; i < temp_rankings[j].size() ; i++) {
			if (CandidateisPresent(K, temp_rankings[j][i])) {
				temp_rankings[j].erase(temp_rankings[j].begin()+(IndexInRankings(temp_rankings[j], temp_rankings[j][i])));
				i = (int)(temp_rankings[j].size());
			}
		}
	}
	
	rankings.clear();
	rankings.resize(temp_rankings.size());
	for (int j = 0 ; j < rankings.size() ; j++) {
		rankings[j].resize(0);
		for (int i = 0 ; i < temp_rankings[j].size() ; i++) {
			rankings[j].push_back(temp_rankings[j][i]);
		}
	}
	
	// reconstructing the set candidates[]
	// those are the admissible candidates after the truncation
	// and the graph
	edge.clear();
	edge.resize(Candidates.size());
	for (int i = 0 ; i < Candidates.size() ; i++) {
		for (int j = 0 ; j < rankings.size(); j++) {
			edge[i].push_back(0);
		}
	}
	
	Candidates.clear();
	Candidates.resize(0);
	for (int j = 1 ; j < rankings.size() ; j++) {
		for (int i = 1 ; i < rankings[j].size() ; i++) {
			bool isPresent = (std::find(Candidates.begin(), Candidates.end(), rankings[j][i]) != Candidates.end());
			if (!isPresent) {
				Candidates.push_back(rankings[j][i]);
			}
			edge[Index(Candidates, rankings[j][i])][j]=1;
		}
	}
	
	
	
	//
	//
	//
	//
	///////////////////
	
	
	//J_0.erase(J_0.begin()+(Index(J_0, cand)));
	
	temp_rankings.clear();
	
}


//////////////////////////////////////////////////////////////////////////
//                                                                      //
//     Enumerating all subsets of size k                                //
//                                                                      //
//                                                                      //


void DUMMY::next_subset(vector<int> &SelectedIndex, int n, int k) {
	int i,j;
	int movable[k];
	for (i = 0 ; i < k ; i++) {
		movable[k]=0;
	}
	for (i = 0 ; i < k ; i++) {
		if (SelectedIndex[i]==n-k+i) {
			movable[i]=0;
			//printf("%d not movable\n", selected_index[i]);
		} else {
			movable[i]=1;
			//printf("%d movable\n", selected_index[i]);
		}
	}
	if (movable[k-1]==1) {
		SelectedIndex[k-1]=SelectedIndex[k-1]+1;
	} else {
		for (i = 0 ; i < k ; i++) {
			if (movable[k-i-1]==1) {
				SelectedIndex[k-i-1]=SelectedIndex[k-i-1]+1;
				for (j = k-i; j < k ; j++) {
					SelectedIndex[j]=SelectedIndex[j-1]+1;
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






