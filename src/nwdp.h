/*
 * nwdp.h
 *
 *  Created on: Oct 28, 2013
 *      Author: Stefan Seemann, seemann@rth.dk
 */

#ifndef NWDP_H_
#define NWDP_H_

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>

#define DEBUG 0

using namespace std;

extern float pnullSgl[];
extern float psubSgl[];
extern float pnullDbl[];
extern float psubDbl[];

extern float maxpnullSgl;
extern float maxpnullDbl;

extern std::map<const char, const int> nucIdx;

/* arguments */
extern float kappa;
extern float alpha;
extern float beta;
extern float deltanull;

extern float INFINITE;

typedef struct tagLocalHit {
	float similarity;		/* Similarity */
	int lstart_1;			/* Start position of LS_a */
	int lend_1;				/* End position of LS_a */
	int lstart_2;			/* Start position of LS_b */
	int lend_2;				/* End position of LS_b */
} LocalHit;

extern int readinput( istream & is, string & name, string & seq, vector<float> & prob );
extern void getunpaired( vector<float> & prob, int len, vector<float> & probSgl );
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);

extern void reducematrix(vector<float> & prob, int len, int prec );
extern float nwdp( string seq_1, vector<float> & probDbl_1, vector<float> & probSgl_1, int relative_startindex_1, int * idx_1_aln, int len_1, string seq_2, vector<float> & probDbl_2, vector<float> & probSgl_2, int relative_startindex_2, int * idx_2_aln, int len_2, float * subprobDbl_1, float * subprobSgl_1, float * subprobDbl_2, float * subprobSgl_2, bool freeEndGaps, bool prm, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP );
extern int nwdp_seed( string seq_1, vector<float> & probDbl_1, vector<float> & probSgl_1, int * idx_1_aln, int len_1, string seq_2, vector<float> & probDbl_2, vector<float> & probSgl_2, int * idx_2_aln, int seedlen, int len_2, float * subprobDbl_1, float * subprobSgl_1, float * subprobDbl_2, float * subprobSgl_2, float ** sim, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP );
extern int nwseq_seed( string seq_1, string seq_2, float ** F, char ** trF );
extern void nwdp_initF( float ** F, int L1, int L2 );
extern void nwdp_initF_affinegaps( float ** F, int L1, int L2, bool local );
extern void nwdp_initGap( float ** Q, int L1, int L2 );
extern float nwdb_align_seq_sim( char nuc_1, float probSgl_1, char nuc_2, float probSgl_2 );
extern float nwdb_global_align_affinegaps( float* probDbl_1, float* probSgl_1, int len_1, float* probDbl_2, float* probSgl_2, int len_2, bool freeEndGaps, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP );
extern float max3( float f1, float f2, float f3, char* ptr );
extern float max( float f1, float f2 );
template <typename T>
extern void print_matrixdp( T ** F, float * prob_1, int len_1, float * prob_2, int len_2 );
extern float simalign_affinegaps( float ** Z, int len_1, int len_2, int * idx_1_aln, int * idx_2_aln, int & len_pair, int & len_aln, bool global2, bool prm, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP );
extern void affinegapcosts( int * idx_1_aln, int * idx_2_aln, int & len_aln, int & open, int & extended );
extern void nwdp_initTB( char ** traceback, int L1, int L2 );
extern void reverse( int * list, int len );
extern float getOverlap2ndInterval( int start_1, int end_1, int start_2, int end_2 );
extern void getlogoddsDblNeighborhood( vector<float> & probDbl, string seq, int len, float pnull, int radius, float theta );
extern void getlogoddsSgl( vector<float> & probSgl, string seq, int len, float pnull );
extern void printalign(string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln );
extern void freeMatrix(float ** matrix, int leny);

extern void usage_dotaligner(char * program);

#endif /* NWDP_H_ */
