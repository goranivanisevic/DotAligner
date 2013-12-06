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

extern std::map<const char, const int> nucIdx;

/* arguments */
extern float gappenalty;
extern float alpha;
extern float extgaps;
extern float stemgaps;


extern int readinput( istream & is, string & name, string & seq, vector<float> & prob );
extern void getunpaired( vector<float> & prob, int len, vector<float> & probSgl );
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);

extern void reducematrix(vector<float> & prob, int len, int prec );
extern float nwdp( string seq_1, vector<float> & probDbl_1, vector<float> & probSgl_1, int startindex_1, int * idx_1_aln, int len_1, string seq_2, vector<float> & probDbl_2, vector<float> & probSgl_2, int startindex_2, int * idx_2_aln, int len_2, bool prm );
extern void nwdp_initF( float ** F, int L1, int L2 );
extern float nwdb_align( float ** F, char nuc_1, float probSgl_1, string seq_1, float* probDbl_1, int len_1, char nuc_2, float probSgl_2, string seq_2, float* probDbl_2, int len_2 );
extern float max3( float f1, float f2, float f3, char* ptr );
extern float min3( float f1, float f2, float f3, char* ptr );
extern float min( float f1, float f2 );
extern void print_matrixdp( float ** F, float * prob_1, int len_1, float * prob_2, int len_2 );
extern float nwdistance( float ** distance, int len_1, int len_2, int * idx_1_aln, int * idx_2_aln, int & len_1_aln, int & len_2_aln, int precision, bool global2, bool prm);
extern void nwdp_initTB( char ** traceback, int L1, int L2 );
extern void reverse( int * list, int len );
extern void printalign(string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln );
extern void freeMatrix(float ** matrix, int leny);

extern int nw(string, string, string&, string&, bool);
extern int nw_align(int **, char **, string, string, string&, string&, int);
extern void dpm_init        ( int **, char **, int, int, int );
extern void print_al        ( string&, string& );
extern void print_matrix    ( int ** const, string, string );
extern void print_traceback ( char ** const, string, string );

extern void usage_dotaligner(char * program);

#endif /* NWDP_H_ */
