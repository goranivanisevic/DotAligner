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
#include <cmath>

#define DEBUG 0
#define DPSUM 1

using namespace std;

extern int readinput ( istream & is, string & name, string & seq, vector<float> & prob );
vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);

extern void reducematrix(vector<float> & prob, int len, int prec );
extern float nwdp( vector<float> & prob_1, int startindex_1, int * idx_1_aln, int len_1, vector<float> & prob_2, int startindex_2, int * idx_2_aln, int len_2, float gappenalty, bool prm );
extern void nwdp_initF( float ** F, int L1, int L2, float gapenalty );
extern float nwdb_align( float ** F, float* prob_1, int len_1, float* prob_2, int len_2, int d );
extern float min3( float f1, float f2, float f3, char* ptr );
extern float min( float f1, float f2 );
extern void print_matrixdp( float ** F, float * prob_1, int len_1, float * prob_2, int len_2 );
extern float nwdistance( float ** distance, int len_1, int len_2, int * idx_1_aln, int * idx_2_aln, int & len_1_aln, int & len_2_aln, int precision, float gappenalty, bool prm);
extern void nwdp_initTB( char ** traceback, int L1, int L2 );
extern void reverse( int * list, int len );
extern void printalign(string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln );
extern void freeMatrix(float ** matrix, int leny);

extern int nw(string, string, string&, string&, bool);
extern int nw_align(int **, char **, string, string, string&, string&, int);
extern void  dpm_init        ( int **, char **, int, int, int );
extern void  print_al        ( string&, string& );
extern void  print_matrix    ( int ** const, string, string );
extern void  print_traceback ( char ** const, string, string );
extern int   max3            ( int, int, int, char * );


#endif /* NWDP_H_ */
