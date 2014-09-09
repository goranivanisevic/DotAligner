/*
 * main_v2.cpp
 *
 *  Created on: Jul 03, 2014
 *      Author: Stefan Seemann, seemann@rth.dk
 *
 *  Implements a heuristic for aligning to base pair probability matrices (dotplots)
 *
 *  Based on Needleman-Wunsch algorithm implementation in http://www.rolfmuertter.com/code/nw.php
 *
 */

#include "nwdp_v2.h"
#include <getopt.h>
#include <map>
#include <ctime>
#include <cstdlib>
#include <algorithm>

using namespace std;

// global arguments
float kappa = 0.5;
float alpha = -0.2;
float beta  = -0.1;
float tau = 0.5;

float INFINITE = -1000;


int main( int argc, char ** argv )
{
		char * program = *argv ;
        vector<float> probDbl_1;
        vector<float> probDbl_2;
        vector<float> probSgl_1;
        vector<float> probSgl_2;
        int len_1 = 0, len_2 = 0, len_aln = 0;
        string name1, name2, seq_1, seq_2;
        int precision = 4;
        float pnull = 0.0005;
        int subopt = 10;

        /* arguments */
        string  filename1, filename2;

        /* boolean flags */
        static int help_flag = 0;

    	/* parsing long options */
    	while(1)
    	{
    		static struct option long_options[] =
    		{
    			{"dotplot", required_argument, 0, 'd'},			// input dot plots
    			{"kappa", required_argument, 0, 'k'},			// weight of sequence / unpaired probability similarity
    			{"alpha", required_argument, 0, 'a'},		 	// affine gap costs = alpha + k * beta (k gaps)
    			{"beta", required_argument, 0, 'b'},			// affine gap costs = alpha + k * beta (k gaps)
    			{"precision", required_argument, 0, 'p'},       // number of digits considered of base pair reliabilities
    			{"pnull", required_argument, 0, 'n'},			// minimal probability
    			{"subopt", required_argument, 0, 's'},			// number of suboptimal sequence/unpaired probability alignments examined
    			{"tau", required_argument, 0, 't'},				// weight of sequence similarity in sequence / unpaired probability similarity
    			{"help", no_argument, &help_flag, 1},
    			{0, 0, 0, 0}
    		};

    		/* print help text */
    		if( help_flag )
    			usage_dotaligner(program);

    		/* getopt_long stores the option index here. */
    		int option_index = 0;

    		int cmd = getopt_long(argc, argv, "d:k:a:b:p:n:s:t:", long_options, &option_index);

    		/* Detect the end of the options. */
    		if (cmd == -1)
    			break;

    		switch(cmd)
    		{
    			case 0:   break;
    			case 'd': if (!filename1.length())
    						  filename1 = optarg;
    			          else if(!filename2.length())
    						  filename2 = optarg;
    			          else
    			        	  usage_dotaligner(program);
    			          break;
    			case 'k': kappa = atof(optarg); break;
    			case 'a': alpha = atof(optarg); break;
    			case 'b': beta = atof(optarg); break;
    			case 'p': precision = atoi(optarg); break;
    			case 'n': pnull = atof(optarg); break;
    			case 's': subopt = atoi(optarg); break;
    			case 't': tau = atof(optarg); break;
    			default:  abort();
    			/* no break */
    		}
    	}

    	/* read input file */
        if( filename1.length() > 0 )
        {
            ifstream inputfile1;
            inputfile1.open(filename1.c_str());
            if( !inputfile1.is_open() ) {
                 cerr << "Could not open file " << filename1 << endl;
                 exit(EXIT_FAILURE);
            }
            len_1 = readinput(inputfile1, name1, seq_1, probDbl_1);
            inputfile1.close();
        }
        else {
        	usage_dotaligner(program);
        }
        if( filename2.length() > 0 )
        {
            ifstream inputfile2;
            inputfile2.open(filename2.c_str());
            if( !inputfile2.is_open() ) {
                 cerr << "Could not open file " << filename2 << endl;
                 exit(EXIT_FAILURE);
            }
            len_2 = readinput(inputfile2, name2, seq_2, probDbl_2);
            inputfile2.close();
        }
        else {
        	usage_dotaligner(program);
        }

        /* get unpaired probabilities */
        getunpaired(probDbl_1, len_1, probSgl_1);
        getunpaired(probDbl_2, len_2, probSgl_2);

        /* log-odds scores */
        getlogoddsSgl( probSgl_1, seq_1, len_1, pnull );
        getlogoddsSgl( probSgl_2, seq_2, len_2, pnull );
        getlogoddsDbl( probDbl_1, seq_1, len_1, pnull );
        getlogoddsDbl( probDbl_2, seq_2, len_2, pnull );

        /* reduce depth of probability matrices */
        reducematrix(probDbl_1, len_1*len_1, precision);
        reducematrix(probDbl_2, len_2*len_2, precision);
        reducematrix(probSgl_1, len_1, precision);
        reducematrix(probSgl_2, len_2, precision);
//		#if DEBUG
//        	for( int i=0; i<len_1; i++ ) { for( int j=0; j<len_1; j++ ) cerr << probDbl_1.at( i*len_1+j ) << "\t"; cerr << endl; }; cerr << endl;
//        	for( int i=0; i<len_2; i++ ) { for( int j=0; j<len_2; j++ ) cerr << probDbl_2.at( i*len_2+j ) << "\t"; cerr << endl; }; cerr << endl;
        	for( int i=0; i<len_1; i++ ) cerr << probSgl_1[ i ] << " "; cerr << endl;
        	for( int i=0; i<len_2; i++ ) cerr << probSgl_2[ i ] << " "; cerr << endl;
//		#endif

        /* initialize arrays of indices */
		int *idx_1_aln = new int[len_1];
		for( int i=0; i<len_1; i++ ) idx_1_aln[ i ] = i;
		int *idx_2_aln = new int[len_2];
		for( int i=0; i<len_2; i++ ) idx_2_aln[ i ] = i;
		int *idx_1_subaln = new int[len_1];
		for( int i=0; i<len_1; i++ ) idx_1_subaln[ i ] = i;
		int *idx_2_subaln = new int[len_2];
		for( int i=0; i<len_2; i++ ) idx_2_subaln[ i ] = i;

	    /* Initialize dynamic programming matrix F (fill in first row and column) */
	    float ** F = new float * [ len_2 + 1  ];
	    for( int j = 0; j <= len_2; j++ )
	    	F[ j ] = new float [ len_1 + 1 ];
	    /* Initialize Q and P matrices for cost of alignments that end with a gap */
	    float ** Q = new float * [ len_2 + 1  ];
	    for( int j = 0; j <= len_2; j++ )
	    	Q[ j ] = new float [ len_1 + 1 ];
	    float ** P = new float * [ len_2 + 1  ];
	    for( int j = 0; j <= len_2; j++ )
	      	P[ j ] = new float [ len_1 + 1 ];
	    /* Initialize traceback probability matrices */
	    char ** trF = new char * [ len_2 + 1 ];
	    for( int j = 0; j <= len_2; j++ )
	    	trF[ j ] = new char [ len_1 + 1 ];
	    char ** trQ = new char * [ len_2 + 1 ];
	    for( int j = 0; j <= len_2; j++ )
	    	trQ[ j ] = new char [ len_1 + 1 ];
	    char ** trP = new char * [ len_2 + 1 ];
	    for( int j = 0; j <= len_2; j++ )
	    	trP[ j ] = new char [ len_1 + 1 ];

		/* STEP 1: global alignment (Needleman-Wunsch algorithm) of sequences by scoring sequence similarity and unpaired probabilities
	       return optimal alignment in idx_1_aln and idx_2_aln */
		float score_seq = nwdp( seq_1, probSgl_1, idx_1_aln, len_1, seq_2, probSgl_2, idx_2_aln, len_2, len_aln, F, Q, P, trF, trQ, trP );

		/* STEP 3: get similarity of pairs of dotplots corresponding to optimal alignment */
		float score_dp = simdp( probDbl_1, len_1, probDbl_2, len_2, idx_1_aln, idx_2_aln, len_aln );
		float sim = kappa * score_seq + (1 - kappa) * score_dp;
		cout << "Score_seq: " << score_seq << endl;
		cout << "Score_dp: " << score_dp << endl;
		cout << "Len_aln: " << len_aln << endl;
		cout << "Score: " << sim << endl;
		printalign(seq_1, idx_1_aln, seq_2, idx_2_aln, len_aln);
		cout << "==========" << endl;

		float tempsim;
		int len_subaln;
   		srand ( time(NULL) );	/* initialize random seed */
		while( subopt ) {

			/* STEP 2: get suboptimal sequence/unpaired probability alignment by probabilistic backtracking through scoring matrix */
			prob_backtracking( idx_1_subaln, len_1, idx_2_subaln, len_2, len_subaln, F);
			score_seq = simbp( seq_1, probSgl_1, len_1, seq_2, probSgl_2, len_2, idx_1_subaln, idx_2_subaln, len_subaln );

			/* STEP 3: get similarity of pairs of dotplots corresponding to suboptimal alignment */
			score_dp = simdp( probDbl_1, len_1, probDbl_2, len_2, idx_1_subaln, idx_2_subaln, len_subaln );
			tempsim = kappa * score_seq + (1 - kappa) * score_dp;
			cout << "Score_seq: " << score_seq << endl;
			cout << "Score_dp: " << score_dp << endl;
			cout << "Len_aln: " << len_subaln << endl;
			cout << "Score: " << tempsim << endl;
			printalign(seq_1, idx_1_subaln, seq_2, idx_2_subaln, len_subaln);
			cout << "==========" << endl;

			if( tempsim > sim ) {
				sim = tempsim;
				len_aln = len_subaln;
				std::copy(idx_1_subaln, idx_1_subaln+len_1, idx_1_aln);
				std::copy(idx_2_subaln, idx_2_subaln+len_2, idx_2_aln);
			}

			/* increment number of suboptimals to be still examined */
			subopt--;
	    }

		/* OUTPUT */

		/* print aligned probabilities and similarity */
		cout << "Similarity = " << sim << ", Length_1 = " << len_1 << ", Length_2 = " << len_2 << ", Length_alignment = " << len_aln << endl;

		/* print sequences aligned by dot plot alignment */
		printalign(seq_1, idx_1_aln, seq_2, idx_2_aln, len_aln);

		/* free memory */
	    delete[] idx_1_aln;
	    delete[] idx_2_aln;
	    delete[] idx_1_subaln;
	    delete[] idx_2_subaln;
	    probDbl_1.clear();
        probDbl_2.clear();
        probSgl_1.clear();
        probSgl_2.clear();
        for( int j = 0; j <= len_2; j++ ) delete[] F[ j ]; delete[] F;
        for( int j = 0; j <= len_2; j++ ) delete[] Q[ j ]; delete[] Q;
        for( int j = 0; j <= len_2; j++ ) delete[] P[ j ]; delete[] P;
        for( int j = 0; j <= len_2; j++ ) delete[] trF[ j ]; delete[] trF;
        for( int j = 0; j <= len_2; j++ ) delete[] trP[ j ]; delete[] trP;
        for( int j = 0; j <= len_2; j++ ) delete[] trQ[ j ]; delete[] trQ;

        return  0;
}
