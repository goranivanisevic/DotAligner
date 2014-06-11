/*
 * main.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Stefan Seemann, seemann@rth.dk
 *
 *  Implements an extension of the Needleman-Wunsch algorithm
 *  for global nucleotide sequence alignment
 *  that aligns two dotplots by taking care of relationships between two indices (basepairs).
 *
 *  Based on http://www.rolfmuertter.com/code/nw.php
 *
 */

#include "nwdp.h"
#include <getopt.h>
#include <map>
#include <ctime>
#include <cstdlib>

using namespace std;

// global arguments
float kappa = 0.5;
float alpha = -0.2;
float beta  = -0.1;
//float deltanull = 0.5;

float INFINITE = -1000;


int main( int argc, char ** argv )
{
		char * program = *argv ;
        vector<float> probDbl_1;
        vector<float> probDbl_2;
        vector<float> probSgl_1;
        vector<float> probSgl_2;
        int len_1 = 0, len_2 = 0, len_pair = 0, len_aln = 0;
        string name1, name2, seq_1, seq_2;
        int precision = 4;
        int maxshift = INFINITE;
        float pnull = 0.0005;
        int radius = 0;
    	float theta = 0.25;
    	int seedlen = 15;
    	float deltanull = 0.5;

        /* arguments */
        string  filename1, filename2;

        /* boolean flags */
        static int setglobal2_flag = 0;
        static int setseqaln_flag = 0;
        static int help_flag = 0;

    	/* parsing long options */
    	while(1)
    	{
    		static struct option long_options[] =
    		{
    			{"global2", no_argument, &setglobal2_flag, 1},  // global alignment in step 2 (default is local alignment in step 2)
    			{"dotplot", required_argument, 0, 'd'},			// input dot plots
    			{"kappa", required_argument, 0, 'k'},			// weight of sequence similarity
    			{"alpha", required_argument, 0, 'a'},		 	// affine gap costs = alpha + k * beta (k gaps)
    			{"beta", required_argument, 0, 'b'},			// affine gap costs = alpha + k * beta (k gaps)
    			{"precision", required_argument, 0, 'p'},       // number of digits considered of base pair reliabilities
    			{"maxshift", required_argument, 0, 'm'},		// maximal shift of two input sequences in the final alignment (default 50% of longer sequence)
    			{"pnull", required_argument, 0, 'n'},			// minimal probability
    			{"radius", required_argument, 0, 'r'},			// consider diagonal (stem) neighborhood of specific radius (default 0)
    			{"theta", required_argument, 0, 't'},			// weight of neighborhood <radius> in base pair probability
    			{"deltanull", required_argument, 0, 'l'},		// value of tau if both probabilities are zero
    			{"seedlen", required_argument, 0, 's'},         // length of seed alignments (default 20 nucleotides)
    			{"seqaln", no_argument, &setseqaln_flag, 1},	// run local sequence alignment to position both sequences
    			{"help", no_argument, &help_flag, 1},
    			{0, 0, 0, 0}
    		};

    		/* print help text */
    		if( help_flag )
    			usage_dotaligner(program);

    		/* getopt_long stores the option index here. */
    		int option_index = 0;

    		int cmd = getopt_long(argc, argv, "d:k:a:b:p:m:n:r:t:l:s:", long_options, &option_index);

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
    			case 'm': maxshift = atoi(optarg); break;
    			case 'n': pnull = atof(optarg); break;
    			case 'r': radius = atoi(optarg); break;
    			case 't': theta = atof(optarg); break;
    			case 'l': deltanull = atof(optarg); break;
    			case 's': seedlen = atoi(optarg); break;
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
            len_1 = len_pair = readinput(inputfile1, name1, seq_1, probDbl_1);
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

        /* maximal sequence shift in alignment is by default 50% of longer sequence */
        if( maxshift == INFINITE )
        	maxshift = abs(max(len_1, len_2) / 2);

        /* get unpaired probabilities */
        getunpaired(probDbl_1, len_1, probSgl_1);
        getunpaired(probDbl_2, len_2, probSgl_2);

        /* log-odds scores */
        getlogoddsSgl( probSgl_1, seq_1, len_1, pnull );
        getlogoddsSgl( probSgl_2, seq_2, len_2, pnull );
        getlogoddsDblNeighborhood( probDbl_1, seq_1, len_1, pnull, radius, theta );
        getlogoddsDblNeighborhood( probDbl_2, seq_2, len_2, pnull, radius, theta );

        /* reduce depth of probability matrices */
        reducematrix(probDbl_1, len_1*len_1, precision);
        reducematrix(probDbl_2, len_2*len_2, precision);
        reducematrix(probSgl_1, len_1, precision);
        reducematrix(probSgl_2, len_2, precision);
		#if DEBUG
        	for( int i=0; i<len_1; i++ ) { for( int j=0; j<len_1; j++ ) cerr << probDbl_1.at( i*len_1+j ) << "\t"; cerr << endl; }; cerr << endl;
        	for( int i=0; i<len_2; i++ ) { for( int j=0; j<len_2; j++ ) cerr << probDbl_2.at( i*len_2+j ) << "\t"; cerr << endl; }; cerr << endl;
        	for( int i=0; i<len_1; i++ ) cerr << probSgl_1[ i ] << " "; cerr << endl;
        	for( int i=0; i<len_2; i++ ) cerr << probSgl_2[ i ] << " "; cerr << endl;
		#endif

		/* semi-local (global) alignment */

        /* initialize arrays of indices */
		int *idx_1_aln = new int[len_1];
		for( int i=0; i<len_1; i++ ) idx_1_aln[ i ] = i;
		int *idx_2_aln = new int[len_2];
		for( int i=0; i<len_2; i++ ) idx_2_aln[ i ] = i;

		/* STEP 1: global alignment (Needleman-Wunsch algorithm) of pairing probabilities of each base in S_a and S_b */
		float **sim = new float*[ len_2 ];
		for( int j=0; j<len_2; j++ )
			sim[j] = new float[ len_1 ];
		for( int j=0; j<len_2; j++ )
			for( int i=0; i<len_1; i++ )
				sim[ j ][ i ] = INFINITE;

		/* vectors keeping temporary probabilities */
		float* subprobDbl_1 = new float[ len_1 ];
		float* subprobDbl_2 = new float[ len_2 ];
		float* subprobSgl_1 = new float[ len_1 ];
		float* subprobSgl_2 = new float[ len_2 ];

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
	    /* Initialize traceback matrices */
	    char ** trF = new char * [ len_2 + 1 ];
	    for( int j = 0; j <= len_2; j++ )
	    	trF[ j ] = new char [ len_1 + 1 ];
	    char ** trQ = new char * [ len_2 + 1 ];
	    for( int j = 0; j <= len_2; j++ )
	    	trQ[ j ] = new char [ len_1 + 1 ];
	    char ** trP = new char * [ len_2 + 1 ];
	    for( int j = 0; j <= len_2; j++ )
	    	trP[ j ] = new char [ len_1 + 1 ];

		/* run structure alignment of seed or local sequence alignment to position both dot plots */
		int offset = 0;
		if( setseqaln_flag )
			offset = nwseq_seed( seq_1, seq_2, F, trF );
		else
			if( seedlen && len_1 > 2*seedlen && len_2 > 2*seedlen )
				offset = nwdp_seed( seq_1, probDbl_1, probSgl_1, idx_1_aln, len_1, seq_2, probDbl_2, probSgl_2, idx_2_aln, len_2, seedlen, subprobDbl_1, subprobSgl_1, subprobDbl_2, subprobSgl_2, sim, F, Q, P, trF, trQ, trP, deltanull );

		/* run all pairs of pairing probabilities */
		int k, l;
		int *max = ( len_2 > len_1 ) ? &l : &k;
		int *min = ( len_2 <= len_1 ) ? &l : &k;
		for( l=0; l<len_2; l++)
			for( k = 0; k<len_1; k++)
				if( offset + *min + maxshift >= *max && offset + *min - maxshift <= *max )
					if( sim[ l ][ k ] == INFINITE )
						sim[ l ][ k ] = nwdp( seq_1, probDbl_1, probSgl_1, k, idx_1_aln, len_1, seq_2, probDbl_2, probSgl_2, l, idx_2_aln, len_2, subprobDbl_1, subprobSgl_1, subprobDbl_2, subprobSgl_2, 1, F, Q, P, trF, trQ, trP, deltanull );

		#if DEBUG
			cout << "Similarity matrix: " << endl;
			for( int j=0; j<len_2; j++) {
				for( int i=0; i<len_1; i++)
					cout << sim[ j ][ i ] << "\t";
				cout << endl;
			}
		#endif

		/* STEP 2: find best local (or global) path through similarity matrix */
		simalign_affinegaps( sim, len_1, len_2, idx_1_aln, idx_2_aln, len_pair, len_aln, setglobal2_flag, F, Q, P, trF, trQ, trP );

		#if DEBUG
			cout << "\tGaps_1 = " << len_1-len_pair << "\tGaps_2 = " << len_2-len_pair << endl;
			for( int i=0; i<len_pair; i++ ) cout << idx_1_aln[ i ] << "\t"; cout << endl;
			for( int j=0; j<len_pair; j++ ) cout << idx_2_aln[ j ] << "\t"; cout << endl;
			cout << "probDbl_1: " << len_pair << endl;
			for( int i=0; i<len_pair; i++) for( int j=0; j<len_pair; j++ ) cout << kappa * probSgl_1[ idx_1_aln[ j ] ] + ( 1 - kappa ) * probDbl_1[ idx_1_aln[ i ]*len_1 + idx_1_aln[ j ] ] << "\t"; cout << endl;
			cout << "probDbl_2: " << len_pair << endl;
			for( int i=0; i<len_pair; i++) for( int j=0; j<len_pair; j++ ) cout << kappa * probSgl_2[ idx_2_aln[ j ] ] + ( 1 - kappa ) * probDbl_2[ idx_2_aln[ i ]*len_2 + idx_2_aln[ j ] ] << "\t"; cout << endl;
		#endif

		/* calculate final similarity */
		float similarity = 0.;
		int open = 0, extended = 0;
		for( int i=0; i<len_pair; i++ )
			similarity += nwdp( seq_1, probDbl_1, probSgl_1, i, idx_1_aln, len_pair, seq_2, probDbl_2, probSgl_2, i, idx_2_aln, len_pair, subprobDbl_1, subprobSgl_1, subprobDbl_2, subprobSgl_2, 0, F, Q, P, trF, trQ, trP, 0.5 );
		affinegapcosts(idx_1_aln, idx_2_aln, len_pair, open, extended);
		//cout << similarity << " " << alpha*open << " " << beta*extended << " " << " " << len_pair << " " << len_aln << " " << extended << endl;
		similarity = ( len_pair ) ? ( similarity + alpha * open + beta * extended ) / len_aln + 0.5 : 0;

		/* OUTPUT */

		/* print aligned probabilities and similarity */
		cout << "Similarity = " << similarity << ", Length_1 = " << len_1 << ", Unaligned_1 = " << len_1-len_pair;
		cout << ", Length_2 = " << len_2 << ", Unaligned_2 = " << len_2-len_pair << endl;

		/* print sequences aligned by dot plot alignment */
		printalign(seq_1, idx_1_aln, seq_2, idx_2_aln, len_pair);

		/* free memory */
		freeMatrix(sim, len_2);
	    delete[] idx_1_aln;
	    delete[] idx_2_aln;
	    delete[] subprobDbl_1;
	    delete[] subprobDbl_2;
	    delete[] subprobSgl_1;
	    delete[] subprobSgl_2;
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
