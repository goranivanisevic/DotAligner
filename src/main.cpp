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

//RNA base statistics from Knudsen et al (1999) Bioinformatics, 15 (6), 446-454.
//unpaired probabilities    A         U         G         C
float pnullSgl[4] = {0.364097, 0.273013, 0.211881, 0.151009};
float maxpnullSgl = pnullSgl[0];

//unpaired rate matrix    A       U       G       C
float psubSgl[16] = {-0.749,  0.263,  0.322,  0.164,
				 	  0.351, -1.050,  0.185,  0.513,
				 	  0.553,  0.239, -0.964,  0.173,
					  0.396,  0.927,  0.242, -1.565};

//paired probabilities       A         U         G         C
float pnullDbl[16] = {0.001167, 0.177977, 0.001058, 0.001806,
					  0.177977, 0.002793, 0.049043, 0.000763,
					  0.001058, 0.049043, 0.000406, 0.266974,
					  0.001806, 0.000763, 0.266974, 0.000391};
float maxpnullDbl = pnullDbl[11];

//paired rate matrix      AA      AU      AG      AC      UA      UU      UG      UC      GA      GU      GG      GC      CA      CU      CG      CC
float psubDbl[256] = {-3.607,  0.617,  0.589,  0.420,  0.617,  0.000,  0.026,  0.019,  0.589,  0.026,  0.000,  0.132,  0.420,  0.019,  0.132,  0.000,
					   0.004, -1.163,  0.002,  0.029,  0.185,  0.016,  0.023,  0.003,  0.002,  0.274,  0.001,  0.501,  0.001,  0.007,  0.115,  0.002,
					   0.650,  0.257, -2.489,  0.116,  0.290,  0.000,  0.237,  0.019,  0.097,  0.024,  0.000,  0.060,  0.006,  0.000,  0.734,  0.000,
					   0.271,  2.861,  0.068, -6.070,  0.124,  0.020,  0.057,  0.008,  0.003,  0.401,  0.000,  2.135,  0.024,  0.010,  0.008,  0.079,
					   0.004,  0.185,  0.002,  0.001, -1.163,  0.016,  0.274,  0.007,  0.002,  0.023,  0.001,  0.115,  0.029,  0.003,  0.501,  0.002,
					   0.000,  0.988,  0.000,  0.013,  0.988, -4.261,  0.663,  0.129,  0.000,  0.663,  0.021,  0.287,  0.013,  0.129,  0.287,  0.079,
					   0.001,  0.084,  0.005,  0.002,  0.996,  0.038, -2.554,  0.005,  0.001,  0.037,  0.005,  0.101,  0.015,  0.001,  1.262,  0.002,
					   0.029,  0.688,  0.026,  0.018,  1.551,  0.473,  0.324, -6.126,  0.000,  0.044,  0.000,  1.265,  0.024,  1.105,  0.473,  0.105,
					   0.650,  0.290,  0.097,  0.006,  0.257,  0.000,  0.024,  0.000, -2.489,  0.237,  0.000,  0.734,  0.116,  0.019,  0.060,  0.000,
					   0.001,  0.996,  0.001,  0.015,  0.084,  0.038,  0.037,  0.001,  0.005, -2.554,  0.005,  1.262,  0.002,  0.005,  0.101,  0.002,
					   0.000,  0.252,  0.000,  0.000,  0.252,  0.145,  0.631,  0.000,  0.000,  0.631, -6.933,  2.511,  0.000,  0.000,  2.511,  0.000,
					   0.001,  0.334,  0.000,  0.014,  0.077,  0.003,  0.019,  0.004,  0.003,  0.232,  0.004, -0.823,  0.000,  0.001,  0.130,  0.001,
					   0.271,  0.124,  0.003,  0.024,  2.861,  0.020,  0.401,  0.010,  0.068,  0.057,  0.000,  0.008, -6.070,  0.008,  2.135,  0.079,
					   0.029,  1.551,  0.000,  0.024,  0.688,  0.473,  0.044,  1.105,  0.026,  0.324,  0.000,  0.473,  0.018, -6.126,  1.265,  0.105,
					   0.001,  0.077,  0.003,  0.000,  0.334,  0.003,  0.232,  0.001,  0.000,  0.019,  0.004,  0.130,  0.014,  0.004, -0.823,  0.001,
					   0.000,  0.799,  0.000,  0.365,  0.799,  0.563,  0.191,  0.204,  0.000,  0.191,  0.000,  0.992,  0.365,  0.204,  0.992, -5.666};
//psubDbl[ 16*(4*nucIdx.at(nuc_2)+nucIdx.at(seq_2.at(j-1))) + 4*nucIdx.at(nuc_1)+nucIdx.at(seq_1.at(i-1)) ];

map<const char, const int> nucIdx = { {'A', 0}, {'U', 1}, {'G', 2}, {'C', 3} };

// global arguments
float kappa = 0.5;
float alpha = -1;
float beta  = -0.1;

float INFINITE = -1000;


int main( int argc, char ** argv )
{
		char * program = *argv ;
        vector<float> probDbl_1;
        vector<float> probDbl_2;
        vector<float> probSgl_1;
        vector<float> probSgl_2;
        int len_1 = 0, len_2 = 0, len_aln = 0, len_1_last = 0, len_2_last = 0;
        string name1, name2, seq_1, seq_2;
        int precision = 4;
        int iter = -1;
        int maxshift = INFINITE;
        int seednr = 5;
        int seedlen = 5;
        float pnull = 0.0005;

        /* arguments */
        string  filename1, filename2;

        /* boolean flags */
        static int setprintmatrix_flag = 0;
        static int setseqalign_flag = 0;
        static int setlogodds_flag = 0;
        static int setlocal1_flag = 0;
        static int setglobal2_flag = 0;
        static int help_flag = 0;
        static int local = 0;

    	/* parsing long options */
    	while(1)
    	{
    		static struct option long_options[] =
    		{
    			{"printdp", no_argument, &setprintmatrix_flag, 1},	// Print matrices
    			{"seq", no_argument, &setseqalign_flag, 1},		// calculate sequence alignment
    			{"logodds", no_argument, &setlogodds_flag, 1},	// replace probabilities by expectation normalized log odds
    			{"local1", no_argument, &setlocal1_flag, 1},	// local alignment in step 1 (default is global in step 1)
    			{"global2", no_argument, &setglobal2_flag, 1},  // global alignment in step 2 (default is local alignment in step 2)
    			{"dotplot", required_argument, 0, 'd'},			// input dotplots
    			{"kappa", required_argument, 0, 'k'},			// weight of sequence similarity
    			{"alpha", required_argument, 0, 'a'},		 	// affine gap costs = alpha + k * beta (k gaps)
    			{"beta", required_argument, 0, 'b'},			// affine gap costs = alpha + k * beta (k gaps)
    			{"precision", required_argument, 0, 'p'},       // number of digits considered of base pair reliabilities
    			{"maxshift", required_argument, 0, 'm'},		// maximal shift of two input sequences in the final alignment (default 50% of longer sequence)
    			{"seednr", required_argument, 0, 's'},			// number of seed alignments (default 3)
    			{"seedlen", required_argument, 0, 'l'},			// length of seed alignments (default 5 nucleotides)
    			{"pnull", required_argument, 0, 'n'},			// minimal probability
    			{"help", no_argument, &help_flag, 1},
    			{0, 0, 0, 0}
    		};

    		/* print help text */
    		if( help_flag )
    			usage_dotaligner(program);

    		/* getopt_long stores the option index here. */
    		int option_index = 0;

    		int cmd = getopt_long(argc, argv, "d:k:a:b:p:m:s:l:n:", long_options, &option_index);

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
    			case 's': seednr = atoi(optarg); break;
    			case 'l': seedlen = atoi(optarg); break;
    			case 'n': pnull = atof(optarg); break;
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
            len_1 = len_aln = readinput(inputfile1, name1, seq_1, probDbl_1);
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

        /* run global sequence alignment */
        if( setseqalign_flag )
        {
            /* aligned sequences */
            string seq_1_al;
            seq_1_al.reserve(2*len_1);
            string seq_2_al;
            seq_2_al.reserve(2*len_2);

            /* get sequence alignment */
            float similarity = (float) nw( seq_1, seq_2, seq_1_al, seq_2_al, setprintmatrix_flag ) ;

            /* print sequence alignment */
            cout << "Similarity = " << similarity << endl;
            print_al( seq_1_al, seq_2_al );

            return 0;
        }

        /* maximal sequence shift in alignment is by default 50% of longer sequence */
        if( maxshift == INFINITE )
        	maxshift = abs(max(len_1, len_2) / 2);

        /* get unpaired probabilities */
        getunpaired(probDbl_1, len_1, probSgl_1);
        getunpaired(probDbl_2, len_2, probSgl_2);

        /* log-odds scores */
        getlogoddsDbl(probDbl_1, seq_1, len_1, pnull );
        getlogoddsDbl(probDbl_2, seq_2, len_2, pnull );
        getlogoddsSgl( probSgl_1, seq_1, len_1, pnull );
        getlogoddsSgl( probSgl_2, seq_2, len_2, pnull );

        /* reduce depth of probability matrices */
        reducematrix(probDbl_1, len_1*len_1, precision);
        reducematrix(probDbl_2, len_2*len_2, precision);
        reducematrix(probSgl_1, len_1, precision);
        reducematrix(probSgl_2, len_2, precision);

		/* local alignment */
		if( setlocal1_flag )
		{
			int nrlocal = 3;   // number of optimal local alignments that should be calculated
			float overlap = 0.5;   // maximal overlap allowed of local hit before they are concatenated
			float minLsim = 0.01;   // minimal local similarity to be considered
			vector<LocalHit> lhit;
			for( int i=0; i<nrlocal; i++ ) {
				lhit.push_back( LocalHit() );
				lhit[ i ].similarity = 0;
			}

	        /* initialize arrays of indices */
			int *idx_1_aln_local = new int[len_1];
			int *idx_2_aln_local = new int[len_2];

			/* local alignment for each pair of pairing probabilities that includes the two reference bases */
			for( int l=0; l<len_2; l++)
				for( int k = 0; k<len_1; k++)
				{
					for( int i=0; i<len_1; i++ ) idx_1_aln_local[ i ] = i;
					for( int i=0; i<len_2; i++ ) idx_2_aln_local[ i ] = i;

					/* locally align pairing probabilities of S_a(k) and S_b(l) */
					float similarity;
					int len_local = 1;
					similarity = nwdp( seq_1, probDbl_1, probSgl_1, idx_1_aln_local[k], idx_1_aln_local, len_1, seq_2, probDbl_2, probSgl_2, idx_2_aln_local[l], idx_2_aln_local, len_2, len_local, setprintmatrix_flag );
					cerr << "Local similarity for pair " << k << " and " << l << " -> " << similarity << " from " << idx_1_aln_local[ 0 ] << " to " << idx_1_aln_local[ len_local-1 ] << " in S_a and from " << idx_2_aln_local[ 0 ] << " to " << idx_2_aln_local[ len_local-1 ] << " in S_b of length " << len_local << endl;
					//cerr << "Local Sim = " << similarity << "; Length = " << len_local;
					//cerr << "; LS_a START = " << idx_1_aln_local[ 0 ] << "; LS_a END = " << idx_1_aln_local[ len_local-1 ];
					//cerr << "; LS_b START = " << idx_2_aln_local[ 0 ] << "; LS_b END = " << idx_2_aln_local[ len_local-1 ] << endl;

					/* for optimal local alignment we test if the aligned local sequences LS_a and LS_b include the base k and base l
					 * and if no more than <nrlocal-1> local alignments with higher similarity exist */
					if( similarity > lhit.front().similarity && k >= idx_1_aln_local[ 0 ] && k <= idx_1_aln_local[ len_local-1 ] && l >= idx_2_aln_local[ 0 ] && l <= idx_2_aln_local[ len_local-1 ] )
					{
						/* if local hit overlaps another hit more then <overlap> in S_a and S_b then concatenate */
						vector<LocalHit>::iterator lhitIt;
						for( lhitIt = lhit.begin(); lhitIt != lhit.end(); lhitIt++ )
						{
							if( ! lhitIt.base()->similarity ) continue;

							float over1 = getOverlap2ndInterval( idx_1_aln_local[ 0 ], idx_1_aln_local[ len_local-1 ], lhitIt.base()->lstart_1, lhitIt.base()->lend_1 );
							float over2 = getOverlap2ndInterval( idx_2_aln_local[ 0 ], idx_2_aln_local[ len_local-1 ], lhitIt.base()->lstart_2, lhitIt.base()->lend_2 );
							float over1r = getOverlap2ndInterval( lhitIt.base()->lstart_1, lhitIt.base()->lend_1, idx_1_aln_local[ 0 ], idx_1_aln_local[ len_local-1 ] );
							float over2r = getOverlap2ndInterval( lhitIt.base()->lstart_2, lhitIt.base()->lend_2, idx_2_aln_local[ 0 ], idx_2_aln_local[ len_local-1 ] );
							if( over1 > overlap || over1r > overlap )
								if( over2 > overlap || over2r > overlap )
								{
									// concatenate
									lhitIt.base()->lstart_1 = min( lhitIt.base()->lstart_1, idx_1_aln_local[ 0 ] );
									lhitIt.base()->lend_1 = max( lhitIt.base()->lend_1, idx_1_aln_local[ len_local-1 ] );
									lhitIt.base()->lstart_2 = min( lhitIt.base()->lstart_2, idx_2_aln_local[ 0 ] );
									lhitIt.base()->lend_2 = max( lhitIt.base()->lend_2, idx_2_aln_local[ len_local-1 ] );
									if( similarity > lhitIt.base()->similarity )
										lhitIt.base()->similarity = similarity;

									break;
								}
						}
						if( lhitIt != lhit.end() ) continue;

						/* store similarity, start and end position of LS_a and LS_b */
						lhit.front().similarity = similarity;
						lhit.front().lstart_1 = idx_1_aln_local[ 0 ];
						lhit.front().lend_1 = idx_1_aln_local[ len_local-1 ];
						lhit.front().lstart_2 = idx_2_aln_local[ 0 ];
						lhit.front().lend_2 = idx_2_aln_local[ len_local-1 ];

						/* sort top <nrlocal> local hits with lowest similarities */
						sort( lhit.begin(), lhit.end(), compareBySimilarity );
					}
				}

			/* run global DotAligner for top <nrlocal> LS_a against S_b */
			for( vector<LocalHit>::iterator lhitIt = lhit.end()-1; lhitIt != lhit.begin()-1; lhitIt-- )
			{
				if( lhitIt.base()->similarity < minLsim )
					continue;
				cerr << "LocalHit: " << lhitIt.base()->similarity << " " << lhitIt.base()->lstart_1 << " " << lhitIt.base()->lend_1 << " " << lhitIt.base()->lstart_2 << " " << lhitIt.base()->lend_2 << endl;

				/* run STEP 1 with global alignments of LS_a against  LS_b */
				int len_local_1 = lhitIt.base()->lend_1 - lhitIt.base()->lstart_1 + 1;
				for( int i=0; i < len_local_1; i++ ) idx_1_aln_local[ i ] = lhitIt.base()->lstart_1 + i;
				int len_local_2 = lhitIt.base()->lend_2 - lhitIt.base()->lstart_2 + 1;
				for( int i=0; i < len_local_2; i++ ) idx_2_aln_local[ i ] = lhitIt.base()->lstart_2 + i;
				float **sim_local = new float*[ len_local_2 ];
				for( int n=0; n<len_local_2; n++ ) {
					sim_local[n] = new float[ len_local_1 ];
					for( int m = 0; m<len_local_1; m++ )
						sim_local[ n ][ m ] = nwdp( seq_1, probDbl_1, probSgl_1, idx_1_aln_local[m], idx_1_aln_local, len_local_1, seq_2, probDbl_2, probSgl_2, idx_2_aln_local[n], idx_2_aln_local, len_local_2, local, setprintmatrix_flag );
				}

				/* run STEP 2 on the similarities of LS_a */
				int *tmp_idx_1_aln_local = new int[len_local_1];
				int *tmp_idx_2_aln_local = new int[len_local_2];
				//for( int i=0; i<len_local; i++ ) cout << idx_1_aln_local[ i ] << ","; cout << endl;
				//for( int j=0; j<len_local; j++ ) cout << idx_2_aln_local[ j ] << ","; cout << endl;
				simalign_affinegaps( sim_local, len_local_1, len_local_2, tmp_idx_1_aln_local, tmp_idx_2_aln_local, len_aln, precision, 0, setprintmatrix_flag);
				for( int i=0; i<len_aln; i++ ) {
					idx_1_aln_local[ i ] = idx_1_aln_local[ tmp_idx_1_aln_local[ i ] ];
					idx_2_aln_local[ i ] = idx_2_aln_local[ tmp_idx_2_aln_local[ i ] ];
				}
				delete[] tmp_idx_1_aln_local;
				delete[] tmp_idx_2_aln_local;

				/* calculate final similarity */
				float similarity = 0.;
				for( int i=0; i<len_aln; i++ )
					similarity += nwdp( seq_1, probDbl_1, probSgl_1, idx_1_aln_local[i], idx_1_aln_local, len_aln, seq_2, probDbl_2, probSgl_2, idx_2_aln_local[i], idx_2_aln_local, len_aln, local, setprintmatrix_flag );
				int open = 0, extended = 0;
				affinegapcosts(idx_1_aln_local, idx_2_aln_local, len_aln, open, extended);
				//similarity = ( len_1_aln ) ? ( similarity + alpha * ( len_1 + len_2 - len_1_aln - len_1_aln ) ) / len_1_aln : 0;
				similarity = ( len_aln ) ? ( similarity + alpha * open + beta * extended ) / ( len_aln + extended ) : 0;

				/* OUTPUT */
				/* print aligned probabilities and similarity */
				cout << "Similarity = " << similarity << ", Length_1 = " << len_1 << ", Unaligned_1 = " << len_1-len_aln;
				cout << ", Length_2 = " << len_2 << ", Unaligned_2 = " << len_2-len_aln << endl;

				/* print sequences aligned by dot plot alignment */
				for( int i=0; i<len_aln; i++ ) cout << idx_1_aln_local[ i ] << ","; cout << endl;
				for( int j=0; j<len_aln; j++ ) cout << idx_2_aln_local[ j ] << ","; cout << endl;
				printalign(seq_1, idx_1_aln_local, seq_2, idx_2_aln_local, len_aln);

				/* free memory */
				freeMatrix(sim_local, len_2);
			}

			/* free memory */
	        delete[] idx_1_aln_local;
	        delete[] idx_2_aln_local;
		}

		/* semi-local (global) alignment */
		else
		{
	        /* initialize arrays of indices */
			int *idx_1_aln = new int[len_1];
			for( int i=0; i<len_1; i++ ) idx_1_aln[ i ] = i;
			int *idx_2_aln = new int[len_2];
			for( int i=0; i<len_2; i++ ) idx_2_aln[ i ] = i;

			/* repeat alignment until no changes in alignment */
			//while( len_1_last != len_aln || len_2_last != len_aln )
			//{
				len_1_last = ( iter == -1 ) ? len_1 : len_aln;
				len_2_last = ( iter == -1 ) ? len_2 : len_aln;

				/* STEP 1: global alignment (Needleman-Wunsch algorithm) of pairing probabilities of each base in S_a and S_b */
				float **sim = new float*[ len_2_last ];
				for( int j=0; j<len_2_last; j++ )
					sim[j] = new float[ len_1_last ];
				for( int j=0; j<len_2_last; j++ )
					for( int i=0; i<len_1_last; i++ )
						sim[ j ][ i ] = INFINITE;

				/* run seed alignments */
				if( iter == -1 && seednr && seedlen && len_1_last > 2*seednr*seedlen && len_2_last > 2*seednr*seedlen ) {
					// structure similarity per base in seed alignment necessary for success
					float threshold = 0.5;
					int *idx_seed_aln = new int[seedlen];
					int len_seed_aln;
					srand((unsigned)time(0));
					int maxlen = ( len_2_last > len_1_last ) ? len_2_last : len_1_last;
					int minlen = ( len_2_last <= len_1_last ) ? len_2_last : len_1_last;
					float **seedsim = new float*[ maxlen ];
					for( int j=0; j<maxlen; j++ )
						seedsim[j] = new float[ seedlen ];
					for( int s=0; s<seednr; s++ )
					{
						for( int j=0; j<maxlen; j++ )
							for( int i=0; i<seedlen; i++ )
								seedsim[ j ][ i ] = INFINITE;
						cerr << "SEED ALIGNMENT NR " << s+1 << endl;
						// get seed sequence from longer sequence
						int seedstart = int( (float)(minlen-seedlen)*rand()/(RAND_MAX + 1.0) );
						//seedstart =+ int( (rand()%minlen)+1 / 2 );
						for( int i=0; i<seedlen; i++ )
							idx_seed_aln[ i ] = seedstart + i;
						// run nwdp
						float tmpkappa = kappa; kappa = 0;
						float tmpalpha = alpha; alpha = 0;
						float tmpbeta  = beta; beta = 0;
						for( int l=0; l<maxlen; l++)
							for( int k = seedstart; k<seedstart+seedlen; k++)
								if( k + maxshift >= l && k - maxshift <= l )
								{
									if( len_2_last > len_1_last )
										sim[ l ][ k ] = seedsim[ l ][ k-seedstart ] = nwdp( seq_1, probDbl_1, probSgl_1, idx_1_aln[k], idx_1_aln, len_1_last, seq_2, probDbl_2, probSgl_2, idx_2_aln[l], idx_2_aln, len_2_last, local, 0 );
									else
										sim[ k ][ l ] = seedsim[ l ][ k-seedstart ] = nwdp( seq_1, probDbl_1, probSgl_1, idx_1_aln[l], idx_1_aln, len_1_last, seq_2, probDbl_2, probSgl_2, idx_2_aln[k], idx_2_aln, len_2_last, local, 0 );
								}
						// run simalign_affinegaps
						int *tmp_idx_seed_aln = new int[seedlen];
						int *tmp_idx_max_aln = new int[maxlen];
						beta = INFINITE;
						float z = simalign_affinegaps( seedsim, seedlen, maxlen, tmp_idx_seed_aln, tmp_idx_max_aln, len_seed_aln, precision, 0, 0);
						kappa = tmpkappa;
						alpha = tmpalpha;
						beta = tmpbeta;
						cerr << "Z " << z << " SeedLen " << len_seed_aln << endl;
						for( int i=0; i<len_seed_aln; i++ ) tmp_idx_seed_aln[ i ] = tmp_idx_seed_aln[ i ] + seedstart;
						if( len_2_last > len_1_last ) {
							for( int i=0; i<len_seed_aln; i++ ) cout << tmp_idx_seed_aln[ i ] << " " << seq_1.at(tmp_idx_seed_aln[ i ]) << "\t"; cout << endl;
							for( int j=0; j<len_seed_aln; j++ ) cout << tmp_idx_max_aln[ j ] << " " << seq_2.at(tmp_idx_max_aln[ j ]) << "\t"; cout << endl;
						}
						else {
							for( int j=0; j<len_seed_aln; j++ ) cout << tmp_idx_max_aln[ j ] << " " << seq_1.at(tmp_idx_max_aln[ j ]) << "\t"; cout << endl;
							for( int i=0; i<len_seed_aln; i++ ) cout << tmp_idx_seed_aln[ i ] << " " << seq_2.at(tmp_idx_seed_aln[ i ]) << "\t"; cout << endl;
						}
						delete[] tmp_idx_seed_aln;
						int targetlen = tmp_idx_max_aln[ len_seed_aln-1 ] - tmp_idx_max_aln[ 0 ] + 1;
						delete[] tmp_idx_max_aln;

						// break if sim >= threshold and gapnr == 0 and alignlength == seedlength
						cerr << "THRESHOLD " << threshold * seedlen << " TARGETLEN " << targetlen << endl;
						if( len_seed_aln == seedlen && targetlen == seedlen && z >= threshold * seedlen )
							break;
						else
							if( s == seednr - 1 ) {
								// exit with similarity = 0 if sim < threshold or gapnr != 0 or alignlength != seedlength
								cout << "Similarity = 0, Length_1 = 0, Unaligned_1 = " << len_1;
								cout << ", Length_2 = 0, Unaligned_2 = " << len_2 << ", Iterations = " << iter << endl;

								/* free memory */
								delete[] idx_1_aln;
								delete[] idx_2_aln;
								probDbl_1.clear();
								probDbl_2.clear();
								probSgl_1.clear();
								probSgl_2.clear();
								freeMatrix(seedsim, maxlen);

								return 0;
							}
					}
					freeMatrix(seedsim, maxlen);
				}

				/* run all pairs of pairing probabilities */
				int k, l;
				int *max = ( len_2_last > len_1_last ) ? &l : &k;
				int *min = ( len_2_last <= len_1_last ) ? &l : &k;
				int count=0;
				for( l=0; l<len_2_last; l++)
					for( k = 0; k<len_1_last; k++)
						if( *min + maxshift >= *max && *min - maxshift <= *max ) {
							if( sim[ l ][ k ] == INFINITE )
								sim[ l ][ k ] = nwdp( seq_1, probDbl_1, probSgl_1, idx_1_aln[k], idx_1_aln, len_1_last, seq_2, probDbl_2, probSgl_2, idx_2_aln[l], idx_2_aln, len_2_last, local, setprintmatrix_flag );
							count++;
						}
				cerr << "ENDNWDP " << maxshift << " " << count << endl;

				if( setprintmatrix_flag )
				{
					cout << "Similarity matrix: " << endl;
					for( int j=0; j<len_2_last; j++) {
						for( int i=0; i<len_1_last; i++)
							cout << sim[ j ][ i ] << "\t";
						cout << endl;
					}
				}

				/* STEP 2: find best local path through similarity matrix */
				int *tmp_idx_1_aln = new int[len_1_last];
				int *tmp_idx_2_aln = new int[len_2_last];
				simalign_affinegaps( sim, len_1_last, len_2_last, tmp_idx_1_aln, tmp_idx_2_aln, len_aln, precision, setglobal2_flag, setprintmatrix_flag);
				for( int i=0; i<len_aln; i++ ) {
					idx_1_aln[ i ] = idx_1_aln[ tmp_idx_1_aln[ i ] ];
					idx_2_aln[ i ] = idx_2_aln[ tmp_idx_2_aln[ i ] ];
				}
				delete[] tmp_idx_1_aln;
				delete[] tmp_idx_2_aln;

				#if DEBUG
					cout << "\tGaps_1 = " << len_1-len_1_aln << "\tGaps_2 = " << len_2-len_2_aln << endl;
					for( int i=0; i<len_1_aln; i++ ) cout << idx_1_aln[ i ] << "\t";
					cout << endl;
					for( int j=0; j<len_2_aln; j++ ) cout << idx_2_aln[ j ] << "\t";
					cout << endl;

					cout << "probDbl_1: " << len_1_aln << endl;
					for( int i=0; i<len_1_aln; i++) for( int j=0; j<len_1_aln; j++ ) cout << probDbl_1[ idx_1_aln[ i ]*len_1 + idx_1_aln[ j ] ] << "\t"; cout << endl;
					cout << "probDbl_2: " << len_2_aln << endl;
					for( int i=0; i<len_2_aln; i++) for( int j=0; j<len_2_aln; j++ ) cout << probDbl_2[ idx_2_aln[ i ]*len_2 + idx_2_aln[ j ] ] << "\t"; cout << endl;
				#endif

				iter++;

				/* free memory */
				freeMatrix(sim, len_2_last);
			//}

			/* calculate final similarity */
			float similarity = 0.;
			for( int i=0; i<len_aln; i++ )
				similarity += nwdp( seq_1, probDbl_1, probSgl_1, idx_1_aln[i], idx_1_aln, len_aln, seq_2, probDbl_2, probSgl_2, idx_2_aln[i], idx_2_aln, len_aln, local, setprintmatrix_flag );
			int open = 0, extended = 0;
			affinegapcosts(idx_1_aln, idx_2_aln, len_aln, open, extended);
			//similarity = ( len_1_aln ) ? ( similarity + alpha * ( len_1 + len_2 - len_1_aln - len_1_aln ) ) / len_1_aln : 0;
			similarity = ( len_aln ) ? ( similarity + alpha * open + beta * extended ) / ( len_aln + extended ) : 0;

			/* OUTPUT */

			/* print aligned probabilities and similarity */
			iter = ( !iter ) ? 1 : iter;
			cout << "Similarity = " << similarity << ", Length_1 = " << len_1 << ", Unaligned_1 = " << len_1-len_aln;
			cout << ", Length_2 = " << len_2 << ", Unaligned_2 = " << len_2-len_aln << ", Iterations = " << iter << endl;

			/* print sequences aligned by dot plot alignment */
			for( int i=0; i<len_aln; i++ ) cout << idx_1_aln[ i ] << ","; cout << endl;
			for( int j=0; j<len_aln; j++ ) cout << idx_2_aln[ j ] << ","; cout << endl;
			printalign(seq_1, idx_1_aln, seq_2, idx_2_aln, len_aln);

			/* free memory */
	        delete[] idx_1_aln;
	        delete[] idx_2_aln;
		}

		/* free memory */
        probDbl_1.clear();
        probDbl_2.clear();
        probSgl_1.clear();
        probSgl_2.clear();

        return  0;
}
