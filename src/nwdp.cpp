/*
 * nwdp.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Stefan Seemann, seemann@rth.dk
 */

#include "nwdp.h"

using namespace std;


void usage_dotaligner(char *program)
{
	cerr << "\nDotAligner v0.2";
	cerr << "\n===============";
	cerr << "\n   by Stefan E Seemann (seemann@rth.dk)\n";
	cerr << "\n   Usage:   " << program << " -d <file> -d <file> [ options ]";
    cerr << "\n   Semi-local pairwise alignment of two base pair probability matrices (dotplots).\n";
    cerr << "\n   -d --dotplot <file>    ... dotplot file (matrix of base pair probabilities)";
    cerr << "\n   Debug modes:";
	cerr << "\n   --printdp              ... print dynamic programming matrices";
    cerr << "\n   --global2              ... global alignment in step 2 (default is local alignment in step 2)";
    cerr << "\n   -k --kappa <float>     ... weight of sequence similarity (0..1; default = 0.5); 1 - kappa is weight of dotplots";
    cerr << "\n   -a --alpha <float>     ... affine gap costs = alpha + k * beta (k gaps; alpha default = -4)";
    cerr << "\n                              set <alpha> equal 0 to make the score becomes larger as a linear function of gap length";
    cerr << "\n   -b --beta <float>      ... affine gap costs = alpha + k * beta (k gaps; beta default = -1)";
    cerr << "\n                              set <beta> equal 0 to keep the score similar regardless of gap length";
    cerr << "\n   -p --precision <int>   ... number of digits considered of log-odds of base pair reliabilities (default = 4)";
    cerr << "\n   -m --maxshift <int>    ... maximal shift of two input sequences in the final alignment (default 50% of longer sequence)";
    cerr << "\n                              speeds up the calculation by ignoring pairwise comparisons of distant nucleotides";
    cerr << "\n                              local alignments may be missed for long sequences if <maxshift> is set too low";
    cerr << "\n   -n --pnull <prob>      ... minimal probability (default = 0.0005 )";
	cerr << "\n   -r --radius <int>      ... consider diagonal (stem) neighborhood of base pair probabilities of radius (default = 0)";
	cerr << "\n   -t --theta <float>     ... weight of neighborhood of radius for base pair probability (default = 0.25)";
	cerr << "\n   -l --deltanull <float> ... value of tau if both compared probabilities are zero (default = 0.5)";
	cerr << "\n   -s --seedlen <int>     ... run structure alignment of seed of specified length at the 5' end of shorter sequence";
	cerr << "\n                              to position both sequences (default = 15)";
	cerr << "\n   --seqaln               ... run local sequence alignment to position both sequences; <seqaln> is preferential to <seedlen>";
	cerr << "\n   --help                 ... this output\n";
    cerr << "\n   Output:   Similarity and alignment\n\n";

    exit( 1 ) ;
}


int readinput ( istream & is, string & name, string & seq, vector<float> & prob )
{
	// read probability matrix
	string line;
	int lnr = 0;
	int i = -1;

	while( !is.eof() ) {
		getline(is, line);

		if( i == -1 )
			name = line.replace(0,1,"");
		else
			if( i == 0 ) {
				seq = line;
				replace( seq.begin(), seq.end(), 'T', 'U' );
				lnr = line.length();
				prob.reserve( lnr * lnr );
			}
			else {
				vector<string> probv = split(line, ' ');
				for( unsigned int k = 0; k < probv.size(); k++ )
					prob.push_back( atof( probv[k].c_str() ) );
			}

		if( i == lnr )
			break;
		i++;
	}

	return lnr;
}


void getunpaired( vector<float> & prob, int len, vector<float> & probSgl )
{
	float psingle;

	probSgl.reserve( len );

	for( int i = 0; i < len; i++ )
	{
		psingle = 1;
		for( int j = 0; j < len; j++ )
			psingle -= prob.at( i*len + j );
		probSgl.push_back( psingle );
	}
}


vector<string> &split(const string &s, char delim, vector<string> &elems)
{
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s, char delim)
{
    vector<string> elems;
    return split(s, delim, elems);
}


void reducematrix(vector<float> & prob, int len, int precision)
{
	 char p[100];
     for( int i = 0; i < len; i++ ) {
    	 	sprintf(p, "%.*lf", precision, prob[i]);
    	 	prob.at(i) = atof(p);
     }
}


/**
   alignment of the pairing probabilities of two reference nucleotides (one line per dot plot)

   \param[seq_1]  sequence 1 ( S_a )
   \param[probDbl_1]  pairing probabilities of S_a (dot plot)
   \param[probSgl_1]  unpaired probabilities of S_a
   \param[rel_idx_1]  index of reference nucleotide of S_a in vector idx_1_aln
   \param[idx_1_aln]  indices of nucleotides in the subsequence of S_a that is considered for alignment ( SS_a )
   \param[len_1] length of considered subsequence SS_a
   \param[seq_2]  sequence 2 ( S_b )
   \param[probDbl_2]  pairing probabilities of S_b (dot plot)
   \param[probSgl_2]  unpaired probabilities of S_b
   \param[rel_idx_2]  index of reference nucleotide of S_b in vector idx_2_aln
   \param[idx_2_aln]  indices of nucleotides in the subsequence of S_b that is considered for alignment ( SS_b )
   \param[len_2] length of considered subsequence SS_b

   \returns  Similarity score ( 0 .. 1 )
*/
float nwdp( string seq_1, vector<float> & probDbl_1, vector<float> & probSgl_1, int rel_idx_1, int * idx_1_aln, int len_1, string seq_2, vector<float> & probDbl_2, vector<float> & probSgl_2, int rel_idx_2, int * idx_2_aln, int len_2, float * subprobDbl_1, float * subprobSgl_1, float * subprobDbl_2, float * subprobSgl_2, bool freeEndGaps, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP )
{
    int idx_1 = idx_1_aln[ rel_idx_1 ];
    int idx_2 = idx_2_aln[ rel_idx_2 ];

    int sidx_1 = idx_1 * seq_1.length();
    int sidx_2 = idx_2 * seq_2.length();

    /* create global alignment */
   	float tau = 0;

   	/* alignment upstream of the constraint pair ( rel_idx_1, rel_idx_2 ) */
   	if( rel_idx_1 > 0 && rel_idx_2 > 0 )
   	{
    	for( int i = 0; i < rel_idx_1; i++ ) {
    		subprobDbl_1[ i ] = probDbl_1.at( sidx_1 + idx_1_aln[ rel_idx_1 - 1 - i ] );
    		subprobSgl_1[ i ] = probSgl_1.at( idx_1_aln[ rel_idx_1 - 1 - i ] );
    	}
   		for( int j = 0; j < rel_idx_2; j++ ) {
   			subprobDbl_2[ j ] = probDbl_2.at( sidx_2 + idx_2_aln[ rel_idx_2 - 1 - j ] );
    		subprobSgl_2[ j ] = probSgl_2.at( idx_2_aln[ rel_idx_2 - 1 - j ] );
    	}
   		tau = nwdb_global_align_affinegaps( subprobDbl_1, subprobSgl_1, rel_idx_1, subprobDbl_2, subprobSgl_2, rel_idx_2, freeEndGaps, F, Q, P, trF, trQ, trP );
    }
    /* alignment downstream of the constraint pair ( rel_idx_1, rel_idx_2 ) */
    if( rel_idx_1<len_1-1 && rel_idx_2<len_2-1 )
    {
		for( int i = 0; i < len_1 - rel_idx_1 - 1; i++ ) {
			subprobDbl_1[ i ] = probDbl_1.at( sidx_1 + idx_1_aln[ rel_idx_1 + i + 1 ] );
			subprobSgl_1[ i ] = probSgl_1.at( idx_1_aln[ rel_idx_1 + i + 1 ] );
		}
		for( int j = 0; j < len_2 - rel_idx_2 - 1; j++ ) {
			subprobDbl_2[ j ] = probDbl_2.at( sidx_2 + idx_2_aln[ rel_idx_2 + j + 1 ] );
    		subprobSgl_2[ j ] = probSgl_2.at( idx_2_aln[ rel_idx_2 + j +1 ] );
		}
		tau += nwdb_global_align_affinegaps( subprobDbl_1, subprobSgl_1, len_1 - rel_idx_1 - 1, subprobDbl_2, subprobSgl_2, len_2 - rel_idx_2 - 1, freeEndGaps, F, Q, P, trF, trQ, trP );
    }

    /* sequence similarity */
   	float sigma = nwdb_align_seq_sim( seq_1.c_str()[idx_1], probSgl_1.at(idx_1), seq_2.c_str()[idx_2], probSgl_2.at(idx_2) );

   	/* final similarity */
    int minlen = ( len_2 > len_1 ) ? len_1 : len_2;
   	//cerr << idx_1 << " " << idx_2 << " " << sigma << " " << tau << " " << minlen << " " << ( tau / minlen ) << " " << kappa * sigma + ( 1 - kappa ) * tau / ( minlen - 1 ) << endl;
   	return kappa * sigma + ( 1 - kappa ) * tau / ( minlen - 1 );
    //return tau / ( minlen - 1 );
}


/*
 * align seed of shorter sequence against longer sequence and return starting position of seed alignment
 */
int nwdp_seed( string seq_1, vector<float> & probDbl_1, vector<float> & probSgl_1, int * idx_1_aln, int len_1, string seq_2, vector<float> & probDbl_2, vector<float> & probSgl_2, int * idx_2_aln, int len_2, int seedlen, float * subprobDbl_1, float * subprobSgl_1, float * subprobDbl_2, float * subprobSgl_2, float ** sim, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP )
{
	int maxlen = ( len_2 > len_1 ) ? len_2 : len_1;
	int minlen = ( len_2 > len_1 ) ? len_1 : len_2;

	float **seedsim = new float*[ maxlen ];
	for( int j=0; j<maxlen; j++ )
		seedsim[j] = new float[ minlen ];
	for( int j=0; j<maxlen; j++ )
    	for( int i=0; i<minlen; i++ )
    		seedsim[ j ][ i ] = INFINITE;

    /* get seed sequence from shorter sequence */
	for( int k = 0; k < seedlen; k++ )
		for( int l = 0; l < maxlen - seedlen - k; l++ )
		{
			if( len_2 > len_1 )
				sim[ l ][ k ] = seedsim[ l ][ k ] = nwdp( seq_1, probDbl_1, probSgl_1, k, idx_1_aln, len_1, seq_2, probDbl_2, probSgl_2, l, idx_2_aln, len_2, subprobDbl_1, subprobSgl_1, subprobDbl_2, subprobSgl_2, 1, F, Q, P, trF, trQ, trP );
			else
				sim[ k ][ l ] = seedsim[ l ][ k ] = nwdp( seq_1, probDbl_1, probSgl_1, l, idx_1_aln, len_1, seq_2, probDbl_2, probSgl_2, k, idx_2_aln, len_2, subprobDbl_1, subprobSgl_1, subprobDbl_2, subprobSgl_2, 1, F, Q, P, trF, trQ, trP );
		}

	/* run simalign_affinegaps */
	int *tmp_idx_seed_aln = new int[seedlen];
	int *tmp_idx_max_aln = new int[maxlen];
	int len_seed_pair, len_seed_aln;
	if( len_2 > len_1 )
		simalign_affinegaps( seedsim, seedlen, maxlen, tmp_idx_seed_aln, tmp_idx_max_aln, len_seed_pair, len_seed_aln, 0, F, Q, P, trF, trQ, trP );
	else
		simalign_affinegaps( seedsim, maxlen, seedlen, tmp_idx_max_aln, tmp_idx_seed_aln, len_seed_pair, len_seed_aln, 0, F, Q, P, trF, trQ, trP );
	int offset = tmp_idx_max_aln[ 0 ];

    /* free memory */
	delete[] tmp_idx_seed_aln;
    delete[] tmp_idx_max_aln;
    freeMatrix(seedsim, maxlen);

    return offset;
}


int nwseq_seed( string seq_1, string seq_2, float ** F, char ** trF )
{
	int x, y;
	float fU, fD, fL;
	char ptr, nuc ;
	int i, j;

	const int a = 2;   // Match like in blastn (optimized for sequence identity of 90%)
	const int b = -3;   // Mismatch like in blastn
	const int d = -5;   // gap penalty like in blastn existence gap costs

	const int  s[ 4 ][ 4 ] = { { a, b, b, b },    /* substitution matrix */
	                           { b, a, b, b },
	                           { b, b, a, b },
	                           { b, b, b, a } } ;

	int L1 = seq_1.length();
	int L2 = seq_2.length();

    /* Initialize dynamic programming matrix */
    nwdp_initF_affinegaps( F, L1, L2, 1 );

    /* Initialize traceback matrix */
    nwdp_initTB( trF, L1, L2 );

	for( j = 1; j <= L2; j++ )
	{
		for( i = 1; i <= L1; i++ )
	    {
			nuc = seq_1[ i-1 ] ;

	        switch( nuc )
	        {
	        	case 'A':  x = 0 ;  break ;
	            case 'C':  x = 1 ;  break ;
	            case 'G':  x = 2 ;  break ;
	            case 'U':  x = 3 ;  break ;
	        }

	        nuc = seq_2[ j-1 ] ;

	        switch( nuc )
	        {
	        	case 'A':  y = 0 ;  break ;
	            case 'C':  y = 1 ;  break ;
	            case 'G':  y = 2 ;  break ;
	            case 'U':  y = 3 ;  break;
	        }

	        fU = F[ j-1 ][ i ] + d ;
	        fD = F[ j-1 ][ i-1 ] + s[ x ][ y ] ;
	        fL = F[ j ][ i-1 ] + d ;

	        F[ j ][ i ] = (int) max( 0, max3( (float) fU, (float) fD, (float) fL, &ptr ) );

	        trF[ j ][ i ] =  ptr ;
	    }
	}

	/* find maximal entry of similarity in matrix F for local alignments */
    float lmax = 0;
    int li = 0;
    int lj = 0;
    for( j = 1; j <= L2; j++ )
    	for( i = 1; i <= L1; i++ ) {
    		if( F[ j ][ i ] >= lmax ) {
    			lmax = F[ j ][ i ];
        		li = i;
        		lj = j;
    		}
    	}
    i = li;
    j = lj;

	while( i > 0 || j > 0 )
	{
		if( F[ j ][ i ] == 0 )
			break;

		switch( trF[ j ][ i ] )
	    {
	    	case '|' : j-- ;
	                   break ;
            case '\\': i-- ; j-- ;
	                   break ;
            case '-' : i-- ;
	                   break;
	    }
	}

	#if DEBUG
		float* tempidx_1 = new float[L1+1];
		float* tempidx_2 = new float[L2+1];
		for( int k = 0; k <= L1; k++ ) tempidx_1[ k ] = (float) k;
		for( int k = 0; k <= L2; k++ ) tempidx_2[ k ] = (float) k;
		cout << "F:" << endl;
		print_matrixdp( F, tempidx_1, L1, tempidx_2, L2);
		print_matrixdp( trF, tempidx_1, L1, tempidx_2, L2);
	#endif
	cerr << "SEQALN: Seq_1 ( " << i << ", " << li << " ); Seq_2 ( " << j << ", " << lj << " )" << endl;

	return abs( i - j );
}


/*
 * initialize matrix of costs for alignment of prefixes (a_1 ... a_i; b_1 ... b_j)
 */
void nwdp_initF( float ** F, int L1, int L2 )
{
	F[ 0 ][ 0 ] =  0. ;

	for( int j = 1; j <= L2; j++ )
		for( int i = 1; i <= L1; i++ )
			F[ j ][ i ] = INFINITE;

    for( int i = 1; i <= L1; i++ )
     	F[ 0 ][ i ] =  i * alpha;
    for( int j = 1; j <= L2; j++ )
      	F[ j ][ 0 ] =  j * alpha;
}


/*
 * initialize matrix of costs for alignment of prefixes (a_1 ... a_i; b_1 ... b_j)
 * for affine gap costs
 */
void nwdp_initF_affinegaps( float ** F, int L1, int L2, bool local )
{
	F[ 0 ][ 0 ] =  0. ;

	for( int j = 1; j <= L2; j++ )
		for( int i = 1; i <= L1; i++ )
			F[ j ][ i ] = INFINITE;

	for( int i = 1; i <= L1; i++ )
		F[ 0 ][ i ] = ( !local ) ? alpha + i * beta : 0;
	for( int j = 1; j <= L2; j++ )
		F[ j ][ 0 ] = ( !local ) ? alpha + j * beta : 0;
}


/*
 * initialize matrix of costs for alignment of prefixes (a_1 ... a_i; b_1 ... b_j)
 * that ends with a gap
 */
void nwdp_initGap( float ** Q, int L1, int L2 )
{
	Q[ 0 ][ 0 ] =  0. ;

	for( int j = 0; j <= L2; j++ )
		for( int i = 0; i <= L1; i++ )
			Q[ j ][ i ] = INFINITE;
}


void nwdp_initTB( char ** traceback, int L1, int L2 )
//void nwdp_initTB( vector< vector<char> > traceback, int L1, int L2 )
{
    traceback[ 0 ][ 0 ] = 'n' ;

    for( int i = 1; i <= L1; i++ )
    	traceback[ 0 ][ i ] =  '-' ;
    for( int j = 1; j <= L2 ; j++ )
    	traceback[ j ][ 0 ] =  '|' ;
}


float nwdb_global_align_affinegaps( float* probDbl_1, float* probSgl_1, int L1, float* probDbl_2, float* probSgl_2, int L2, bool freeEndGaps, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP )
{
    float tau, sigma;;
    char ptr;

    /* opening gap penalty */
    float gapopen = alpha + beta;

    /* Initialize dynamic programming matrix F (fill in first row and column) */
    nwdp_initF_affinegaps( F, L1, L2, 1 );

    /* Initialize Q and P matrices for cost of alignments that end with a gap */
    nwdp_initGap( Q, L1, L2 );
   	nwdp_initGap( P, L1, L2 );

    /* Initialize traceback matrices */
    nwdp_initTB( trF, L1, L2 );
    nwdp_initTB( trQ, L1, L2 );
    nwdp_initTB( trP, L1, L2 );

   	/* base pair similarity of bases pairing with nuc_1 and bases pairing with nuc_2 */
   	int j, i;
    for( j = 1; j <= L2; j++ )
    	for( i = 1; i <= L1; i++ )
        {
			// calculate P
			P[ j ][ i ] = max3( P[ j-1 ][ i ] + beta, F[ j-1 ][ i ] + gapopen, INFINITE, &ptr );
            trP[ j ][ i ] =  ptr;

			// calculate Q
			Q[ j ][ i ] = max3( INFINITE, F[ j ][ i-1 ] + gapopen, Q[ j ][ i-1 ] + beta, &ptr );
            trQ[ j ][ i ] =  ptr;

			// calculate F
			tau = ( probDbl_1[ i-1 ]==0 && probDbl_2[ j-1 ]==0 ) ? deltanull : 0.5 - abs( probDbl_1[ i-1 ] - probDbl_2[ j-1 ] );
    		tau *= ( 1 - kappa );
    		sigma = ( probSgl_1[ i-1 ]==0 && probSgl_2[ j-1 ]==0 ) ? deltanull : 0.5 - abs( probSgl_1[ i-1 ] - probSgl_2[ j-1 ] );
    		tau += kappa * sigma;
			F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + tau, Q[ j ][ i ], &ptr );
            trF[ j ][ i ] =  ptr;
        }
	j--; i--;

    if( freeEndGaps )
    {
		/* remove 3' tail gaps */
		bool tail = 0;
		int maxmat = 0;  // 0 .. F; 1 .. P; 2 .. Q
		while( i > 0 || j > 0 )
		{
			if( tail )
				break;

			switch( maxmat )
			{
				case 0 : switch( trF[ j ][ i ] )
						 {
							 case '|' : maxmat = 1;
										break;
							 case '\\': tail = 1;
										break ;
							 case '-' : maxmat = 2;
										break;
						 }
						 break;
				 case 1: switch( trP[ j ][ i ] )
						 {
							 case '|' : j--;
										break;
							 case '\\': maxmat = 0;
										j--;
										break;
						 }
						 break;
				 case 2: switch( trQ[ j ][ i ] )
						 {
							 case '\\': maxmat = 0;
										i--;
										break;
							 case '-' : i--;
										break;
						 }
						 break;
			}
		}
    }

	/* score with or without 3' tail gaps */
    return F[ j ][ i ];
}


/*
 * sequence similarity of nuc_1 and nuc_2
 */
float nwdb_align_seq_sim( char nuc_1, float probSgl_1, char nuc_2, float probSgl_2 )
{
	float sigma = 0;
	if( nuc_1 == nuc_2 )
		sigma = ( probSgl_1==0 && probSgl_2==0 ) ? deltanull : 0.5 - abs( probSgl_1 - probSgl_2 );

	return sigma;
}


float max3( float f1, float f2, float f3, char* ptr )
{
	float  max = 0;

    if( f1 >= f2 && f1 >= f3 )
    {
    	max = f1;
        *ptr = '|';
    }
    else if( f2 > f3 )
    {
    	max = f2;
    	*ptr = '\\';
    }
    else
    {
    	max = f3;
    	*ptr = '-';
    }

    return max;
}


float max( float f1, float f2 )
{
    return ( f1 >= f2 ) ? f1 : f2;
}


template <typename T>
void print_matrixdp( T ** F, float * prob_1, int L1, float * prob_2, int L2 )
{
    cout << "        ";
    for( int i = 0; i < L1; i++ )
    	cout << prob_1[ i ] << "   ";
    cout << "\n  ";

    for( int j = 0; j <= L2; j++ )
    {
    	if( j > 0 )
    		cout << prob_2[ j-1 ] << " ";
    	for( int i = 0; i <= L1; i++ )
        {
    		cout.width( 3 );
            cout << F[ j ][ i ] << " ";
        }
        cout << endl;
    }
}


float simalign_affinegaps( float ** Z, int L1, int L2, int * idx_1_aln, int * idx_2_aln, int & Lpair, int & Laln, bool global, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP )
{
    int i = 0, j = 0;
    float sim;

    /* opening gap penalty */
    float gapopen = alpha + beta;

    /* Initialize dynamic programming matrix */
    nwdp_initF_affinegaps( F, L1, L2, !global );

    /* Initialize Q and P matrices for cost of alignments that end with a gap */
    nwdp_initGap( Q, L1, L2 );
   	nwdp_initGap( P, L1, L2 );

    /* Initialize traceback matrices */
    nwdp_initTB( trF, L1, L2 );
    nwdp_initTB( trQ, L1, L2 );
    nwdp_initTB( trP, L1, L2 );

    /* base pair similarity of bases pairing with nuc_1 and bases pairing with nuc_2 */
    /* create alignment */
    char ptr;
    for( j = 1; j <= L2; j++ )
    	for( i = 1; i <= L1; i++ )
        {
    	    // calculate P
   			P[ j ][ i ] = max3( P[ j-1 ][ i ] + beta, F[ j-1 ][ i ] + gapopen, INFINITE, &ptr );
            if( !global && P[ j ][ i ] < 0 )
            	P[ j ][ i ] = 0;
            trP[ j ][ i ] =  ptr;

    		// calculate Q
           	Q[ j ][ i ] = max3( INFINITE, F[ j ][ i-1 ] + gapopen, Q[ j ][ i-1 ] + beta, &ptr );
            if( !global && Q[ j ][ i ] < 0 )
            	Q[ j ][ i ] = 0;
            trQ[ j ][ i ] =  ptr;

            // calculate F
            F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + Z[ j-1 ][ i-1 ], Q[ j ][ i ], &ptr );
            if( !global && F[ j ][ i ] < 0 )
            	F[ j ][ i ] = 0;
            trF[ j ][ i ] =  ptr;
        }
    sim = F[ L2 ][ L1 ];
    i--; j--;

    /* find maximal entry of similarity in matrix F, P or Q for local alignments */
    int maxmat = 0;  // 0 .. F; 1 .. P; 2 .. Q
    int li = L1;
    int lj = L2;
    if( !global ) {
    	float lmax = 0;
        for( j = 1; j <= L2; j++ )
        	for( i = 1; i <= L1; i++ ) {
        		if( F[ j ][ i ] > lmax ) {
        			lmax = F[ j ][ i ];
        			li = i;
        			lj = j;
        			maxmat = 0;
        		}
        		if( P[ j ][ i ] > lmax ) {
        			lmax = P[ j ][ i ];
        			li = i;
        			lj = j;
        			maxmat = 1;
        		}
        		if( Q[ j ][ i ] > lmax ) {
        			lmax = Q[ j ][ i ];
        			li = i;
        			lj = j;
        			maxmat = 2;
        		}
        	}
        i = li;
        j = lj;
        sim = lmax;
    }

    /* backtracking */
    int k = 0;
    int w = 0;
    while( i > 0 || j > 0 )
    {
    	if( !global )
    		if( (maxmat == 0 && F[ j ][ i ] == 0) || (maxmat == 1 && P[ j ][ i ] == 0) || (maxmat == 2 && Q[ j ][ i ] == 0) )
    			break;

    	switch( maxmat )
    	{
    		case 0 : switch( trF[ j ][ i ] )
    				 {
        			 	 case '|' : maxmat = 1;
        			 	 	 	 	break;
        			 	 case '\\': idx_1_aln[ k ] = i-1;
        			 	 	 	 	idx_2_aln[ k++ ] = j-1;
        			 	 	 	 	i--; j--;
        			 	 	 	 	w++;
        			 	 	 	 	break ;
        			 	 case '-' : maxmat = 2;
        			 	 	 	 	break;
    				 }
    				 break;
    		 case 1: switch( trP[ j ][ i ] )
    				 {
    		 	 	 	 case '|' : j--;
    		 	 	 	 	 	 	w++;
    		 	 	 	 	 	 	break;
    		 	 	 	 case '\\': maxmat = 0;
    		 	 	 	 	 	 	j--;
    		 	 	 	 	 	 	w++;
    		 	 	 	 	 	 	break;
    				 }
    				 break;
    		 case 2: switch( trQ[ j ][ i ] )
    		 	 	 {
 	 	 	 	 	 	 case '\\': maxmat = 0;
 	 	 	 	 	 	 	 	 	i--;
 	 	 	 	 	 	 	 	 	w++;
 	 	 	 	 	 	 	 	 	break;
 	 	 	 	 	 	 case '-' : i--;
 	 	 	 	 	 	 	 	 	w++;
 	 	 	 	 	 	 	 	 	break;
    		 	 	 }
    		 	 	 break;
    	}
    }

    Laln = w;
    Lpair = k;
    reverse( idx_1_aln, Lpair );
    reverse( idx_2_aln, Lpair );

	#if DEBUG
        cout << "\nDynamic programming matrix: " << endl;
        float* tempidx_1 = new float[L1+1];
        float* tempidx_2 = new float[L2+1];
        for( int k = 0; k <= L1; k++ ) tempidx_1[ k ] = (float) k;
        for( int k = 0; k <= L2; k++ ) tempidx_2[ k ] = (float) k;
        cout << "F:" << endl;
    	print_matrixdp( F, tempidx_1, L1, tempidx_2, L2);
    	print_matrixdp( trF, tempidx_1, L1, tempidx_2, L2);
        cout << "P:" << endl;
    	print_matrixdp( P, tempidx_1, L1, tempidx_2, L2);
    	print_matrixdp( trP, tempidx_1, L1, tempidx_2, L2);
        cout << "Q:" << endl;
    	print_matrixdp( Q, tempidx_1, L1, tempidx_2, L2);
    	print_matrixdp( trQ, tempidx_1, L1, tempidx_2, L2);
        free(tempidx_1);
        free(tempidx_2);
        for(int i=0; i<Lpair; i++){cerr << idx_1_aln[i] << " ";} cerr << endl;
        for(int i=0; i<Lpair; i++){cerr << idx_2_aln[i] << " ";} cerr << endl;
	#endif

	cerr << "STRUCTALN: Seq_1 ( " << i << ", " << li << " ); Seq_2 ( " << j << ", " << lj << " )" << endl;

 	return sim;
}


void affinegapcosts( int * idx_1_aln, int * idx_2_aln, int & len_aln, int & open, int & extended )
{
	int gaplen;

	for( int i=1; i<len_aln; i++ ) {
		gaplen = idx_1_aln[ i ] - idx_1_aln[ i-1 ] - 1;
		if( gaplen ) {
			open++;
			extended += gaplen;
		}

		gaplen = idx_2_aln[ i ] - idx_2_aln[ i-1 ] - 1;
		if( gaplen ) {
			open++;
			extended += gaplen;
		}
	}
}


void freeMatrix(float ** matrix, int leny)
{
	for(int i = 0 ; i < leny ; i++) {
		delete[] matrix[i];
	}
	delete[] matrix;
}


void reverse( int * list, int len )
{
	int temp;
	for( int i = 0; i < len/2; i++ )
	{
		temp = list[ i ];
	    list[ i ] = list[ len-i-1 ];
	    list[ len-i-1 ] = temp;
	}
}


/*
 * return overlap of the second interval by the first interval
 */
float getOverlap2ndInterval( int start_1, int end_1, int start_2, int end_2 )
{
	int len_2 = end_2 - start_2 + 1;

	if( end_1 >= start_2 )
		if( start_1 <= start_2 )
			if( end_1 <= end_2 )
				return (float) (end_1 - start_2 + 1) / len_2;
			else
				return 1.;
		else
			if( end_1 <= end_2)
				return (float) (end_1 - start_2 + 1 - start_1 + start_2) / len_2;
			else
				return (float) (end_2 - start_1 + 1) / len_2;
	else
		return 0.;
}


/*
 * return normalized log odds of base pair probabilities weighted by minimal considered paired probabilities
 */
void getlogoddsDblNeighborhood( vector<float> & probDbl, string seq, int len, float pnull, int radius, float theta )
{
	float * tmpprob = new float[ len * len ];
	float adjprob;
	int offset, offsetu, offsetd;
	int norm;

    for( int i = 0; i < len; i++ )
    {
    	offset = i * len;
    	for( int j = 0; j < len; j++ )
    		tmpprob[ offset + j ] = max( 0, log( probDbl[ offset + j ] / pnull ) / log( 1 / pnull ) );
    		//tmpprob[ offset + j ] = probDbl[ offset + j ];
    }

    for( int i = 0; i < len; i++ )
    {
    	offset = i * len;
    	for( int j = 0; j < len; j++ )
    	{
    		norm = 2 * radius * radius;
    		adjprob = 0;
    		for( int ry = 1; ry <= radius; ry++ )
    		{
				offsetu = ( i - ry ) * len;
				offsetd = ( i + ry ) * len;
    			for( int rx = 1; rx <= radius; rx++ )
    			{
    				if( offsetu < 0 || j + rx >= len )
    					norm --;
    				else
    					adjprob += tmpprob[ offsetu + j + rx ];

    				if( j - rx < 0 || offsetd >= len*len )
    					norm --;
    				else
    					adjprob += tmpprob[ offsetd + j - rx ];

    			}
    		}

        	probDbl[ offset + j ] = ( adjprob > 0 ) ? ( 1 - theta ) * tmpprob[ offset + j ] + theta * adjprob / norm : tmpprob[ offset + j ];
    	}
    }

    delete[] tmpprob;
}


/*
 * return normalized log odds of unpaired probabilities weighted by minimal considered unpaired probability
 */
void getlogoddsSgl( vector<float> & probSgl, string seq, int len, float pnull )
{
    for( int i = 0; i < len; i++ )
    {
    	//float pnulli = pnullSgl[ nucIdx.at( seq.c_str()[ i ] ) ];
    	//probSgl[ i ] = max( 0, log( probSgl[ i ] / pnulli ) / log( 1 / pnulli ) );
    	//probSgl[ i ] = ( probSgl[ i ] <= pnull ) ? log( pnull / pnulli ) / log( 1 / pnulli ) : log( probSgl[ i ] / pnulli ) / log( 1 / pnulli );
    	//probSgl[ i ] = ( probSgl[ i ] <= pnulli ) ? 0 : log( probSgl[ i ] / pnulli ) / log( 1 / pnulli );
    	probSgl[ i ] = max( 0, log( probSgl[ i ] / pnull ) / log( 1 / pnull ) );
    }
}


void printalign( string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln )
{
	int len_1 = seq_1.length();
	int len_2 = seq_2.length();

	char *seqar_1 = new char[len_1+1];
	strcpy(seqar_1, seq_1.c_str());
	char *seqar_2 = new char[len_2+1];
	strcpy(seqar_2, seq_2.c_str());
	string seq_1_aln, seq_2_aln;

	int offset_1 = 0, offset_2 = 0;
	for (int k=0; k < len_aln; k++) {
		while( idx_1_aln[ k ] != offset_1 ) { seq_1_aln.push_back( seqar_1[ offset_1 ] ); seq_2_aln.push_back( '-' ); offset_1++; }
		while( idx_2_aln[ k ] != offset_2 ) { seq_2_aln.push_back( seqar_2[ offset_2 ] ); seq_1_aln.push_back( '-' ); offset_2++; }
		seq_1_aln.push_back( seqar_1[ idx_1_aln[ k ] ] ); offset_1++;
		seq_2_aln.push_back( seqar_2[ idx_2_aln[ k ] ] ); offset_2++;
	}
	while( offset_1 < len_1 ) { seq_1_aln.push_back( seqar_1[ offset_1 ] ); seq_2_aln.push_back( '-' ); offset_1++; }
	while( offset_2 < len_2 ) { seq_2_aln.push_back( seqar_2[ offset_2 ] ); seq_1_aln.push_back( '-' ); offset_2++; }

	cout << seq_1_aln << endl;
	cout << seq_2_aln << endl;

    delete[] seqar_1;
    delete[] seqar_2;
}
