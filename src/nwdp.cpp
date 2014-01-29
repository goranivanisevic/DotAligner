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
    cerr << "\n   -d --dotplot <file>   ... dotplot file (matrix of base pair probabilities)";
    cerr << "\n   Debug modes:";
	cerr << "\n   --printdp             ... print dynamic programming matrices";
    cerr << "\n   --global2             ... global alignment in step 2 (default is local alignment in step 2)";
    cerr << "\n   -k --kappa <float>    ... weight of sequence similarity (0..1; default = 0.5); 1 - kappa is weight of dotplots";
    cerr << "\n   -a --alpha <float>    ... affine gap costs = alpha + k * beta (k gaps; alpha default = -4)";
    cerr << "\n                             set <alpha> equal 0 to make the score becomes larger as a linear function of gap length";
    cerr << "\n   -b --beta <float>     ... affine gap costs = alpha + k * beta (k gaps; beta default = -1)";
    cerr << "\n                             set <beta> equal 0 to keep the score similar regardless of gap length";
    cerr << "\n   -p --precision <int>  ... number of digits considered of log-odds of base pair reliabilities (default = 4)";
    cerr << "\n   -m --maxshift <int>   ... maximal shift of two input sequences in the final alignment (default 50% of longer sequence)";
    cerr << "\n                             speeds up the calculation by ignoring pairwise comparisons of distant nucleotides";
    cerr << "\n                             local alignments may be missed for long sequences if <maxshift> is set too low";
    cerr << "\n   -s --seednr <int>     ... number of seed alignments (default = 3)";
    cerr << "\n   -l --seedlen <int>    ... length of seed alignments (default = 5 nucleotides)\n";
    cerr << "\n   -n --pnull <prob>     ... Minimal probability (default = 0.0005 )\n";
	cerr << "\n   --help                ... this output\n";
    cerr << "\n   Output:   similarity and alignment\n\n";

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
   \param[prm]  boolean if dynamic programming matrices are printed

   \returns  Similarity score ( 0 .. 1 )
*/
float nwdp( string seq_1, vector<float> & probDbl_1, vector<float> & probSgl_1, int rel_idx_1, int * idx_1_aln, int len_1, string seq_2, vector<float> & probDbl_2, vector<float> & probSgl_2, int rel_idx_2, int * idx_2_aln, int len_2, bool prm )
{
    float sim;
    float* subprob_1 = new float[len_1];
    float* subprob_2 = new float[len_2];
    string subseq_1, subseq_2;
    int idx_1 = idx_1_aln[ rel_idx_1 ];
    int idx_2 = idx_2_aln[ rel_idx_2 ];

    int sidx_1 = idx_1 * seq_1.length();
    int sidx_2 = idx_2 * seq_2.length();

    // get basepair probabilities that should be aligned
    for( int i = 0; i < len_1; i++ )
    	subprob_1[ i ] = probDbl_1.at( sidx_1 + idx_1_aln[ i ] );
    for( int j = 0; j < len_2; j++ )
    	subprob_2[ j ] = probDbl_2.at( sidx_2 + idx_2_aln[ j ] );
    // get subsequence that should be aligned
    for( int i = 0; i < len_1; i++ )
    	subseq_1 += seq_1.at( idx_1_aln[ i ] );
    for( int j = 0; j < len_2; j++ )
    	subseq_2 += seq_2.at( idx_2_aln[ j ] );

    /* create global alignment */
   	float tau = 0;

   	//tau = nwdb_constraint_global_align_affinegaps( subprob_1, probSgl_1, len_1, rel_idx_1, subprob_2, probSgl_2, len_2, rel_idx_2 );

   	/* alignment upstream of the constraint pair ( rel_idx_1, rel_idx_2 ) */
   	if( rel_idx_1 > 0 && rel_idx_2 > 0 )
   	{
    	float* subprobDbl_1u = new float[ rel_idx_1 ];
    	float* subprobDbl_2u = new float[ rel_idx_2 ];
    	float* subprobSgl_1u = new float[ rel_idx_1 ];
    	float* subprobSgl_2u = new float[ rel_idx_2 ];

    	for( int i = 0; i < rel_idx_1; i++ ) {
    		subprobDbl_1u[ i ] = probDbl_1.at( sidx_1 + idx_1_aln[ rel_idx_1 - 1 - i ] );
    		subprobSgl_1u[ i ] = probSgl_1.at( idx_1_aln[ rel_idx_1 - 1 - i ] );
    	}
   		for( int j = 0; j < rel_idx_2; j++ ) {
   			subprobDbl_2u[ j ] = probDbl_2.at( sidx_2 + idx_2_aln[ rel_idx_2 - 1 - j ] );
    		subprobSgl_2u[ j ] = probSgl_2.at( idx_2_aln[ rel_idx_2 - 1 - j ] );
    	}
   		tau = nwdb_global_align_affinegaps( subprobDbl_1u, subprobSgl_1u, rel_idx_1, subprobDbl_2u, subprobSgl_2u, rel_idx_2 );

   		delete[] subprobDbl_1u;
        delete[] subprobDbl_2u;
        delete[] subprobSgl_1u;
        delete[] subprobSgl_2u;
    }
    /* alignment downstream of the constraint pair ( rel_idx_1, rel_idx_2 ) */
    if( rel_idx_1<len_1-1 && rel_idx_2<len_2-1 )
    {
		float* subprobDbl_1d = new float[ len_1 - rel_idx_1 - 1 ];
		float* subprobDbl_2d = new float[ len_2 - rel_idx_2 - 1 ];
    	float* subprobSgl_1d = new float[ len_1 - rel_idx_1 - 1 ];
    	float* subprobSgl_2d = new float[ len_2 - rel_idx_2 - 1 ];

		for( int i = 0; i < len_1 - rel_idx_1 - 1; i++ ) {
			subprobDbl_1d[ i ] = probDbl_1.at( sidx_1 + idx_1_aln[ rel_idx_1 + i + 1 ] );
			subprobSgl_1d[ i ] = probSgl_1.at( idx_1_aln[ rel_idx_1 + i + 1 ] );
		}
		for( int j = 0; j < len_2 - rel_idx_2 - 1; j++ ) {
			subprobDbl_2d[ j ] = probDbl_2.at( sidx_2 + idx_2_aln[ rel_idx_2 + j + 1 ] );
    		subprobSgl_2d[ j ] = probSgl_2.at( idx_2_aln[ rel_idx_2 + j +1 ] );
		}
		tau += nwdb_global_align_affinegaps( subprobDbl_1d, subprobSgl_1d, len_1 - rel_idx_1 - 1, subprobDbl_2d, subprobSgl_2d, len_2 - rel_idx_2 - 1 );

		delete[] subprobDbl_1d;
	    delete[] subprobDbl_2d;
        delete[] subprobSgl_1d;
        delete[] subprobSgl_2d;
    }

    /* length of shorter sequence */
    int minlen = ( len_2 > len_1 ) ? len_1 : len_2;

    /* sequence similarity */
   	float sigma = nwdb_align_seq_sim( seq_1.c_str()[idx_1], probSgl_1.at(idx_1), seq_2.c_str()[idx_2], probSgl_2.at(idx_2) );

   	sim = kappa * sigma + ( 1 - kappa ) * tau / ( minlen - 1 );
   	//sim = kappa * sigma + ( 1 - kappa ) * tau;
   	//cerr << idx_1 << " " << idx_2 << " " << sigma << " " << tau << " " << minlen << " " << ( tau / minlen ) << " " << sim << endl;

   	delete[] subprob_1;
   	delete[] subprob_2;

    return ( sim > 0 ) ? sim : 0;
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

    /*for( int i = 1; i <= L1; i++ )
     	Q[ 0 ][ i ] =  INFINITE;
    for( int j = 1; j <= L2; j++ )
      	Q[ j ][ 0 ] =  INFINITE;*/
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


float nwdb_global_align_affinegaps( float* probDbl_1, float* probSgl_1, int L1, float* probDbl_2, float* probSgl_2, int L2 )
{
    float tau, sigma;;
    char ptr;
    float ** F, ** Q, ** P;
    float score;

    /* opening gap penalty */
    float gapopen = alpha + beta;

    /* Initialize dynamic programming matrix F (fill in first row and column) */
    F = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
    	F[ j ] = new float [ L1 + 1 ];
    nwdp_initF_affinegaps( F, L1, L2, 1 );

    /* Initialize Q and P matrices for cost of alignments that end with a gap */
    Q = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
    	Q[ j ] = new float [ L1 + 1 ];
    nwdp_initGap( Q, L1, L2 );

    P = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
      	P[ j ] = new float [ L1 + 1 ];
   	nwdp_initGap( P, L1, L2 );

    /* Traceback matrices */
    char ** trF = new char * [ L2+1 ];
    for( int j = 0; j <= L2; j++ )
    	trF[ j ] = new char [ L1+1 ];
    nwdp_initTB( trF, L1, L2 );
    char ** trQ = new char * [ L2+1 ];
    for( int j = 0; j <= L2; j++ )
    	trQ[ j ] = new char [ L1+1 ];
    nwdp_initTB( trQ, L1, L2 );
    char ** trP = new char * [ L2+1 ];
    for( int j = 0; j <= L2; j++ )
    	trP[ j ] = new char [ L1+1 ];
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
			tau = ( probDbl_1[ i-1 ]==0 && probDbl_2[ j-1 ]==0 ) ? 0 : 0.5 - abs( probDbl_1[ i-1 ] - probDbl_2[ j-1 ] );
    		tau *= ( 1 - kappa );
    		sigma = ( probSgl_1[ i-1 ]==0 && probSgl_2[ j-1 ]==0 ) ? 0 : 0.5 - abs( probSgl_1[ i-1 ] - probSgl_2[ j-1 ] );
    		tau += kappa * sigma;
			F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + tau, Q[ j ][ i ], &ptr );
            trF[ j ][ i ] =  ptr;
        }
	j--; i--;
    score = F[ j ][ i ];

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
	/* update score without 3' tail gaps */
    score = F[ j ][ i ];
    //cerr << "trF " << i << " " << j << " " << score << " " << F[ L2 ][ L1 ] << endl;

    /* free memory */
    for( int j = 0; j <= L2; j++ )  delete[] F[ j ];
  	delete[] F;
    for( int j = 0; j <= L2; j++ )  delete[] Q[ j ];
    delete[] Q;
    for( int j = 0; j <= L2; j++ )  delete[] P[ j ];
    delete[] P;
    for( int j = 0; j <= L2; j++ )  delete trF[ j ];
    delete[] trF;
    for( int j = 0; j <= L2; j++ )  delete trP[ j ];
    delete[] trP;
    for( int j = 0; j <= L2; j++ )  delete trQ[ j ];
    delete[] trQ;

    return score;
}


float nwdb_constraint_global_align_affinegaps( float* probDbl_1, vector<float> & probSgl_1, int L1, int constr_1, float* probDbl_2, vector<float> & probSgl_2, int L2, int constr_2 )
{
    float tau, sigma;
    char ptr;
    float ** F, ** Q, ** P;
    float score;

    // opening gap penalty
    float gapopen = alpha + beta;

    // Initialize dynamic programming matrix F (fill in first row and column)
    F = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
    	F[ j ] = new float [ L1 + 1 ];
    nwdp_initF_affinegaps( F, L1, L2, 1 );

    // Initialize Q and P matrices for cost of alignments that end with a gap
    Q = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
    	Q[ j ] = new float [ L1 + 1 ];
    nwdp_initGap( Q, L1, L2 );

    P = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
      	P[ j ] = new float [ L1 + 1 ];
   	nwdp_initGap( P, L1, L2 );

    /* Traceback matrices */
    char ** trF = new char * [ L2+1 ];
    for( int j = 0; j <= L2; j++ )
    	trF[ j ] = new char [ L1+1 ];
    nwdp_initTB( trF, L1, L2 );
    char ** trQ = new char * [ L2+1 ];
    for( int j = 0; j <= L2; j++ )
    	trQ[ j ] = new char [ L1+1 ];
    nwdp_initTB( trQ, L1, L2 );
    char ** trP = new char * [ L2+1 ];
    for( int j = 0; j <= L2; j++ )
    	trP[ j ] = new char [ L1+1 ];
    nwdp_initTB( trP, L1, L2 );

   	// base pair similarity of bases pairing with nuc_1 and bases pairing with nuc_2
   	int j = 1, i = 1;
    for( j = 1; j <= constr_2; j++ )
    	for( i = 1; i <= constr_1; i++ )
        {
    		// calculate P
    		P[ j ][ i ] = max3( P[ j-1 ][ i ] + beta, F[ j-1 ][ i ] + gapopen, INFINITE, &ptr );
            trP[ j ][ i ] =  ptr;

    		// calculate Q
            Q[ j ][ i ] = max3( INFINITE, F[ j ][ i-1 ] + gapopen, Q[ j ][ i-1 ] + beta, &ptr );
            trQ[ j ][ i ] =  ptr;

            tau = ( probDbl_1[ i-1 ]==0 && probDbl_2[ j-1 ]==0 ) ? 0 : 0.5 - abs( probDbl_1[ i-1 ] - probDbl_2[ j-1 ] );
    		tau *= ( 1 - kappa );
    		sigma = ( probSgl_1[ i-1 ]==0 && probSgl_2[ j-1 ]==0 ) ? 0 : 0.5 - abs( probSgl_1[ i-1 ] - probSgl_2[ j-1 ] );
    		tau += kappa * sigma;

            // calculate F
            F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + tau, Q[ j ][ i ], &ptr );
            trF[ j ][ i ] =  ptr;
        }
    if( constr_2+1 <= L2 && constr_1+1 <= L1 ) {
    	j = constr_2+1;
    	i = constr_1+1;
    	F[ j ][ i ] = ( constr_2 > 0 && constr_1 > 0 ) ? F[ constr_2 ][ constr_1 ] + 1 : 1;
    	trF[ j ][ i ] = '\\';
    	trQ[ j ][ i ] = '-';
    	trP[ j ][ i ] = '|';
    }
    if( constr_2+2 <= L2 && constr_1+2 <= L1 )
    {
		for( j = constr_2+2; j <= L2; j++ )
			for( i = constr_1+2; i <= L1; i++ )
			{
				// calculate P
				P[ j ][ i ] = max3( P[ j-1 ][ i ] + beta, F[ j-1 ][ i ] + gapopen, INFINITE, &ptr );
	            trP[ j ][ i ] =  ptr;

				// calculate Q
				Q[ j ][ i ] = max3( INFINITE, F[ j ][ i-1 ] + gapopen, Q[ j ][ i-1 ] + beta, &ptr );
	            trQ[ j ][ i ] =  ptr;

				tau = ( probDbl_1[ i-1 ]==0 && probDbl_2[ j-1 ]==0 ) ? 0 : 0.5 - abs( probDbl_1[ i-1 ] - probDbl_2[ j-1 ] );
	    		tau *= ( 1 - kappa );
	    		sigma = ( probSgl_1[ i-1 ]==0 && probSgl_2[ j-1 ]==0 ) ? 0 : 0.5 - abs( probSgl_1[ i-1 ] - probSgl_2[ j-1 ] );
	    		tau += kappa * sigma;

				// calculate F
				F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + tau, Q[ j ][ i ], &ptr );
	            trF[ j ][ i ] =  ptr;
			}
		j--; i--;
    }
    score = F[ j ][ i ];

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
	/* update score without 3' tail gaps */
    score = F[ j ][ i ];
    //cerr << "trF " << constr_1 << " " << constr_2 << " " << score << " " << F[ L2 ][ L1 ] << endl;

    // free memory
    for( int j = 0; j <= L2; j++ )  delete[] F[ j ];
  	delete[] F;
    for( int j = 0; j <= L2; j++ )  delete[] Q[ j ];
    delete[] Q;
    for( int j = 0; j <= L2; j++ )  delete[] P[ j ];
    delete[] P;
    for( int j = 0; j <= L2; j++ )  delete trF[ j ];
    delete[] trF;
    for( int j = 0; j <= L2; j++ )  delete trP[ j ];
    delete[] trP;
    for( int j = 0; j <= L2; j++ )  delete trQ[ j ];
    delete[] trQ;

    return score;
}


/*
 * sequence similarity of nuc_1 and nuc_2
 */
float nwdb_align_seq_sim( char nuc_1, float probSgl_1, char nuc_2, float probSgl_2 )
{
	float sigma = 0;
	if( nuc_1 == nuc_2 )
		sigma = ( probSgl_1==0 && probSgl_2==0 ) ? 0.9 : 1 - abs( probSgl_1 - probSgl_2 );

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


float simalign_affinegaps( float ** Z, int L1, int L2, int * idx_1_aln, int * idx_2_aln, int & Laln, int precision, bool global, bool prm)
{
    int i = 0, j = 0;
    float sim;

    // opening gap penalty
    float gapopen = alpha + beta;

    /* Dynamic programming matrix */
    float ** F = new float * [ L2+1 ];
    for( j = 0; j <= L2; j++ )
    	F[ j ] = new float [ L1+1 ];
    nwdp_initF_affinegaps( F, L1, L2, !global );

    // Initialize Q and P matrices for cost of alignments that end with a gap
    float ** Q = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
    	Q[ j ] = new float [ L1 + 1 ];
    nwdp_initGap( Q, L1, L2 );
    float ** P = new float * [ L2 + 1  ];
    for( int j = 0; j <= L2; j++ )
      	P[ j ] = new float [ L1 + 1 ];
   	nwdp_initGap( P, L1, L2 );

    /* Traceback matrices */
    char ** trF = new char * [ L2+1 ];
    for( j = 0; j <= L2; j++ )
    	trF[ j ] = new char [ L1+1 ];
    nwdp_initTB( trF, L1, L2 );
    char ** trQ = new char * [ L2+1 ];
    for( j = 0; j <= L2; j++ )
    	trQ[ j ] = new char [ L1+1 ];
    nwdp_initTB( trQ, L1, L2 );
    char ** trP = new char * [ L2+1 ];
    for( j = 0; j <= L2; j++ )
    	trP[ j ] = new char [ L1+1 ];
    nwdp_initTB( trP, L1, L2 );

    /* base pair similarity of bases pairing with nuc_1 and bases pairing with nuc_2 */
    /* create alignment */
    char ptr;
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
            F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + Z[ j-1 ][ i-1 ], Q[ j ][ i ], &ptr );
            if( !global && F[ j ][ i ] < 0 )
            	F[ j ][ i ] = 0;
            trF[ j ][ i ] =  ptr;
        }
    sim = F[ L2 ][ L1 ];
    i--; j--;

    /* find maximal entry of similarity in matrix F, P or Q for local alignments */
    int maxmat = 0;  // 0 .. F; 1 .. P; 2 .. Q
    if( !global ) {
    	float lmax = 0;
    	int li = 0;
    	int lj = 0;
        for( j = 1; j <= L2; j++ )
        	for( i = 1; i <= L1; i++ ) {
        		if( F[ j ][ i ] >= lmax ) {
        			lmax = F[ j ][ i ];
        			li = i;
        			lj = j;
        			maxmat = 0;
        		}
        		if( P[ j ][ i ] >= lmax ) {
        			lmax = P[ j ][ i ];
        			li = i;
        			lj = j;
        			maxmat = 1;
        		}
        		if( Q[ j ][ i ] >= lmax ) {
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

    Laln = k;
    reverse( idx_1_aln, Laln );
    reverse( idx_2_aln, Laln );

    if( prm ) {
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
    }

    /* free memory */
    for( int j = 0; j <= L2; j++ )  delete F[ j ];
    delete[] F;
    for( int j = 0; j <= L2; j++ )  delete P[ j ];
    delete[] P;
    for( int j = 0; j <= L2; j++ )  delete Q[ j ];
    delete[] Q;
    for( int j = 0; j <= L2; j++ )  delete trF[ j ];
    delete[] trF;
    for( int j = 0; j <= L2; j++ )  delete trP[ j ];
    delete[] trP;
    for( int j = 0; j <= L2; j++ )  delete trQ[ j ];
    delete[] trQ;

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
void getlogoddsDbl( vector<float> & probDbl, string seq, int len, float pnull )
{
    for( int i = 0; i < len; i++ )
    {
    	int offset = i * len;
    	for( int j = 0; j < len; j++ )
    	{
    		//float pnullij = pnullDbl[ 4*nucIdx.at( seq.c_str()[ i ] ) + nucIdx.at( seq.c_str()[ j ] ) ];
    		//probDbl[ offset + j ] = ( probDbl[ offset + j ] <= pnull ) ? log( pnull / pnullij ) / log( 1 / pnullij ) : log( probDbl[ offset + j ] / pnullij ) / log( 1 / pnullij );
    		//probDbl[ offset + j ] = ( probDbl[ offset + j ] <= pnullij ) ? 0 : log( probDbl[ offset + j ] / pnullij ) / log( 1 / pnullij );
    		probDbl[ offset + j ] = max( 0, log( probDbl[ offset + j ] / pnull ) / log( 1 / pnull ) );
    	}
    }
}


void getlogoddsDblNeighborhood( vector<float> & probDbl, string seq, int len, float pnull, int radius, float theta )
{
	vector<float> tmpprob;
	tmpprob.reserve( len * len );
	float adjprob;
	int offset, offsetu, offsetd;
	int norm;

    for( int i = 0; i < len; i++ )
    {
    	offset = i * len;
    	for( int j = 0; j < len; j++ ) {
    		//float pnullij = pnullDbl[ 4*nucIdx.at( seq.c_str()[ i ] ) + nucIdx.at( seq.c_str()[ j ] ) ];
    		//tmpprob[ offset + j ] = max( 0, log( probDbl[ offset + j ] / pnullij ) / log( 1 / pnullij ) );
    		tmpprob[ offset + j ] = max( 0, log( probDbl[ offset + j ] / pnull ) / log( 1 / pnull ) );
    	}
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

    		//cerr << i << " " << j << " " << tmpprob[ offset + j ] << " " << adjprob << " " << norm << endl;
        	probDbl[ offset + j ] = ( adjprob > 0 ) ? ( 1 - theta ) * tmpprob[ offset + j ] + theta * adjprob / norm : tmpprob[ offset + j ];
    	}
    }
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
