/*
 * nwdp_v2.cpp
 *
 *  Created on: Jul 03, 2014
 *      Author: Stefan Seemann, seemann@rth.dk
 */

#include "nwdp_v2.h"
#include <math.h>

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
    cerr << "\n   -k --kappa <float>     ... weight of sequence / unpaired probability similarity (0..1; default = 0.5); 1 - kappa is weight of dotplots";
    cerr << "\n   -t --tau <float>       ... weight of sequence similarity in sequence / unpaired probability similarity (0..1; default = 0.5)";
    cerr << "\n   -a --alpha <float>     ... affine gap costs = alpha + k * beta (k gaps; alpha default = -4)";
    cerr << "\n                              set <alpha> equal 0 to make the score becomes larger as a linear function of gap length";
    cerr << "\n   -b --beta <float>      ... affine gap costs = alpha + k * beta (k gaps; beta default = -1)";
    cerr << "\n                              set <beta> equal 0 to keep the score similar regardless of gap length";
    cerr << "\n   -s --subopt <int>      ... number of suboptimal sequence / unpaired probability alignments examined (default = 10)";
    cerr << "\n   -p --precision <int>   ... number of digits considered of log-odds of base pair reliabilities (default = 4)";
    cerr << "\n   -n --pnull <prob>      ... minimal probability (default = 0.0005 )";
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
   \param[probSgl_1]  unpaired probabilities of S_a
   \param[idx_1_aln]  indices of nucleotides in the subsequence of S_a that is considered for alignment ( SS_a )
   \param[len_1] length of considered subsequence SS_a
   \param[seq_2]  sequence 2 ( S_b )
   \param[probSgl_2]  unpaired probabilities of S_b
   \param[idx_2_aln]  indices of nucleotides in the subsequence of S_b that is considered for alignment ( SS_b )
   \param[len_2]  length of considered subsequence SS_b
   \param[len_pair]  number of aligned bases

   \returns  Similarity score ( 0 .. 1 )
*/
float nwdp( string seq_1, vector<float> & probSgl_1, int * idx_1_aln, int L1, string seq_2, vector<float> & probSgl_2, int * idx_2_aln, int L2, int & Laln, float ** F, float ** Q, float ** P, char ** trF, char ** trQ, char ** trP )
{
	int global = 1;
    int i = 0, j = 0;
    float score, sim;

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
    int x, y;
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
            switch( seq_1[ i-1 ] )
            {
            	case 'A':  x = 0 ;  break ;
                case 'C':  x = 1 ;  break ;
                case 'G':  x = 2 ;  break ;
                case 'T':  x = 3 ;  break ;
                case 'U':  x = 3 ;  break ;
            }
            switch( seq_2[ j-1 ] )
            {
            	case 'A':  y = 0 ;  break ;
                case 'C':  y = 1 ;  break ;
                case 'G':  y = 2 ;  break ;
                case 'T':  y = 3 ;  break ;
                case 'U':  y = 3 ;  break ;
            }

            score = tau * submatrix[ x ][ y ] + (1 - tau) * ( 1 - abs( probSgl_1[ i-1 ] - probSgl_2[ j-1 ] ) );
            //cerr << i << " " << j << " " << score << endl;
            F[ j ][ i ] = max3( P[ j ][ i ], F[ j-1 ][ i-1 ] + score, Q[ j ][ i ], &ptr );
            if( !global && F[ j ][ i ] < 0 )
            	F[ j ][ i ] = 0;
            trF[ j ][ i ] =  ptr;
        }
    sim = F[ L2 ][ L1 ];
    i--; j--;

    /* backtracking */
    int k = 0;
    int w = 0;
    int maxmat = 0;
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

    Laln = k;
    reverse( idx_1_aln, Laln );
    reverse( idx_2_aln, Laln );

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
        for(int i=0; i<Laln; i++){cerr << idx_1_aln[i] << " ";} cerr << endl;
        for(int i=0; i<Laln; i++){cerr << idx_2_aln[i] << " ";} cerr << endl;
	#endif

 	return sim/Laln;
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


float simdp( vector<float> & probDbl_1, int len_1, vector<float> & probDbl_2, int len_2, int * idx_1_aln, int * idx_2_aln, int len_aln )
{
	int offset_1, offset_2;
	float diff = 0;
	for( int i = 0; i < len_aln; i++ ) {
		offset_1 = idx_1_aln[ i ] * len_1;
		offset_2 = idx_2_aln[ i ] * len_2;
		for( int j = 0; j < len_aln; j++ ) {
			diff += abs( probDbl_1[ offset_1 + idx_1_aln[ j ] ] - probDbl_2[ offset_2 + idx_2_aln[ j ] ] );
		}
	}
	//cerr << diff << " " << diff/len_aln/len_aln << endl;

	return ( 1 - diff/len_aln/len_aln );
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
 * return normalized log odds of base pair probabilities weighted by minimal considered paired probabilities
 */
void getlogoddsDbl( vector<float> & probDbl, string seq, int len, float pnull )
{
	int offset;

    for( int i = 0; i < len; i++ )
    {
    	offset = i * len;
    	for( int j = 0; j < len; j++ )
    		probDbl[ offset + j ] = max( 0, log( probDbl[ offset + j ] / pnull ) / log( 1 / pnull ) );
    }
}


/*
 * return normalized log odds of unpaired probabilities weighted by minimal considered unpaired probability
 */
void getlogoddsSgl( vector<float> & probSgl, string seq, int len, float pnull )
{
    for( int i = 0; i < len; i++ )
    	probSgl[ i ] = max( 0, log( probSgl[ i ] / pnull ) / log( 1 / pnull ) );
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


/*
 * generate a random number from 0 to sum of entries in sam_prob
 */
int statistical_sampling(unsigned int *sam_prob, int sam_len)
{
	long sum = 0;
	int i, sample;

	for( i=0; i<sam_len; i++ ) {
		sum += sam_prob[i];
	}

	sample = rand()%sum;
	//printf("Sample: %i\n",sample);

	sum = 0;
	for( i=0; i<sam_len; i++ )
	{
		sum += sam_prob[i];
		if (sample<=sum) {
			return i;
		}
	}

	return 0;
}


/*
 *
 */
void prob_backtracking( int * idx_1_aln, int len_1, int * idx_2_aln, int len_2, int & Laln, float ** F )
{
	/* run through dynamic programming matrix using probabilistic sampling */
	int i = len_1;
	int j = len_2;
    int k = 0;
    int w = 0;
    int samprob;
    unsigned int btscores[3];

    while( i > 0 && j > 0 )
    {
   		btscores[0] = round(F[ j ][ i ]*10);
   		btscores[1] = ( round(F[ j-1 ][ i ]) > 0 ) ? round(F[ j-1 ][ i ]*10) : 0;
   		btscores[2] = ( round(F[ j ][ i-1 ]) > 0 ) ? round(F[ j ][ i-1 ]*10) : 0;

   		samprob = statistical_sampling(btscores, 3);
   		//cerr << i << " " << j << " " << btscores[0] << " " << btscores[1] << " " << btscores[2] << " " << samprob << endl;

   		switch( samprob )
   		{
   			case 0:	/* base pair */
   					idx_1_aln[ k ] = i-1;
   					idx_2_aln[ k++ ] = j-1;
   					i--; j--;
   					w++;
   					break;
   			case 1: /* base i unmatched */
   					j--;
   					w++;
   					break;
   			case 2: /* base j unmatched */
   					i--;
   					w++;
   					break;
   			default:  abort();
   		}
    }

    Laln = k;
    reverse( idx_1_aln, Laln );
    reverse( idx_2_aln, Laln );
}


/*
 * calculate sequence / unpaired similarity score of suboptimal alignment
 */
float simbp( string seq_1, vector<float> & probSgl_1, int len_1, string seq_2, vector<float> & probSgl_2, int len_2, int * idx_1_aln, int * idx_2_aln, int Laln )
{
	float sim = 0;

	/* sum scores of aligned bases */
    int x, y;
    for( int i = 0; i < Laln; i++ )
    {
    	switch( seq_1[ idx_1_aln[ i ] ] )
    	{
        	case 'A':  x = 0 ;  break ;
            case 'C':  x = 1 ;  break ;
            case 'G':  x = 2 ;  break ;
            case 'T':  x = 3 ;  break ;
            case 'U':  x = 3 ;  break ;
    	}
        switch( seq_2[ idx_2_aln[ i ] ] )
        {
          	case 'A':  y = 0 ;  break ;
          	case 'C':  y = 1 ;  break ;
            case 'G':  y = 2 ;  break ;
            case 'T':  y = 3 ;  break ;
            case 'U':  y = 3 ;  break ;
        }
        //cerr << idx_1_aln[ i ] << " " << idx_2_aln[ i ] << " " << seq_1[ idx_1_aln[ i ] ] << " " << seq_2[ idx_2_aln[ i ] ] << " " << submatrix[ x ][ y ] << endl;
        sim += tau * submatrix[ x ][ y ] + (1 - tau) * ( 1 - abs( probSgl_1[ idx_1_aln[ i ] ] - probSgl_2[ idx_2_aln[ i ] ] ) );
    }

    /* add scores of gaps */
    sim += gappenalty( idx_1_aln, Laln, len_1 );
    sim += gappenalty( idx_2_aln, Laln, len_2 );

    return sim/Laln;
}


float gappenalty( int * idx_aln, int Laln, int len )
{
    int diff;
    float sim = 0;

    /* opening gap penalty */
    float gapopen = alpha + beta;

    /* 5' gaps */
    diff = idx_aln[ 0 ];
    if( diff )
    	sim = gapopen + ( diff - 1 ) * beta;

    for( int i = 0; i < Laln-1; i++ )
    {
    	diff = idx_aln[ i+1 ] - idx_aln[ i ] - 1;
    	if( diff )
    		sim += gapopen + ( diff - 1 ) * beta;
    }

    /* 3' gaps */
    diff = len - idx_aln[ Laln-1 ] - 1;
    if( diff )
    	sim += gapopen + ( diff - 1 ) * beta;

    return sim;
}
