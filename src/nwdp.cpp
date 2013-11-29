/*
 * nwdp.cpp
 *
 *  Created on: Oct 28, 2013
 *      Author: Stefan Seemann, seemann@rth.dk
 */

#include "nwdp.h"

using namespace std;


int readinput ( istream & is, string & name, string & seq, vector<float> & prob )
{
	// read probability matrix
	string line;
	int lnr = 0;
	int i = -1;

	while( !is.eof() ) {
		getline(is, line);

		if( i == -1 ) {
			name = line.replace(0,1,"");
			//prob = (float) realloc (prob, lnr * lnr * sizeof(float));
		}
		else {
			if( i == 0 ) {
				seq = line;
				lnr = line.length();
				prob.reserve( lnr * lnr );
			}
			else {
				vector<string> probv = split(line, ' ');
				//int offset = (i-1)*lnr;
				for( unsigned int k = 0; k < probv.size(); k++ )
				{
					//prob[offset + k] = atof( probv[k].c_str() );
					prob.push_back( atof( probv[k].c_str() ) );
					//cout << offset+k << " " << prob[offset + k] << endl;
				}
			}
		}

		if( i == lnr )
			break;
		i++;
	}

	return lnr;
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


/*
 * adapted Needleman-Wunsch algorithm for global alignment of two probability dotplots
 */
float nwdp( vector<float> & prob_1, int sidx_1, int* idx_1_aln, int len_1, vector<float> & prob_2, int sidx_2, int* idx_2_aln, int len_2, float gappenalty, bool prm )
{
    float* subprob_1 = new float[len_1];
    float* subprob_2 = new float[len_2];

    // get basepair probabilities that should be aligned
    for( int i = 0; i < len_1; i++ )
    	subprob_1[ i ] = prob_1.at( sidx_1 + idx_1_aln[ i ] );
    for( int j = 0; j < len_2; j++ )
    	subprob_2[ j ] = prob_2.at( sidx_2 + idx_2_aln[ j ] );

    // Dynamic programming matrix
    float ** F = new float * [ len_2 + 1  ];
    for( int j = 0; j <= len_2; j++ )
    	F[ j ] = new float [ len_1 + 1 ];

    // Initialize traceback and F matrix (fill in first row and column)
    nwdp_initF( F, len_1, len_2, gappenalty );

    // Create alignment
    float distance = nwdb_align( F, subprob_1, len_1, subprob_2, len_2, gappenalty );

/*    if( prm )
    {
            cout << "\nDynamic programming matrix: " << endl;
            print_matrixdp( F, subprob_1, len_1, subprob_2, len_2 );
            cout << endl;
    }
*/
    for( int j = 0; j <= len_2; j++ )  delete[] F[ j ];
	delete[] F;
	delete[] subprob_1;
	delete[] subprob_2;

    return distance;
}


void nwdp_initF( float ** F, int L1, int L2, float d )
{
	F[ 0 ][ 0 ] =  0. ;

    for( int i = 1; i <= L1; i++ )
     	F[ 0 ][ i ] =  i * d ;
    for( int j = 1; j <= L2; j++ )
      	F[ j ][ 0 ] =  j * d ;
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


float nwdb_align( float ** F, float* prob_1, int L1, float* prob_2, int L2, int gap )
{
    float fU, fD, fL, bpdist;
    char ptr;

    for( int j = 1; j <= L2; j++ )
    {
            for( int i = 1; i <= L1; i++ )
            {
                    bpdist = abs( prob_1[ i-1 ] - prob_2[ j-1 ] );

                    fU = F[ j-1 ][ i ] + gap ;
                    fD = F[ j-1 ][ i-1 ] + bpdist ;
                    fL = F[ j ][ i-1 ] + gap ;

                    F[ j ][ i ] = min3( fU, fD, fL, &ptr ) ;
            }
    }

    return F[ L2 ][ L1 ];
}


float  min3( float f1, float f2, float f3, char* ptr )
{
        float  min = 0 ;

        if( f1 <= f2 && f1 <= f3 )
        {
                min = f1 ;
                *ptr = '|' ;
        }
        else if( f2 < f3 )
        {
                min = f2 ;
                *ptr = '\\' ;
        }
        else
        {
                min = f3 ;
                *ptr = '-' ;
        }

        return  min ;
}


float min( float f1, float f2 )
{
    return ( f1 <= f2 ) ? f1 : f2;
}


void print_matrixdp( float ** F, float * prob_1, int L1, float * prob_2, int L2 )
{
    cout << "        ";
    for( int i = 0; i < L1; i++ )
    {
            cout << prob_1[ i ] << "   ";
    }
    cout << "\n  ";

    for( int j = 0; j <= L2; j++ )
    {
            if( j > 0 )
            {
                    cout << prob_2[ j-1 ] << " ";
            }
            for( int i = 0; i <= L1; i++ )
            {
                    cout.width( 3 );
                    cout << F[ j ][ i ] << " ";
            }
            cout << endl;
    }
}


float nwdistance( float ** D, int L1, int L2, int * idx_1_aln, int * idx_2_aln, int & L1_aln, int & L2_aln, int precision, float gappenalty, bool prm)
{
    int i = 0, j = 0;
    float distance;

    /* get gap penalty as mean of all distances */
    float gap = 0;
    for( j = 0; j < L2; j++ )
        for( i = 0; i < L1; i++ )
        	gap += D[ j ][ i ];
    gap = gappenalty * gap / ( L1 * L2 );
    char p[10];
    sprintf(p, "%.*lf", precision, gap);
    gap = ( atof(p) > 0. ) ? atof(p) : 0.01;
	#if DEBUG
    	cout << "gap: " << gap << endl;
	#endif

    /* Dynamic programming matrix */
    float ** F = new float * [ L2+1 ];
    for( j = 0; j <= L2; j++ )
    	F[ j ] = new float [ L1+1 ];

    /* Traceback matrix */
    char ** traceback = new char * [ L2+1 ];
    for( j = 0; j <= L2; j++ )
    	traceback[ j ] = new char [ L1+1 ];

    /* Initialize traceback and F matrix (fill in first row and column) */
    for( i = 0; i <= L1; i++ )
     	F[ 0 ][ i ] = i*gap;
    for( j = 1; j <= L2; j++ )
      	F[ j ][ 0 ] = j*gap;
    nwdp_initTB( traceback, L1, L2 );

    /* create alignment */
    float fU, fD, fL;
    char ptr;

    for( j = 1; j <= L2; j++ )
    {
    	for( i = 1; i <= L1; i++ )
        {
    		fU = F[ j-1 ][ i ] + gap;
            fD = F[ j-1 ][ i-1 ] + D[ j-1 ][ i-1 ];
            fL = F[ j ][ i-1 ] + gap;

            F[ j ][ i ] = min3( fU, fD, fL, &ptr );
            traceback[ j ][ i ] =  ptr;
        }
    }
    distance = F[ L2 ][ L1 ];
    i--; j--;

    /* backtracking */
    int k = 0;
    while( i > 0 || j > 0 )
    {
    	//cout << j << " " << i << " " << traceback[ j ][ i ] << endl;
    	switch( traceback[ j ][ i ] )
        {
        	case '|' : j--;
                       break;
            case '\\': idx_1_aln[ k ] = i-1;
                       idx_2_aln[ k++ ] = j-1;
                       i--; j--;
                       break ;
            case '-' : i--;
                       break;
        }
    }

    L1_aln = L2_aln = k;
    reverse( idx_1_aln, L1_aln );
    reverse( idx_2_aln, L2_aln );

    if( prm )
    {
            cout << "\nDynamic programming matrix: " << endl;
            float* tempidx_1 = new float[L1+1];
            float* tempidx_2 = new float[L2+1];
            for( int k = 0; k <= L1; k++ ) tempidx_1[ k ] = (float) k;
            for( int k = 0; k <= L2; k++ ) tempidx_2[ k ] = (float) k;
            cout << "    ";
            for( i = 0; i <= L1; i++ ) cout << tempidx_1[ i ] << "   ";
            cout << endl;
            for( j = 0; j <= L2; j++ )
            {
                    cout << tempidx_2[ j ] << " ";
                    for( i = 0; i <= L1; i++ )
                    {
                            cout.width( 3 );
                            cout << F[ j ][ i ] << " ";
                    }
                    cout << endl;
            }
            cout << endl;

            cout << "  ";
            for( i = 0; i <= L1; i++ ) cout << tempidx_1[ i ] << " ";
            cout << endl;
            for( j = 0; j <= L2; j++ )
            {
                    cout << tempidx_2[ j ] << " ";
                    for( int i = 0; i <= L1; i++ )
                    {
                            cout << traceback[ j ][ i ] << " ";
                    }
                    cout << endl;
            }

            free(tempidx_1);
            free(tempidx_2);
    }

    /* free memory */
    for( int j = 0; j <= L2; j++ )  delete F[ j ];
    delete[] F;
    for( int j = 0; j <= L2; j++ )  delete traceback[ j ];
    delete[] traceback;

    return distance;
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


void printalign(string & seq_1, int * idx_1_aln, string & seq_2, int * idx_2_aln, int len_aln )
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
 * Needleman-Wunsch algorithm for global alignment of nt sequence
 */
int nw( string seq_1, string seq_2, string& seq_1_al, string& seq_2_al, bool prm)
{
        int  d = 5;                 /* gap penalty like in blastn existence gap costs */

        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        // Dynamic programming matrix
        int ** F = new int * [ L2+1 ];
        for( int i = 0; i <= L2; i++ )  F[ i ] = new int [ L1+1 ];

        // Traceback matrix
        char ** traceback = new char * [ L2+1 ];
        for( int i = 0; i <= L2; i++ )  traceback[ i ] = new char [ L1+1 ];

        // Initialize traceback and F matrix (fill in first row and column)
        dpm_init( F, traceback, L1, L2, d );

        // Create alignment
        int dist = nw_align( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, d );

        #if DEBUG
            int  L_al = seq_1_al.length();
            cout << "Length after alignment: " << L_al << endl;
        #endif

        if( prm )
        {
                cout << "\nDynamic programming matrix: " << "\n\n";
                print_matrix( F, seq_1, seq_2 );

                cout << "\nTraceback matrix: " << "\n\n";
                print_traceback( traceback, seq_1, seq_2 );

                cout << endl;
        }

        for( int i = 0; i <= L2; i++ )  delete F[ i ];
        delete[] F;
        for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
        delete[] traceback;

        return dist;
}


void  dpm_init( int ** F, char ** traceback, int L1, int L2, int d )
{
        F[ 0 ][ 0 ] =  0 ;
        traceback[ 0 ][ 0 ] = 'n' ;

        int i=0, j=0;

        for( j = 1; j <= L1; j++ )
        {
                F[ 0 ][ j ] =  -j * d ;
                traceback[ 0 ][ j ] =  '-' ;
        }
        for( i = 1; i <= L2; i++ )
        {
                F[ i ][ 0 ] =  -i * d ;
                traceback[ i ][ 0 ] =  '|' ;
        }
}


/*
 * Needleman-Wunsch algorithm
 */
int nw_align( int **F, char **traceback, string seq_1, string seq_2, string& seq_1_al, string& seq_2_al, int d )
{
        int        k = 0, x = 0, y = 0;
        int        fU, fD, fL ;
        char       ptr, nuc ;
        int        i = 0, j = 0;

        const int a = 2;   // Match like in blastn (optimized for sequence identity of 90%)
        const int b = -3;   // Mismatch like in blastn

        const int  s[ 4 ][ 4 ] = { { a, b, b, b },    /* substitution matrix */
                                   { b, a, b, b },
                                   { b, b, a, b },
                                   { b, b, b, a } } ;

        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

		for( i = 1; i <= L2; i++ )
        {
                for( j = 1; j <= L1; j++ )
                {
                        nuc = seq_1[ j-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  x = 0 ;  break ;
                                case 'C':  x = 1 ;  break ;
                                case 'G':  x = 2 ;  break ;
                                case 'T':  x = 3 ;  break ;
                        }

                        nuc = seq_2[ i-1 ] ;

                        switch( nuc )
                        {
                                case 'A':  y = 0 ;  break ;
                                case 'C':  y = 1 ;  break ;
                                case 'G':  y = 2 ;  break ;
                                case 'T':  y = 3 ;  break;
                        }

                        fU = F[ i-1 ][ j ] - d ;
                        fD = F[ i-1 ][ j-1 ] + s[ x ][ y ] ;
                        fL = F[ i ][ j-1 ] - d ;

                        F[ i ][ j ] = max3( fU, fD, fL, &ptr ) ;

                        traceback[ i ][ j ] =  ptr ;
                }
        }
        i-- ; j-- ;

        while( i > 0 || j > 0 )
        {
                switch( traceback[ i ][ j ] )
                {
                        case '|' :      seq_1_al += '-' ;
                                        seq_2_al += seq_2[ i-1 ] ;
                                        i-- ;
                                        break ;

                        case '\\':     	seq_1_al += seq_1[ j-1 ] ;
                                        seq_2_al += seq_2[ i-1 ] ;
                                        i-- ;  j-- ;
                                        break ;

                        case '-' :      seq_1_al += seq_1[ j-1 ] ;
                                        seq_2_al += '-' ;
                                        j-- ;
                                        break;
                }
                k++ ;
        }

        reverse( seq_1_al.begin(), seq_1_al.end() );
        reverse( seq_2_al.begin(), seq_2_al.end() );

        return  F[ L2 ][ L1 ] ;
}


int  max3( int f1, int f2, int f3, char * ptr )
{
        int  max = 0 ;

        if( f1 >= f2 && f1 >= f3 )
        {
                max = f1 ;
                *ptr = '|' ;
        }
        else if( f2 > f3 )
        {
                max = f2 ;
                *ptr = '\\' ;
        }
        else
        {
                max = f3 ;
                *ptr = '-' ;
        }

        return  max ;
}


void  print_matrix( int ** F, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        cout << "        ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << "   ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout.width( 3 );
                        cout << F[ i ][ j ] << " ";
                }
                cout << endl;
        }
}


void  print_traceback( char ** traceback, string seq_1, string seq_2 )
{
        int  L1 = seq_1.length();
        int  L2 = seq_2.length();

        cout << "    ";
        for( int j = 0; j < L1; j++ )
        {
                cout << seq_1[ j ] << " ";
        }
        cout << "\n  ";

        for( int i = 0; i <= L2; i++ )
        {
                if( i > 0 )
                {
                        cout << seq_2[ i-1 ] << " ";
                }
                for( int j = 0; j <= L1; j++ )
                {
                        cout << traceback[ i ][ j ] << " ";
                }
                cout << endl;
        }
}


void  print_al( string& seq_1_al, string& seq_2_al )
{
        cout << seq_1_al << endl;
        cout << seq_2_al << endl;
}




