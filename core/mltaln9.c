#include "mltaln.h"

#define DEBUG 0

#ifdef PCALLS
#define CALLS 1
#else
#define CALLS 0
#endif

int seqlen( char *seq )
{
	int val = 0;
	if( *newgapstr == '-' )
	{
		while( *seq )
			if( *seq++ != '-' ) val++;
	}
	else
	{
		while( *seq )
		{
			if( *seq != '-' && *seq != *newgapstr ) val++;
			seq++;
		}
	}
	return( val );
}

int intlen( int *num )
{
	int value = 0;
	while( *num++ != -1 ) value++;
	return( value );
}

char seqcheck( char **seq )
{
	CALLS && printf("called %s:seqcheck()\n", __FILE__);
	int i, len;
	char **seqbk = seq;
	while( *seq )	
	{
		len = strlen( *seq );
		for( i=0; i<len; i++ ) 
		{
			if( amino_n[(int)(*seq)[i]] == -1 ) 
			{

				reporterr(       "========================================================================= \n" );
				reporterr(       "========================================================================= \n" );
				reporterr(       "=== \n" );
				reporterr(       "=== Alphabet '%c' is unknown.\n", (*seq)[i] );
				reporterr(       "=== Please check site %d in sequence %d.\n", i+1, (int)(seq-seqbk+1) );
				reporterr(       "=== \n" );
				reporterr(       "=== To make an alignment having unusual characters (U, @, #, etc), try\n" );
				reporterr(       "=== %% mafft --anysymbol input > output\n" );
				reporterr(       "=== \n" );
				reporterr(       "========================================================================= \n" );
				reporterr(       "========================================================================= \n" );
				return( (int)(*seq)[i] );
			}
		}
		seq++;
	}
	return( 0 );
}
	
void scmx_calc( int icyc, char **aseq, double *effarr, float **scmx )
{
	int  i, j, lgth;
	 
	lgth = strlen( aseq[0] );
	for( j=0; j<lgth; j++ )
	{
		for( i=0; i<nalphabets; i++ )
		{
			scmx[i][j] = 0;
		}
	}
	for( i=0; i<icyc+1; i++ )
	{
		int id;
		id = amino_n[(int)aseq[i][0]];
		scmx[id][0] += (float)effarr[i];
	}
	for( j=1; j<lgth-1; j++ )
	{
		for( i=0; i<icyc+1; i++ )
		{
			int id;
			id = amino_n[(int)aseq[i][j]];
			scmx[id][j] += (float)effarr[i];
		}
	}
	for( i=0; i<icyc+1; i++ )
	{
		int id;
		id = amino_n[(int)aseq[i][lgth-1]];
		scmx[id][lgth-1] += (float)effarr[i];
	}
}

void exitall( char arr[] )
{
	reporterr(       "%s\n", arr );
	exit( 1 );
}

void display( char **seq, int nseq )
{
	int i, imax;
	char b[121];

	if( !disp ) return;
		if( nseq > DISPSEQF ) imax = DISPSEQF;
		else                  imax = nseq;
		reporterr(       "    ....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+\n" );
		for( i=0; i<+imax; i++ )
		{
			strncpy( b, seq[i]+DISPSITEI, 120 );
			b[120] = 0;
			reporterr(       "%3d %s\n", i+1, b );
		}
}


void intergroup_score_consweight( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	int ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];



	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* なぜか配列を使わないとおかしくなる, 多分バグ */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (int)mseq1[k];
				ms2 = (int)mseq2[k];
				if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
				tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
	
				if( ms1 == (int)'-' ) 
				{
					tmpscore += (double)penalty;
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms1=(int)mseq1[++k]) == (int)'-' )
						;
//						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == (int)'-' )
				{
					tmpscore += (double)penalty;
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms2=(int)mseq2[++k]) == (int)'-' )
						;
//						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
//			reporterr(       "val in _gapnomi = %f\n", *value );
		}
	}
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}
void intergroup_score_gapnomi( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	int ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];



	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* なぜか配列を使わないとおかしくなる, 多分バグ */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (int)mseq1[k];
				ms2 = (int)mseq2[k];
				if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
//				tmpscore += (double)amino_dis[ms1][ms2];
	
				if( ms1 == (int)'-' ) 
				{
					tmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					while( (ms1=(int)mseq1[++k]) == (int)'-' )
						;
//						tmpscore += (double)amino_dis[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == (int)'-' )
				{
					tmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					while( (ms2=(int)mseq2[++k]) == (int)'-' )
						;
//						tmpscore += (double)amino_dis[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
//			reporterr(       "val in _gapnomi = %f\n", *value );
		}
	}
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}

void intergroup_score_multimtx( int **whichmtx, double ***scoringmatrices, char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k, c;
	int len2 = len - 2;
	int mn1, mn2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;
	int gapnum = amino_n['-'];

	double gaptmpscore;
	double gapscore = 0.0;

//	reporterr(       "#### in intergroup_score\n" );

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

//	reporterr(       "\n intergroup_score_multimtx ..." );
	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* なぜか配列を使わないとおかしくなる, 多分バグ */
			c = whichmtx[i][j];
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			gaptmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				mn1 = amino_n[(int)(mseq1[k])];
				mn2 = amino_n[(int)(mseq2[k])];
				if( mn1 == gapnum && mn2 == gapnum ) continue;
				tmpscore += (double)scoringmatrices[c][mn1][mn2];
//				tmpscore += (double)scoringmtx[mn1][mn2];
	
				if( mn1 == gapnum ) 
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)scoringmtx[mn1][mn2];
					tmpscore += (double)scoringmatrices[c][mn1][mn2];
					while( (mn1=amino_n[(int)mseq1[++k]]) == gapnum )
						tmpscore += (double)scoringmatrices[c][mn1][mn2];
//						tmpscore += (double)scoringmtx[mn1][mn2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( mn2 == gapnum )
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
					tmpscore += (double)scoringmatrices[c][mn1][mn2];
//					tmpscore += (double)scoringmtx[mn1][mn2];
					while( (mn2=amino_n[(int)mseq2[++k]]) == gapnum )
						tmpscore += (double)scoringmatrices[c][mn1][mn2];
//						tmpscore += (double)scoringmtx[mn1][mn2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
			gapscore += (double)gaptmpscore * (double)efficient;
		}
	}
//	reporterr(       "done." );
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}
void intergroup_score_dynmtx( double **offsetmtx, int scoringmtx[0x80][0x80], char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	int ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

	double gaptmpscore;
	double gapscore = 0.0;

//	reporterr(       "#### in intergroup_score\n" );

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

	reporterr(       "\n intergroup_score_dynmtx ..." );
	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* なぜか配列を使わないとおかしくなる, 多分バグ */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			gaptmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (int)mseq1[k];
				ms2 = (int)mseq2[k];
				if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
				tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//				tmpscore += (double)scoringmtx[ms1][ms2];
	
				if( ms1 == (int)'-' ) 
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)scoringmtx[ms1][ms2];
					tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;;
					while( (ms1=(int)mseq1[++k]) == (int)'-' )
						tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//						tmpscore += (double)scoringmtx[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == (int)'-' )
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
					tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//					tmpscore += (double)scoringmtx[ms1][ms2];
					while( (ms2=(int)mseq2[++k]) == (int)'-' )
						tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//						tmpscore += (double)scoringmtx[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
			gapscore += (double)gaptmpscore * (double)efficient;
		}
	}
	reporterr(       "done." );
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}
void intergroup_score( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	int ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

	double gaptmpscore;
	double gapscore = 0.0;

//	reporterr(       "#### in intergroup_score\n" );

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* なぜか配列を使わないとおかしくなる, 多分バグ */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			gaptmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (int)mseq1[k];
				ms2 = (int)mseq2[k];
				if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
//				tmpscore += (double)amino_dis[ms1][ms2];
				tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
	
				if( ms1 == (int)'-' ) 
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms1=(int)mseq1[++k]) == (int)'-' )
//						tmpscore += (double)amino_dis[ms1][ms2];
						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == (int)'-' )
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms2=(int)mseq2[++k]) == (int)'-' )
//						tmpscore += (double)amino_dis[ms1][ms2];
						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
			gapscore += (double)gaptmpscore * (double)efficient;
		}
	}
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}
void intergroup_score_new( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	int ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	static double efficient[1];

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			*efficient = eff1[i] * eff2[j]; /* なぜか配列を使わないとおかしくなる, 多分バグ */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (int)mseq1[k];
				ms2 = (int)mseq2[k];
				if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
				tmpscore += (double)amino_dis[ms1][ms2];
	
				if( ms1 == (int)'-' ) 
				{
					tmpscore += (double)penalty;
					tmpscore += (double)amino_dis[ms1][ms2];
					while( (ms1=(int)mseq1[++k]) == (int)'-' )
						tmpscore += (double)amino_dis[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == (int)'-' )
				{
					tmpscore += (double)penalty;
					tmpscore += (double)amino_dis[ms1][ms2];
					while( (ms2=(int)mseq2[++k]) == (int)'-' )
						tmpscore += (double)amino_dis[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)*efficient;
		}
	}
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}


double score_calc5( char **seq, int s, double **eff, int ex )  /* method 3 deha nai */
{
    int i, j, k;
    double c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;
    double efficient;
#if DEBUG
	FILE *fp;
#endif

    score = 0.0;
    c = 0.0;

	for( i=0; i<s; i++ ) 
	{
		
			if( i == ex ) continue;
            efficient = eff[i][ex];
            mseq1 = seq[i];
            mseq2 = seq[ex];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore * efficient;
/*
			fprintf( stdout, "%d-%d tmpscore = %f, eff = %f, tmpscore*eff = %f\n", i, ex, tmpscore, efficient, tmpscore*efficient );
*/
	}
	/*
	fprintf( stdout, "total score = %f\n", score );
	*/

    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
			if( i == ex || j == ex ) continue;

            efficient = eff[i][j];
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore * efficient;
        }
    }
/*
	reporterr(       "score in score_calc5 = %f\n", score );
*/
    return( (double)score );
/*

fprintf( trap_g, "score by fast = %f\n", (float)score );

tmpscore = score = 0.0;
	for( i=0; i<s; i++ ) 
	{
		if( i == ex ) continue;
		tmpscore = Cscore_m_1( seq, i, eff );
		fprintf( stdout, "%d %f\n", i, tmpscore );

		score += tmpscore;
	}
	tmpscore = Cscore_m_1( seq, ex, eff );
	fprintf( stdout, "ex%d %f\n", i, tmpscore );
	score += tmpscore;

	return( score );
*/
}


	
double score_calc4( char **seq, int s, double **eff, int ex )  /* method 3 deha nai */
{
    int i, j, k;
	double c;
    int len = strlen( seq[0] );
    double score;
    long tmpscore;
	char *mseq1, *mseq2;
	double efficient;

    score = 0.0;
	c = 0.0;
/*
	printf( "in score_calc4\n" );
	for( i=0; i<s; i++ )
	{
		for( j=0; j<s; j++ ) 
		{
			printf( "% 5.3f", eff[i][j] ); 
		}
		printf( "\n" );
		
	}
*/
    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
			efficient = eff[i][j];
			if( mix == 1 ) efficient = 1.0;
			/*
			printf( "weight for %d v.s. %d = %f\n", i, j, efficient );
			*/
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]] + 400 * !scoremtx ;

				c += efficient;

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty - n_dis[24][0];
                    while( mseq1[++k] == '-' )
						;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty - n_dis[24][0];
                    while( mseq2[++k] == '-' )
						;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
			/*
			if( x == 65 ) printf( "i=%d j=%d tmpscore=%d l=%d\n", i, j, tmpscore, len );
			*/
            score += (double)tmpscore * efficient;
        }
    }
    score /= c;
    return( (double)score );
}



void upg2( int nseq, double **eff, int ***topol, double **len )
{
    int i, j, k;
	double tmplen[M];

    static char **pair = NULL;

	if( !pair )
	{
		pair = AllocateCharMtx( njob, njob );
	}

	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        float minscore = 9999.0;
        int im = -1, jm = -1;
        int count;

        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
        {
            if( eff[i][j] < minscore )
            {
                minscore = eff[i][j];
                im = i; jm = j;
            }
        }
        for( i=0, count=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
            {
                topol[k][0][count] = i;
                count++;
            }
        topol[k][0][count] = -1;
        for( i=0, count=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
            {
                topol[k][1][count] = i;
                count++;
            }
        topol[k][1][count] = -1;

		len[k][0] = minscore / 2.0 - tmplen[im];
		len[k][1] = minscore / 2.0 - tmplen[jm];

		tmplen[im] = minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
                eff[MIN(i,im)][MAX(i,im)] =
                ( eff[MIN(i,im)][MAX(i,im)] + eff[MIN(i,jm)][MAX(i,jm)] ) / 2.0;
                eff[MIN(i,jm)][MAX(i,jm)] = 9999.0;
            }
            eff[im][jm] = 9999.0;
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

static void setnearest( int nseq, Bchain *acpt, float **eff, float *mindisfrompt, int *nearestpt, int pos )
{
	CALLS && printf("called %s:freeconstants()\n", __FILE__);
	int j;
	float tmpfloat;
	float mindisfrom;
	int nearest;
//	float **effptpt;
	Bchain *acptj;

	mindisfrom = 999.9;
	nearest = -1;

//	if( (acpt+pos)->next ) effpt = eff[pos]+(acpt+pos)->next->pos-pos;

//	for( j=pos+1; j<nseq; j++ )
	for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpfloat=*effpt++) < *mindisfrompt )
		if( (tmpfloat=eff[pos][j-pos]) < mindisfrom )
		{
			mindisfrom = tmpfloat;
			nearest = j;
		}
	}
//	effptpt = eff;
//	for( j=0; j<pos; j++ )
	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpfloat=(*effptpt++)[pos-j]) < *mindisfrompt )
		if( (tmpfloat=eff[j][pos-j]) < mindisfrom )
		{
			mindisfrom = tmpfloat;
			nearest = j;
		}
	}

	*mindisfrompt = mindisfrom;
	*nearestpt = nearest;
}

static void setnearest_double_fullmtx( int nseq, Bchain *acpt, double **eff, double *mindisfrompt, int *nearestpt, int pos )
{
	int j;
	double tmpfloat;
	double **effptpt;
	Bchain *acptj;

	*mindisfrompt = 999.9;
	*nearestpt = -1;

//	if( (acpt+pos)->next ) effpt = eff[pos]+(acpt+pos)->next->pos-pos;

//	for( j=pos+1; j<nseq; j++ )
	for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpfloat=*effpt++) < *mindisfrompt )
		if( (tmpfloat=eff[pos][j]) < *mindisfrompt )
		{
			*mindisfrompt = tmpfloat;
			*nearestpt = j;
		}
	}
	effptpt = eff;
//	for( j=0; j<pos; j++ )
	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpfloat=(*effptpt++)[pos-j]) < *mindisfrompt )
		if( (tmpfloat=eff[j][pos]) < *mindisfrompt )
		{
			*mindisfrompt = tmpfloat;
			*nearestpt = j;
		}
	}
}



static void loadtreeoneline( int *ar, float *len, FILE *fp )
{
	static char gett[1000];
	int res;
	char *p;

	p = fgets( gett, 999, fp );
	if( p == NULL )
	{
		reporterr(       "\n\nFormat error (1) in the tree?  It has to be a bifurcated and rooted tree.\n" );
		reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
		exit( 1 );
	}


	res = sscanf( gett, "%d %d %f %f", ar, ar+1, len, len+1 );
	if( res != 4 )
	{
		reporterr(       "\n\nFormat error (2) in the tree?  It has to be a bifurcated and rooted tree.\n" );
		reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
		exit( 1 );
	}

	ar[0]--;
	ar[1]--;

	if( ar[0] >= ar[1] )
	{
		reporterr(       "\n\nIncorrect guide tree\n" );
		reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
		exit( 1 );
	}


//	reporterr(       "ar[0] = %d, ar[1] = %d\n", ar[0], ar[1] );
//	reporterr(       "len[0] = %f, len[1] = %f\n", len[0], len[1] );
}

void loadtop( int nseq, float **mtx, int ***topol, float **len, char **name, int *nlen, Treedep *dep )
{
	int i, j, k, minijm, maxijm;
	int *intpt, *intpt2;
	int *hist = NULL;
	Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar;
	int nmemim, nmemjm;
	char **tree;
	char *treetmp;
	char *nametmp, *nameptr, *tmpptr; 
	char namec;
	FILE *fp;
	int node[2];
	float *height;
	float clusterdist;
	int mpair, mi, mj;

	fp = fopen( "_guidetree", "r" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( nseq );
		ac = (Bchain *)malloc( nseq * sizeof( Bchain ) );
		nmemar = AllocateIntVec( nseq );
//		treetmp = AllocateCharVec( nseq*50 );
		treetmp = NULL;
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( nseq, nseq*50 );
		tree = AllocateCharMtx( nseq, 0 );
		height = AllocateFloatVec( nseq );
	}

	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}


	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;


	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );
		len[k][0] = len[k][1] = -1.0;
		loadtreeoneline( node, len[k], fp );
		im = node[0];
		jm = node[1];

		if( im > nseq-1 || jm > nseq-1 || tree[im] == NULL || tree[jm] == NULL )
		{
			reporterr(       "\n\nCheck the guide tree.\n" );
			reporterr(       "im=%d, jm=%d\n", im+1, jm+1 );
			reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
			exit( 1 );
		}

		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];

//		reporterr(       "prevnode = %d, nmemim = %d\n", prevnode, nmemim );

		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


		nmemjm = nmemar[jm];
		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;

//		reporterr(       "prevnode = %d, nmemjm = %d\n", prevnode, nmemjm );

		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


//		len[k][0] = ( minscore - tmptmplen[im] );
//		len[k][1] = ( minscore - tmptmplen[jm] );
//		len[k][0] = -1;
//		len[k][1] = -1;


		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;


		if( len[k][0] == -1 || len[k][1] == -1 )
		{
			reporterr( "Re-computing the length of branch %d..\n", k );
			clusterdist = 0.0;
			mpair = 0;
			for( i=0; (mi=topol[k][0][i])>-1; i++ ) for( j=0; (mj=topol[k][1][j])>-1; j++ ) 
			{
				minijm = MIN(mi,mj);
				maxijm = MAX(mi,mj);
				clusterdist += mtx[minijm][maxijm-minijm];
				mpair += 1;
			}
			clusterdist /= (float)mpair;
			reporterr( "clusterdist = %f\n", clusterdist );
			if( len[k][0] == -1 ) len[k][0] = clusterdist/2.0 - height[im];
			if( len[k][1] == -1 ) len[k][1] = clusterdist/2.0 - height[im];
	
			fprintf( stderr, "len0 = %f\n", len[k][0] );
			fprintf( stderr, "len1 = %f\n\n", len[k][1] );
		}

		height[im] += len[k][0]; // for ig tree, 2015/Dec/25
		dep[k].distfromtip = height[im]; // for ig tree, 2015/Dec/25
//		reporterr( "##### dep[%d].distfromtip = %f\n", k, height[im] );



		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

//		reporterr(       "im,jm=%d,%d\n", im, jm );
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;


    }
	fclose( fp );
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
		fprintf( fp, "#by loadtop\n" );
	fclose( fp );

	FreeCharMtx( tree );
	free( treetmp );
	free( nametmp );
	free( hist );
	free( (char *)ac );
	free( (void *)nmemar );
	free( height );

}

void loadtree( int nseq, int ***topol, float **len, char **name, int *nlen, Treedep *dep )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	int *hist = NULL;
	Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar;
	int nmemim, nmemjm;
	char **tree;
	char *treetmp;
	char *nametmp, *nameptr, *tmpptr; 
	char namec;
	FILE *fp;
	int node[2];
	float *height;

	fp = fopen( "_guidetree", "r" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( nseq );
		ac = (Bchain *)malloc( nseq * sizeof( Bchain ) );
		nmemar = AllocateIntVec( nseq );
//		treetmp = AllocateCharVec( nseq*50 );
		treetmp = NULL;
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( nseq, nseq*50 );
		tree = AllocateCharMtx( nseq, 0 );
		if( dep ) height = AllocateFloatVec( nseq );
	}

	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}


	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;


	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );
		len[k][0] = len[k][1] = -1.0;
		loadtreeoneline( node, len[k], fp );
		im = node[0];
		jm = node[1];

		if( im > nseq-1 || jm > nseq-1 || tree[im] == NULL || tree[jm] == NULL )
		{
			reporterr(       "\n\nCheck the guide tree.\n" );
			reporterr(       "im=%d, jm=%d\n", im+1, jm+1 );
			reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
			exit( 1 );
		}


		if( len[k][0] == -1.0 || len[k][1] == -1.0 )
		{
			reporterr(       "\n\nERROR: Branch length is not given.\n" );
			exit( 1 );
		}

		if( len[k][0] < 0.0 ) len[k][0] = 0.0;
		if( len[k][1] < 0.0 ) len[k][1] = 0.0;


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];

//		reporterr(       "prevnode = %d, nmemim = %d\n", prevnode, nmemim );

		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


		nmemjm = nmemar[jm];
		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;

//		reporterr(       "prevnode = %d, nmemjm = %d\n", prevnode, nmemjm );

		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


//		len[k][0] = ( minscore - tmptmplen[im] );
//		len[k][1] = ( minscore - tmptmplen[jm] );
//		len[k][0] = -1;
//		len[k][1] = -1;


		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

//		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
            }
        }


		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

//		reporterr(       "im,jm=%d,%d\n", im, jm );
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;

		if( dep ) 
		{
			height[im] += len[k][0]; // for ig tree, 2015/Dec/25
			dep[k].distfromtip = height[im]; // for ig tree, 2015/Dec/25
//			reporterr(       "##### dep[%d].distfromtip = %f\n\n", k, height[im] );
		}
    }
	fclose( fp );
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
		fprintf( fp, "#by loadtree\n" );
	fclose( fp );

	FreeCharMtx( tree );
	free( treetmp );
	free( nametmp );
	free( hist );
	free( (char *)ac );
	free( (void *)nmemar );
	if( dep ) free( height );

}

static float sueff1, sueff05;
static double sueff1_double, sueff05_double;

static float cluster_mix_float( float d1, float d2 )
{
	return( MIN( d1, d2 ) * sueff1 + ( d1 + d2 ) * sueff05 ); 
}
static float cluster_average_float( float d1, float d2 )
{
	return( ( d1 + d2 ) * 0.5 ); 
}
static float cluster_minimum_float( float d1, float d2 )
{
	return( MIN( d1, d2 ) ); 
}
static double cluster_mix_double( double d1, double d2 )
{
	return( MIN( d1, d2 ) * sueff1_double + ( d1 + d2 ) * sueff05_double ); 
}
static double cluster_average_double( double d1, double d2 )
{
	return( ( d1 + d2 ) * 0.5 ); 
}
static double cluster_minimum_double( double d1, double d2 )
{
	return( MIN( d1, d2 ) ); 
}


void fixed_supg_float_realloc_nobk_halfmtx_treeout_constrained( int nseq, float **eff, int ***topol, float **len, char **name, int *nlen, Treedep *dep, int ngroup, int **groups, int efffree )
{
	CALLS && printf("called %s:fixed_supg_float_realloc_nobk_halfmtx_treeout_constrained()\n", __FILE__);
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	float tmpfloat;
	float eff1, eff0;
	float *tmptmplen = NULL; //static?
	int *hist = NULL; //static?
	Bchain *ac = NULL; //static?
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti, *acptj;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar; //static?
	int nmemim, nmemjm;
	float minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	float *mindisfrom = NULL; // by D.Mathog, a guess
	char **tree; //static?
	char *treetmp; //static?
	char *nametmp, *nameptr, *tmpptr; //static?
	FILE *fp;
	float (*clusterfuncpt[1])(float,float);
	char namec;
	int *testtopol, **inconsistent;
	int **inconsistentpairlist;
	int ninconsistentpairs;
	int *warned;
	int allinconsistent;
	int firsttime;


	sueff1 = 1 - (float)sueff_global;
	sueff05 = (float)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_float;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_float;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_float;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
		testtopol = AllocateIntVec( njob + 1 );
		inconsistent = AllocateIntMtx( njob, njob ); // muda
		inconsistentpairlist = AllocateIntMtx( njob*(njob-1)/2+1, 2 ); // muda
		warned = AllocateIntVec( ngroup );
	}

	
	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	ninconsistentpairs = 0;
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		for( i=0; i<ninconsistentpairs; i++ ) inconsistent[inconsistentpairlist[i][0]][inconsistentpairlist[i][1]] = 0;
//		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next ) inconsistent[acpti->pos][acptj->pos] = 0; // osoi!!!
		ninconsistentpairs = 0;
		firsttime = 1;
		while( 1 )
		{
			if( firsttime )
			{
				firsttime = 0;
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					if( mindisfrom[i] < minscore ) // muscle
					{
						im = i;
						minscore = mindisfrom[i];
					}
				}
				jm = nearest[im];
				if( jm < im ) 
				{
					j=jm; jm=im; im=j;
				}
			}
			else
			{
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[i][j] && (tmpfloat=eff[i][j-i]) < minscore )
						{
							minscore = tmpfloat;
							im = i; jm = j;
						}
					}
					for( acptj=ac; (acptj&&acptj->pos!=i); acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[j][i] && (tmpfloat=eff[j][i-j]) < minscore )
						{
							minscore = tmpfloat;
							im = j; jm = i;
						}
					}
				}
			}


			allinconsistent = 1;
			for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
			{
				for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
				{
					if( inconsistent[acpti->pos][acptj->pos] == 0 )
					{
						allinconsistent = 0;
						goto exitloop_f;
					}
				}
			}
			exitloop_f:

			if( allinconsistent )
			{
				reporterr(       "\n\n\nPlease check whether the grouping is possible.\n\n\n" );
				exit( 1 );
			}
#if 1
			intpt = testtopol;
			prevnode = hist[im];
			if( prevnode == -1 )
			{
				*intpt++ = im;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			
			prevnode = hist[jm];
			if( prevnode == -1 )
			{
				*intpt++ = jm;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			*intpt = -1;
//			reporterr(       "testtopol = \n" );
//       	for( i=0; testtopol[i]>-1; i++ ) reporterr(       " %03d", testtopol[i]+1 );
//			reporterr(       "\n" );
#endif
			for( i=0; i<ngroup; i++ )
			{
//				reporterr(       "groups[%d] = \n", i );
//  		     	for( j=0; groups[i][j]>-1; j++ ) reporterr(       " %03d", groups[i][j]+1 );
//				reporterr(       "\n" );
				if( overlapmember( groups[i], testtopol ) )
				{
					if( !includemember( testtopol, groups[i] ) && !includemember( groups[i], testtopol ) )
					{
						if( !warned[i] )
						{
							warned[i] = 1;
							reporterr(       "\n###################################################################\n" );
							reporterr(       "# WARNING: Group %d is forced to be a monophyletic cluster.\n", i+1 );
							reporterr(       "###################################################################\n" );
						}
						inconsistent[im][jm] = 1;
						inconsistentpairlist[ninconsistentpairs][0] = im;
						inconsistentpairlist[ninconsistentpairs][1] = jm;
						ninconsistentpairs++;
						break;
					}
				}
			}
			if( i == ngroup )
			{
//				reporterr(       "OK\n" );
				break;
			}
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );
		if( len[k][0] < 0.0 ) len[k][0] = 0.0;
		if( len[k][1] < 0.0 ) len[k][1] = 0.0;

		if( dep ) dep[k].distfromtip = minscore;
//		reporterr(       "\n##### dep[%d].distfromtip = %f\n", k, minscore );

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		eff[im][jm-im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
                tmpfloat = eff[miniim][maxiim-miniim] =
				(clusterfuncpt[0])( eff0, eff1 );
				if( tmpfloat < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpfloat;
					nearest[i] = im;
				}
				if( tmpfloat < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpfloat;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }

		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL;
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif

    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );

	free( tree[0] );
	free( tree );
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
	free( testtopol );
	FreeIntMtx( inconsistent );
	FreeIntMtx( inconsistentpairlist );
	free( warned );
}

void fixed_musclesupg_float_realloc_nobk_halfmtx_treeout( int nseq, float **eff, int ***topol, float **len, char **name, int *nlen, Treedep *dep, int efffree )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	float tmpfloat;
	float eff1, eff0;
	float *tmptmplen = NULL; //static?
	int *hist = NULL; //static?
	Bchain *ac = NULL; //static?
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar; //static?
	int nmemim, nmemjm;
	float minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	float *mindisfrom = NULL; // by D.Mathog, a guess
	char **tree; //static?
	char *treetmp; //static?
	char *nametmp, *nameptr, *tmpptr; //static?
	FILE *fp;
	float (*clusterfuncpt[1])(float,float);
	char namec;


	sueff1 = 1 - (float)sueff_global;
	sueff05 = (float)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_float;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_float;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_float;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
	}

	
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );

		if( dep ) dep[k].distfromtip = minscore;
//		reporterr(       "\n##### dep[%d].distfromtip = %f\n", k, minscore );

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
                tmpfloat = eff[miniim][maxiim-miniim] =
				(clusterfuncpt[0])( eff0, eff1 );
				if( tmpfloat < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpfloat;
					nearest[i] = im;
				}
				if( tmpfloat < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpfloat;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }

		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL;
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif

    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );

	free( tree[0] );
	free( tree );
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}

void fixed_musclesupg_double_treeout( int nseq, double **eff, int ***topol, double **len, char **name )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpfloat;
	double eff1, eff0;
	static double *tmptmplen = NULL;
	static int *hist = NULL;
	static Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	static int *nmemar;
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	double *mindisfrom = NULL; // by D.Mathog, a guess
	static char **tree;
	static char *treetmp;
	static char *nametmp, *nameptr, *tmpptr;
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;


	sueff1_double = 1.0 - sueff_global;
	sueff05_double = sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}


	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateDoubleVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateDoubleVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
	}

	
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}


	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		nmemim = nmemar[im];
//		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		intpt = topol[k][0];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		nmemjm = nmemar[jm];
//		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		intpt = topol[k][1];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );


		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                tmpfloat = eff[miniim][maxiim] =
				(clusterfuncpt[0])( eff0, eff1 );
				if( tmpfloat < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpfloat;
					nearest[i] = im;
				}
				if( tmpfloat < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpfloat;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }
		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim] > mindisfrom[i] )
					setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif

    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );
	free( tree[0] );
	free( tree );
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}

void fixed_supg_double_treeout_constrained( int nseq, double **eff, int ***topol, double **len, char **name, int ngroup, int **groups )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpfloat;
	double eff1, eff0;
	static double *tmptmplen = NULL;
	static int *hist = NULL;
	static Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti, *acptj;
	int *pt1, *pt2, *pt11, *pt22;
	static int *nmemar;
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	double *mindisfrom = NULL; // by D.Mathog, a guess
	static char **tree;
	static char *treetmp;
	static char *nametmp, *nameptr, *tmpptr;
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;
	int *testtopol, **inconsistent;
	int **inconsistentpairlist;
	int ninconsistentpairs;
	int *warned;
	int allinconsistent;
	int firsttime;


	sueff1_double = 1 - sueff_global;
	sueff05_double = sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}


	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateDoubleVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateDoubleVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
		testtopol = AllocateIntVec( njob + 1 );
		inconsistent = AllocateIntMtx( njob, njob ); // muda
		inconsistentpairlist = AllocateIntMtx( njob*(njob-1)/2+1, 2 ); // muda
		warned = AllocateIntVec( ngroup );
	}

	
	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}


	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	ninconsistentpairs = 0;
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );



//		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next ) inconsistent[acpti->pos][acptj->pos] = 0;
		for( i=0; i<ninconsistentpairs; i++ ) inconsistent[inconsistentpairlist[i][0]][inconsistentpairlist[i][1]] = 0;
		ninconsistentpairs = 0;
		firsttime = 1;
		while( 1 )
		{
			if( firsttime )
			{
				firsttime = 0;
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					if( mindisfrom[i] < minscore ) // muscle
					{
						im = i;
						minscore = mindisfrom[i];
					}
				}
				jm = nearest[im];
				if( jm < im ) 
				{
					j=jm; jm=im; im=j;
				}
			}
			else
			{
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[i][j] && (tmpfloat=eff[i][j]) < minscore )
						{
							minscore = tmpfloat;
							im = i; jm = j;
						}
					}
					for( acptj=ac; (acptj&&acptj->pos!=i); acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[j][i] && (tmpfloat=eff[j][i]) < minscore )
						{
							minscore = tmpfloat;
							im = j; jm = i;
						}
					}
				}
			}

			allinconsistent = 1;
			for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
			{
				for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
				{
					if( inconsistent[acpti->pos][acptj->pos] == 0 )
					{
						allinconsistent = 0;
						goto exitloop_d;
					}
				}
			}
			exitloop_d:

			if( allinconsistent )
			{
				reporterr(       "\n\n\nPlease check whether the grouping is possible.\n\n\n" );
				exit( 1 );
			}
#if 1
			intpt = testtopol;
			prevnode = hist[im];
			if( prevnode == -1 )
			{
				*intpt++ = im;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			
			prevnode = hist[jm];
			if( prevnode == -1 )
			{
				*intpt++ = jm;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			*intpt = -1;
//			reporterr(       "testtopol = \n" );
//       	for( i=0; testtopol[i]>-1; i++ ) reporterr(       " %03d", testtopol[i]+1 );
//			reporterr(       "\n" );
#endif
			for( i=0; i<ngroup; i++ )
			{
//				reporterr(       "groups[%d] = \n", i );
//  		     	for( j=0; groups[i][j]>-1; j++ ) reporterr(       " %03d", groups[i][j]+1 );
//				reporterr(       "\n" );
				if( overlapmember( testtopol, groups[i] ) )
				{
					if( !includemember( testtopol, groups[i] ) && !includemember( groups[i], testtopol ) )
					{
						if( !warned[i] )
						{
							warned[i] = 1;
							reporterr(       "\n###################################################################\n" );
							reporterr(       "# WARNING: Group %d is forced to be a monophyletic cluster.\n", i+1 );
							reporterr(       "###################################################################\n" );
						}
						inconsistent[im][jm] = 1;
						inconsistentpairlist[ninconsistentpairs][0] = im;
						inconsistentpairlist[ninconsistentpairs][1] = jm;
						ninconsistentpairs++;
						break;
					}
				}
			}
			if( i == ngroup )
			{
//				reporterr(       "OK\n" );
				break;
			}
		}






		prevnode = hist[im];
		nmemim = nmemar[im];
//		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		intpt = topol[k][0];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		nmemjm = nmemar[jm];
//		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		intpt = topol[k][1];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );
		if( len[k][0] < 0.0 ) len[k][0] = 0.0;
		if( len[k][1] < 0.0 ) len[k][1] = 0.0;


		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		eff[im][jm] = 999.9;
//		eff[im][jm-im] = 999.9; // bug??

		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                tmpfloat = eff[miniim][maxiim] =
				(clusterfuncpt[0])( eff0, eff1 );

				if( tmpfloat < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpfloat;
					nearest[i] = im;
				}
				if( tmpfloat < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpfloat;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }
		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim] > mindisfrom[i] )
					setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );
	free( tree[0] );
	free( tree );
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
	free( testtopol );
	FreeIntMtx( inconsistent );
	FreeIntMtx( inconsistentpairlist );
	free( warned );
}

void fixed_musclesupg_float_realloc_nobk_halfmtx( int nseq, float **eff, int ***topol, float **len, Treedep *dep, int progressout, int efffree )
{
	CALLS && printf("called %s:fixed_musclesupg_float_realloc_nobk_halfmtx()\n", __FILE__);
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	float tmpfloat;
	float eff1, eff0;
	float *tmptmplen = NULL; // static TLS -> local, 2012/02/25
	int *hist = NULL; // static TLS -> local, 2012/02/25 
	Bchain *ac = NULL; // static TLS -> local, 2012/02/25
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar; // static TLS -> local, 2012/02/25
	int nmemim, nmemjm;
	float minscore;
	int *nearest = NULL; // by Mathog, a guess
	float *mindisfrom = NULL; // by Mathog, a guess
	float (*clusterfuncpt[1])(float,float);


	sueff1 = 1 - (float)sueff_global;
	sueff05 = (float)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_float;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_float;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_float;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
	}

	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	if( progressout ) reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( progressout && k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );

		if( dep ) dep[k].distfromtip = minscore;

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
				tmpfloat = eff[miniim][maxiim-miniim] =
				(clusterfuncpt[0])( eff0, eff1 );
				if( tmpfloat < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpfloat;
					nearest[i] = im;
				}
				if( tmpfloat < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpfloat;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }

//		reporterr(       "im,jm=%d,%d\n", im, jm );
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL;
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif

    }
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}









void veryfastsupg_double_loadtree( int nseq, double **eff, int ***topol, double **len, char **name )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double eff1, eff0;
    int *hist = NULL;
	Achain *ac = NULL;
	double minscore;
	char **tree;
	char *treetmp;
	int im = -1, jm = -1;
	int prevnode, acjmnext, acjmprev;
	int *pt1, *pt2, *pt11, *pt22;
	FILE *fp;
	int node[2];
	float lenfl[2];
	char *nametmp, *nameptr, *tmpptr; //static?
	char namec;

	fp = fopen( "_guidetree", "r" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}


	if( !hist )
	{
//		treetmp = AllocateCharVec( njob*50 );
		treetmp = NULL;
//		tree = AllocateCharMtx( njob, njob*50 );
		tree = AllocateCharMtx( njob, 0 );
		nametmp = AllocateCharVec( 1000 ); // nagasugi
		hist = AllocateIntVec( njob );
		ac = (Achain *)malloc( njob * sizeof( Achain ) );
	}

	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}
	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

    for( i=0; i<nseq; i++ ) hist[i] = -1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		lenfl[0] = lenfl[1] = -1.0;
		loadtreeoneline( node, lenfl, fp );
		im = node[0];
		jm = node[1];
		minscore = eff[im][jm];

		if( im > nseq-1 || jm > nseq-1 || tree[im] == NULL || tree[jm] == NULL )
		{
			reporterr(       "\n\nCheck the guide tree.\n" );
			reporterr(       "im=%d, jm=%d\n", im+1, jm+1 );
			reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
			exit( 1 );
		}


//		reporterr(       "im=%d, jm=%d, minscore = %f\n", im, jm, minscore );


		if( lenfl[0] == -1.0 || lenfl[1] == -1.0 )
		{
			reporterr(       "\n\nWARNING: Branch length is not given.\n" );
			exit( 1 );
		}

		if( lenfl[0] < 0.0 ) lenfl[0] = 0.0;
		if( lenfl[1] < 0.0 ) lenfl[1] = 0.0;

//		reporterr(       "im=%d, jm=%d\n", im, jm );

		intpt = topol[k][0];
		prevnode = hist[im];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		intpt = topol[k][1];
		prevnode = hist[jm];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = lenfl[0];
		len[k][1] = lenfl[1];


		hist[im] = k;

		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) +
				( eff0 + eff1 ) * 0.5 * sueff_global;
            }
        }
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		ac[acjmprev].next = acjmnext;
		if( acjmnext != -1 )
			ac[acjmnext].prev = acjmprev;


		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

//		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
//		strcpy( tree[im], treetmp );

    }
	fclose( fp );


	fp = fopen( "infile.tree", "w" );
	fprintf( fp, "%s\n", treetmp );
//	fprintf( fp, "by veryfastsupg_double_loadtree\n" );
	fclose( fp );

	reporterr(       "\n" );
	free( hist );
	free( (char *)ac );
	FreeCharMtx( tree );
	free( treetmp );
	free( nametmp );

}


void veryfastsupg_double_outtree( int nseq, double **eff, int ***topol, double **len, char **name ) // not used
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	static double *tmptmplen = NULL;
    static int *hist = NULL;
	static Achain *ac = NULL;
	double minscore;
	static char **tree;
	static char *treetmp;
	static char *nametmp;
	FILE *fpout;
	int im = -1, jm = -1;
	int prevnode, acjmnext, acjmprev;
	int *pt1, *pt2, *pt11, *pt22;
	double (*clusterfuncpt[1])(double,double);


	sueff1_double = 1 - sueff_global;
	sueff05_double = sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		treetmp = AllocateCharVec( njob*50 );
		tree = AllocateCharMtx( njob, njob*50 );
		hist = AllocateIntVec( njob );
		tmptmplen = (double *)malloc( njob * sizeof( double ) );
		ac = (Achain *)malloc( njob * sizeof( Achain ) );
		nametmp = AllocateCharVec( 31 );
	}

//	for( i=0; i<nseq; i++ ) sprintf( tree[i], "%d", i+1 );
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<30; j++ ) nametmp[j] = 0;
		for( j=0; j<30; j++ ) 
		{
			if( isalnum( name[i][j] ) )
				nametmp[j] = name[i][j];
			else
				nametmp[j] = '_';
		}
		nametmp[30] = 0;
		sprintf( tree[i], "%d_%.20s", i+1, nametmp+1 );
	}
	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) hist[i] = -1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		minscore = 99999.9;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
	        {
				tmpdouble = eff[i][j];
				if( tmpdouble < minscore )
				{
					minscore = tmpdouble;
					im = i; jm = j;
				}
			}
		}

//		reporterr(       "im=%d, jm=%d\n", im, jm );

		intpt = topol[k][0];
		prevnode = hist[im];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		intpt = topol[k][1];
		prevnode = hist[jm];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = minscore - tmptmplen[im];
		len[k][1] = minscore - tmptmplen[jm];

		tmptmplen[im] = minscore;

		hist[im] = k;

		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
				(clusterfuncpt[0])( eff0, eff1 );
            }
        }
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		ac[acjmprev].next = acjmnext;
		if( acjmnext != -1 )
			ac[acjmnext].prev = acjmprev;

		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		strcpy( tree[im], treetmp );
    }
	fpout = fopen( "infile.tree", "w" );
	fprintf( fpout, "%s\n", treetmp );
//	fprintf( fpout, "by veryfastsupg_double_outtree\n" );
	fclose( fpout );
	reporterr(       "\n" );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	FreeCharMtx( tree );
	free( treetmp );
	free( nametmp );
}

void veryfastsupg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	int tmpint;
	int eff1, eff0;
	static double *tmptmplen = NULL;
	static int **eff = NULL;
    static int *hist = NULL;
	static Achain *ac = NULL;
	int minscore;
	double minscoref;
	int im = -1, jm = -1;
	int prevnode, acjmnext, acjmprev;
	int *pt1, *pt2, *pt11, *pt22;
	if( !eff )
	{
		eff = AllocateIntMtx( njob, njob );
		hist = AllocateIntVec( njob );
		tmptmplen = (double *)malloc( njob * sizeof( double ) );
		ac = (Achain *)malloc( njob * sizeof( Achain ) );
	}
	
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			eff[i][j] = (int)( oeff[i][j] * INTMTXSCALE + 0.5 );
		}
	}

	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) hist[i] = -1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		minscore = INTMTXSCALE*4;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
	        {
				tmpint = eff[i][j];
				if( tmpint < minscore )
				{
					minscore = tmpint;
					im = i; jm = j;
				}
			}
		}
		minscoref = (double)minscore * 0.5 / ( INTMTXSCALE );

//		reporterr(       "im=%d, jm=%d\n", im, jm );

#if 1
		intpt = topol[k][0];
		prevnode = hist[im];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		intpt = topol[k][1];
		prevnode = hist[jm];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}
#else
		intpt = topol[k][0];
        for( i=0; i<nseq; i++ )
            if( pair[im][i] > -2 )
				*intpt++ = i;
		*intpt = -1;

		intpt = topol[k][1];
        for( i=0; i<nseq; i++ )
            if( pair[jm][i] > -2 )
				*intpt++ = i;
		*intpt = -1;
#endif

		len[k][0] = minscoref - tmptmplen[im];
		len[k][1] = minscoref - tmptmplen[jm];

		tmptmplen[im] = minscoref;

		hist[im] = k;

		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) + // int??
				( eff0 + eff1 ) * 0.5 * sueff_global;        // int??
            }
        }
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		ac[acjmprev].next = acjmnext;
		if( acjmnext != -1 )
			ac[acjmnext].prev = acjmprev;
    }
	FreeIntMtx( eff ); eff = NULL;
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
}

void fastsupg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	static float *tmplen;
	int *intpt;
	float tmpfloat;
	float eff1, eff0;
	static float **eff = NULL;
    static char **pair = NULL;
	static Achain *ac;
	float minscore;
	int im = -1, jm = -1;
	if( !eff )
	{
		eff = AllocateFloatMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
		tmplen = AllocateFloatVec( njob );
		ac = (Achain *)calloc( njob, sizeof( Achain ) );
	}
	
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			eff[i][j] = (float)oeff[i][j];
		}
	}

	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		minscore = 9999.0;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
//		for( i=0; i<nseq-1; i++ ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
//			for( j=i+1; j<nseq; j++ ) 
	        {
				tmpfloat = eff[i][j];
				if( tmpfloat < minscore )
				{
					minscore = tmpfloat;
					im = i; jm = j;
				}
			}
		}

//		reporterr(       "im=%d, jm=%d\n", im, jm );

		intpt = topol[k][0];
        for( i=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		intpt = topol[k][1];
        for( i=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		minscore /= 2.0;

		len[k][0] = (double)minscore - tmplen[im];
		len[k][1] = (double)minscore - tmplen[jm];

		tmplen[im] = (double)minscore;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

//		for( i=0; i<nseq; i++ )
		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) +
				( eff0 + eff1 ) * 0.5 * sueff_global;
//        		eff[minijm][maxijm] = 9999.0;
            }
        }
		ac[ac[jm].prev].next = ac[jm].next;
		ac[ac[jm].next].prev = ac[jm].prev;
//		eff[im][jm] = 9999.0;
    }
	reporterr(       "\n" );

//	FreeFloatMtx( eff );
//	FreeCharMtx( pair );
//	FreeFloatVec( tmplen );
//	free( ac );
}

void supg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	static float *tmplen;
	int *intpt;
	float **floatptpt;
	float *floatpt;
	float tmpfloat;
	float eff1, eff0;
	static float **eff = NULL;
    static char **pair = NULL;
	if( !eff )
	{
		eff = AllocateFloatMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
		tmplen = AllocateFloatVec( njob );
	}

	
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			eff[i][j] = (float)oeff[i][j];
		}
	}
	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        float minscore = 9999.0;
        int im = -1, jm = -1;


		floatptpt = eff;
        for( i=0; i<nseq-1; i++ ) 
		{
			floatpt = *floatptpt++ + i + 1;
			for( j=i+1; j<nseq; j++ )
	        {
				tmpfloat = *floatpt++;
				if( tmpfloat < minscore )
				{
					minscore = tmpfloat;
					im = i; jm = j;
				}
			}
		}
		intpt = topol[k][0];
        for( i=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		intpt = topol[k][1];
        for( i=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		len[k][0] = (double)minscore / 2.0 - tmplen[im];
		len[k][1] = (double)minscore / 2.0 - tmplen[jm];

		tmplen[im] = (double)minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
#if 1
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
#else
				miniim = MIN( i, im );
				maxiim = MAX( i, im );
				minijm = MIN( i, jm );
				maxijm = MAX( i, jm );
#endif
#if 1
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) +
				( eff0 + eff1 ) * 0.5 * sueff_global;
#else
                MIN( eff[miniim][maxiim], eff[minijm][maxijm] ) * ( 1.0 - sueff_global ) +
				( eff[miniim][maxiim] + eff[minijm][maxijm] ) * 0.5 * sueff_global;
#endif
                eff[minijm][maxijm] = 9999.0;
            	eff[im][jm] = 9999.0;
            }
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

void spg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k;
	double tmplen[M];
	double **eff = NULL;
    char **pair = NULL;
	if( !eff )
	{
		eff = AllocateDoubleMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
	}
	
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) eff[i][j] = oeff[i][j];
	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        float minscore = 9999.0;
        int im = -1, jm = -1;
        int count;

        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
        {
            if( eff[i][j] < minscore )
            {
                minscore = eff[i][j];
                im = i; jm = j;
            }
        }
        for( i=0, count=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
            {
                topol[k][0][count] = i;
                count++;
            }
        topol[k][0][count] = -1;
        for( i=0, count=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
            {
                topol[k][1][count] = i;
                count++;
            }
        topol[k][1][count] = -1;

		len[k][0] = minscore / 2.0 - tmplen[im];
		len[k][1] = minscore / 2.0 - tmplen[jm];

		tmplen[im] = minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
                eff[MIN(i,im)][MAX(i,im)] =
                MIN( eff[MIN(i,im)][MAX(i,im)], eff[MIN(i,jm)][MAX(i,jm)] );
                eff[MIN(i,jm)][MAX(i,jm)] = 9999.0;
            }
            eff[im][jm] = 9999.0;
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

double ipower( double x, int n )    /* n > 0  */
{
    double r;

    r = 1;
    while( n != 0 )
    {
        if( n & 1 ) r *= x;
        x *= x; n >>= 1;
    }
    return( r );
}

void countnode( int nseq, int ***topol, double **node ) /* node[j][i] != node[i][j] */
{
    int i, j, k, s1, s2;
    static double rootnode[M];

    if( nseq-2 < 0 )
	{
		reporterr(       "Too few sequence for countnode: nseq = %d\n", nseq );
		exit( 1 );
    }

    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    for( i=0; i<nseq-2; i++ )
    {
        for( j=0; topol[i][0][j]>-1; j++ )
            rootnode[topol[i][0][j]]++;
        for( j=0; topol[i][1][j]>-1; j++ )
            rootnode[topol[i][1][j]]++;
        for( j=0; topol[i][0][j]>-1; j++ )
        {
            s1 = topol[i][0][j];
            for( k=0; topol[i][1][k]>-1; k++ )
            {
                s2 = topol[i][1][k];
                node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            }
        }
    }
    for( j=0; topol[nseq-2][0][j]>-1; j++ )
    {
        s1 = topol[nseq-2][0][j];
        for( k=0; topol[nseq-2][1][k]>-1; k++ )
        {
            s2 = topol[nseq-2][1][k];
            node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        }
    }
}

void countnode_int( int nseq, int ***topol, int **node )  /* node[i][j] == node[j][i] */
{
    int i, j, k, s1, s2;
    int rootnode[M];

    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    for( i=0; i<nseq-2; i++ )
    {
        for( j=0; topol[i][0][j]>-1; j++ )
            rootnode[topol[i][0][j]]++;
        for( j=0; topol[i][1][j]>-1; j++ )
            rootnode[topol[i][1][j]]++;
        for( j=0; topol[i][0][j]>-1; j++ )
        {
            s1 = topol[i][0][j];
            for( k=0; topol[i][1][k]>-1; k++ )
            {
                s2 = topol[i][1][k];
                node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            }
        }
    }
    for( j=0; topol[nseq-2][0][j]>-1; j++ )
    {
        s1 = topol[nseq-2][0][j];
        for( k=0; topol[nseq-2][1][k]>-1; k++ )
        {
            s2 = topol[nseq-2][1][k];
            node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        }
    }
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) 
		node[j][i] = node[i][j];
#if DEBUG
	reporterr(       "node[][] in countnode_int" );
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			reporterr(       "%#3d", node[i][j] );
		}
		reporterr(       "\n" );
	}
#endif
}

void counteff_simple_float( int nseq, int ***topol, float **len, double *node )
{
    int i, j, s1, s2;
	double total;
	static double rootnode[M];
	static double eff[M];

#if DEBUG
	for( i=0; i<nseq; i++ ){
		reporterr(       "len0 = %f\n", len[i][0] );
		reporterr(       "len1 = %f\n", len[i][1] );
	}
#endif
    for( i=0; i<nseq; i++ )
	{
		rootnode[i] = 0.0;
		eff[i] = 1.0;
/*
		rootnode[i] = 1.0;
*/
	}
   	for( i=0; i<nseq-1; i++ )
   	{
       	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
		{
           	rootnode[s1] += (double)len[i][0] * eff[s1];
			eff[s1] *= 0.5;
/*
           	rootnode[s1] *= 0.5;
*/
			
		}
       	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
		{
           	rootnode[s2] +=  (double)len[i][1] * eff[s2];
			eff[s2] *= 0.5;
/*
           	rootnode[s2] *= 0.5;
*/
				
		}
	}
	for( i=0; i<nseq; i++ ) 
	{
#if 1 /* 97.9.29 */
		rootnode[i] += GETA3;
#endif
	}
	total = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		total += rootnode[i];
	}
		
	for( i=0; i<nseq; i++ ) 
	{
		node[i] = rootnode[i] / total;
	}

}


void counteff_simple_float_nostatic( int nseq, int ***topol, float **len, double *node )
{
	CALLS && printf("called %s:counteff_simple_float_nostatic()\n", __FILE__);
    int i, j, s1, s2;
	double total;
	double *rootnode;
	double *eff;

	rootnode = AllocateDoubleVec( nseq );
	eff = AllocateDoubleVec( nseq );

	for( i=0; i<nseq; i++ ) // 2014/06/07, fu no eff wo sakeru.
	{
		if( len[i][0] < 0.0 ) 
		{
			reporterr( "WARNING: negative branch length %f, step %d-0\n", len[i][0], i );
			len[i][0] = 0.0;
		}
		if( len[i][1] < 0.0 ) 
		{
			reporterr( "WARNING: negative branch length %f, step %d-1\n", len[i][1], i );
			len[i][1] = 0.0;
		}
	}
#if DEBUG
	for( i=0; i<nseq-1; i++ )
	{
		reporterr( "\nstep %d, group 0\n", i );
		for( j=0; topol[i][0][j]!=-1; j++) reporterr( "%3d ", topol[i][0][j] );
		reporterr( "\n", i );
		reporterr( "step %d, group 1\n", i );
		for( j=0; topol[i][1][j]!=-1; j++) reporterr( "%3d ", topol[i][1][j] );
		reporterr( "\n", i );
		reporterr(       "len0 = %f\n", len[i][0] );
		reporterr(       "len1 = %f\n", len[i][1] );
	}
#endif
    for( i=0; i<nseq; i++ )
	{
		rootnode[i] = 0.0;
		eff[i] = 1.0;
/*
		rootnode[i] = 1.0;
*/
	}
   	for( i=0; i<nseq-1; i++ )
   	{
       	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
		{
           	rootnode[s1] += (double)len[i][0] * eff[s1];
			eff[s1] *= 0.5;
/*
           	rootnode[s1] *= 0.5;
*/
			
		}
       	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
		{
           	rootnode[s2] +=  (double)len[i][1] * eff[s2];
			eff[s2] *= 0.5;
/*
           	rootnode[s2] *= 0.5;
*/
				
		}
	}
	for( i=0; i<nseq; i++ ) 
	{
#if 1 /* 97.9.29 */
		rootnode[i] += GETA3;
#endif
	}
#if 1
	total = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		total += rootnode[i];
	}
#else
	total = 1.0;
#endif
		
	for( i=0; i<nseq; i++ ) 
	{
		node[i] = rootnode[i] / total;
	}

	free( rootnode );
	free( eff );
}

void counteff_simple( int nseq, int ***topol, double **len, double *node )
{
    int i, j, s1, s2;
	double total;
	double *rootnode;
	double *eff;
	rootnode = AllocateDoubleVec( nseq );
	eff = AllocateDoubleVec( nseq );

#if DEBUG
	for( i=0; i<nseq; i++ ){
		reporterr(       "len0 = %f\n", len[i][0] );
		reporterr(       "len1 = %f\n", len[i][1] );
	}
#endif
    for( i=0; i<nseq; i++ )
	{
		rootnode[i] = 0.0;
		eff[i] = 1.0;
/*
		rootnode[i] = 1.0;
*/
	}
   	for( i=0; i<nseq-1; i++ )
   	{
       	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
		{
           	rootnode[s1] += len[i][0] * eff[s1];
			eff[s1] *= 0.5;
/*
           	rootnode[s1] *= 0.5;
*/
			
		}
       	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
		{
           	rootnode[s2] +=  len[i][1] * eff[s2];
			eff[s2] *= 0.5;
/*
           	rootnode[s2] *= 0.5;
*/
				
		}
	}
	for( i=0; i<nseq; i++ ) 
	{
#if 1 /* 97.9.29 */
		rootnode[i] += GETA3;
#endif
	}
#if 1
	total = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		total += rootnode[i];
	}
#else
	total = 1.0;
#endif
		
	for( i=0; i<nseq; i++ ) 
	{
		node[i] = rootnode[i] / total;
	}

#if 1
	free( rootnode );
	free( eff );
#endif
}


void counteff( int nseq, int ***topol, double **len, double **node )
{
    int i, j, k, s1, s2;
	double rootnode[M];
	double eff[M];

	if( mix ) 
	{
		switch( weight )
		{
			case( 2 ): 
				weight = 3;
				break;
			case( 3 ): 
				weight = 2;
				break;
			default: 
				ErrorExit( "mix error" );
				break;
		}
	}

	if( weight == 2 )
	{
	    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    	for( i=0; i<nseq-2; i++ )
    	{
        	for( j=0; topol[i][0][j]>-1; j++ )
            	rootnode[topol[i][0][j]]++;
        	for( j=0; topol[i][1][j]>-1; j++ )
            	rootnode[topol[i][1][j]]++;
        	for( j=0; topol[i][0][j]>-1; j++ )
        	{
            	s1 = topol[i][0][j];
            	for( k=0; topol[i][1][k]>-1; k++ )
            	{
                	s2 = topol[i][1][k];
                	node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            	}
        	}
    	}
    	for( j=0; topol[nseq-2][0][j]>-1; j++ )
    	{
        	s1 = topol[nseq-2][0][j];
        	for( k=0; topol[nseq-2][1][k]>-1; k++ )
        	{
            	s2 = topol[nseq-2][1][k];
            	node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        	}
    	}
   		for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
   	   		node[i][j] = ipower( 0.5, (int)node[i][j] ) + geta2;
		for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) 
			node[j][i] = node[i][j];
	}

	if( weight == 3 )
	{
#if DEBUG
		for( i=0; i<nseq; i++ ){
			reporterr(       "len0 = %f\n", len[i][0] );
			reporterr(       "len1 = %f\n", len[i][1] );
		}
#endif
	    for( i=0; i<nseq; i++ )
		{
			rootnode[i] = 0.0;
			eff[i] = 1.0;
/*
			rootnode[i] = 1.0;
*/
		}
    	for( i=0; i<nseq-1; i++ )
    	{
        	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
			{
   	        	rootnode[s1] += len[i][0] * eff[s1];
				eff[s1] *= 0.5;
/*
   	        	rootnode[s1] *= 0.5;
*/
				
			}
   	    	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
			{
   	        	rootnode[s2] +=  len[i][1] * eff[s2];
				eff[s2] *= 0.5;
/*
   	        	rootnode[s2] *= 0.5;
*/
				
			}
		}
		for( i=0; i<nseq; i++ ) 
		{
#if 1 /* 97.9.29 */
			rootnode[i] += GETA3;
#endif
#if DEBUG
			reporterr(       "rootnode for %d = %f\n", i, rootnode[i] );
#endif
		}
		for( i=0; i<nseq; i++ ) 
		{
			for( j=0; j<nseq; j++ ) 
				if( j != i )
					node[i][j] = (double)rootnode[i] * rootnode[j];
				else node[i][i] = rootnode[i];
		}
	}

}

float score_calcp( char *seq1, char *seq2, int len )
{
	int k;
	int ms1, ms2;
	float tmpscore;
	int len2 = len - 2;

	tmpscore = 0.0;
	for( k=0; k<len; k++ )
	{
		ms1 = (int)seq1[k];
		ms2 = (int)seq2[k];
		if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
		tmpscore += (float)amino_dis[ms1][ms2];
	
		if( ms1 == (int)'-' ) 
		{
			tmpscore += (float)penalty;
			tmpscore += (float)amino_dis[ms1][ms2];
			while( (ms1=(int)seq1[++k]) == (int)'-' )
				tmpscore += (float)amino_dis[ms1][ms2];
			k--;
			if( k >len2 ) break;
			continue;
		}
		if( ms2 == (int)'-' )
		{
			tmpscore += (float)penalty;
			tmpscore += (float)amino_dis[ms1][ms2];
			while( (ms2=(int)seq2[++k]) == (int)'-' )
				tmpscore += (float)amino_dis[ms1][ms2];
			k--;
			if( k > len2 ) break;
			continue;
		}
	}
	return( tmpscore );
}

float score_calc1( char *seq1, char *seq2 )   /* method 1 */
{
	int k;
	float score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (float)amino_dis[(int)seq1[k]][(int)seq2[k]];
			count++;
		}
	}
	if( count ) score /= (float)count;
	else score = 1.0;
	return( score );
}

float substitution_nid( char *seq1, char *seq2 )
{
	int k;
	float s12;
	int len = strlen( seq1 );
	
	s12 = 0.0;
	for( k=0; k<len; k++ )
		if( seq1[k] != '-' && seq2[k] != '-' )
			s12 += ( seq1[k] == seq2[k] );

//	fprintf( stdout, "s12 = %f\n", s12 );
	return( s12 );
}

float substitution_score( char *seq1, char *seq2 )
{
	int k;
	float s12;
	int len = strlen( seq1 );
	
	s12 = 0.0;
	for( k=0; k<len; k++ )
		if( seq1[k] != '-' && seq2[k] != '-' )
			s12 += amino_dis[(int)seq1[k]][(int)seq2[k]];

//	fprintf( stdout, "s12 = %f\n", s12 );
	return( s12 );
}

float substitution_hosei( char *seq1, char *seq2 )   /* method 1 */
{
	int count = 0;
	float score;
	int iscore = 0;
	char s1, s2;

	while( (s1=*seq1++) )
	{
		s2 = *seq2++;
		if( s1 == '-' ) continue;
		if( s2 == '-' ) continue;
		iscore += ( s1 != s2 );
		count++;
	}
	if( count ) score = (float)iscore / count;
	else score = 1.0;
	if( score < 0.95 ) score = - log( 1.0 - score );
	else score = 3.0;
	return( score );
}

float substitution( char *seq1, char *seq2 )   /* method 1 */
{
	int k;
	float score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (float)( seq1[k] != seq2[k] );
			count++;
		}
	}
	if( count ) score /= (float)count;
	else score = 1.0;
	return( score );
}


void treeconstruction( char **seq, int nseq, int ***topol, double **len, double **eff )
{
    int i, j;

	if( weight > 1 )
	{
		if( utree == 0 )
		{
	    	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
   		 	{
/*
		       	 eff[i][j] = (double)score_calc1( seq[i], seq[j] );
*/
		       	 eff[i][j] = (double)substitution_hosei( seq[i], seq[j] );
 /*
				 reporterr(       "%f\n", eff[i][j] );
 */
   		 	}
/*
			reporterr(       "distance matrix\n" );
			for( i=0; i<nseq; i++ )
			{
				for( j=0; j<nseq; j++ ) 
				{
					reporterr(       "%f ", eff[i][j] );
				}
				reporterr(       "\n" );
			}
*/
/*
   			upg( nseq, eff, topol, len );
   			upg2( nseq, eff, topol, len );
*/
   			spg( nseq, eff, topol, len );
   			counteff( nseq, topol, len, eff );
		}
	}
	else
	{
		for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) 
			eff[i][j] = 1.0;
	}
/*
reporterr(       "weight matrix\n" );
for( i=0; i<nseq; i++ )
{
	for( j=0; j<nseq; j++ ) 
	{
		reporterr(       "%f ", eff[i][j] );
	}
	reporterr(       "\n" );
}
*/
}

float bscore_calc( char **seq, int s, double **eff )  /* algorithm B */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    long score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 + !gb1  * !gc1 
                 * !gb2  *  gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
				 + gb1   * !gc1
				 * gb2   *  gc2      *BEFF

				 + gb1   *  gc1
				 * gb2   * !gc2      *BEFF
                 ;
			score += (long)cob * penalty * efficient;
			score += (long)amino_dis[(int)seq[i][k]][(int)seq[j][k]] * efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (float)score / nglen + 400.0 * !scoremtx );
}

void AllocateTmpSeqs( char ***mseq2pt, char **mseq1pt, int locnlenmax )
{
	*mseq2pt = AllocateCharMtx( njob, locnlenmax+1 );
	*mseq1pt = AllocateCharVec( locnlenmax+1 );
}

void FreeTmpSeqs( char **mseq2, char *mseq1 )
{
	FreeCharMtx( mseq2 );
	free( (char *)mseq1 );
}


void gappick0( char *aseq, char *seq )
{
	CALLS && printf("called %s:gappick0()\n", __FILE__);
	for( ; *seq != 0; seq++ )
	{
		if( *seq != '-' )
			*aseq++ = *seq;
	}
	*aseq = 0;

}

void gappick( int nseq, int s, char **aseq, char **mseq2, 
			  double **eff, double *effarr )
{
	int i, j, count, countjob, len, allgap;
	len = strlen( aseq[0] );
	for( i=0, count=0; i<len; i++ ) 
	{
		allgap = 1;
		for( j=0; j<nseq; j++ ) if( j != s ) allgap *= ( aseq[j][i] == '-' );
        if( allgap == 0 )
		{
			for( j=0, countjob=0; j<nseq; j++ ) 
			{
				if( j != s )
				{
					mseq2[countjob][count] = aseq[j][i];
					countjob++;
				}
			}
			count++;
		}
	}
	for( i=0; i<nseq-1; i++ ) mseq2[i][count] = 0;

	for( i=0, countjob=0; i<nseq; i++ ) 
	{
		if( i != s )
		{
			effarr[countjob] = eff[s][i];
			countjob++;
		}
	}
/*
fprintf( stdout, "effarr in gappick s = %d\n", s+1 );
for( i=0; i<countjob; i++ ) 
	fprintf( stdout, " %f", effarr[i] );
printf( "\n" );
*/
}

void commongappick_record( int nseq, char **seq, int *map )
{
	int i, j, count;
	int len = strlen( seq[0] );


	for( i=0, count=0; i<=len; i++ ) 
	{
	/*
		allgap = 1;
		for( j=0; j<nseq; j++ ) 
			allgap *= ( seq[j][i] == '-' );
		if( !allgap )
	*/
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			for( j=0; j<nseq; j++ )
			{
				seq[j][count] = seq[j][i];
			}
			map[count] = i;
			count++;
	 	}
	}
}

void commongappick( int nseq, char **seq )
{
	CALLS && printf("called %s:commongappick()\n", __FILE__);
	int i, j, count;
	int len = strlen( seq[0] );
#if 1

	int *mapfromnewtoold;

	mapfromnewtoold = calloc( len+1, sizeof( int ) );

	for( i=0, count=0; i<=len; i++ ) 
	{
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			mapfromnewtoold[count++] = i;
	 	}
	}
//	mapfromnewtoold[count] = -1; // iranai
	for( j=0; j<nseq; j++ )
	{
		for( i=0; i<count; i++ )
		{
			seq[j][i] = seq[j][mapfromnewtoold[i]];
		}
	}
	free( mapfromnewtoold );
#else
--------------------------

	int *mapfromoldtonew;
	int pos;

	mapfromoldtonew = calloc( len+1, sizeof( int ) );
	for( i=0; i<=len; i++ ) mapfromoldtonew[i] = -1;

	for( i=0, count=0; i<=len; i++ ) 
	{
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			mapfromoldtonew[i] = count;
			count++;
	 	}
	}
	for( j=0; j<nseq; j++ )
	{
		for( i=0; i<=len; i++ ) 
		{
			if( (pos=mapfromoldtonew[i]) != -1 )
				seq[j][pos] = seq[j][i];
		}
	}
	free( mapfromoldtonew );
--------------------------

	for( i=0, count=0; i<=len; i++ ) 
	{
	/*
		allgap = 1;
		for( j=0; j<nseq; j++ ) 
			allgap *= ( seq[j][i] == '-' );
		if( !allgap )
	*/
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			for( j=0; j<nseq; j++ )
			{
				seq[j][count] = seq[j][i];
			}
			count++;
	 	}
	}

#endif
}
		
double score_calc0( char **seq, int s, double **eff, int ex )
{
	double tmp;

	if( scmtd == 4 ) tmp = score_calc4( seq, s, eff, ex );
	if( scmtd == 5 ) tmp = score_calc5( seq, s, eff, ex );
	else             tmp = score_calc5( seq, s, eff, ex );

	return( tmp );

}



void strins( char *str1, char *str2 )
{
	char *bk;
	int len1 = strlen( str1 );
	int len2 = strlen( str2 );

	bk = str2;
	str2 += len1+len2;
	str1 += len1-1;

	while( str2 >= bk+len1 ) { *str2 = *(str2-len1); str2--;} // by D.Mathog
	while( str2 >= bk ) { *str2-- = *str1--; }
}

int isaligned( int nseq, char **seq )
{
	int i;
	int len = strlen( seq[0] );
	for( i=1; i<nseq; i++ ) 
	{
		if( strlen( seq[i] ) != len ) return( 0 );
	}
	return( 1 );
}

double score_calc_for_score( int nseq, char **seq )
{
    int i, j, k, c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;

    score = 0.0;
    for( i=0; i<nseq-1; i++ )
    {
        for( j=i+1; j<nseq; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            c = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                c++;
                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq1[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq2[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore / (double)c;
#if DEBUG
			printf( "tmpscore in mltaln9.c = %f\n", tmpscore );
			printf( "tmpscore / c          = %f\n", tmpscore/(double)c );
#endif
        }
    }
	reporterr(       "raw score = %f\n", score );
	score /= (double)nseq * ( nseq-1.0 ) / 2.0;
	score += 400.0;
#if DEBUG
	printf( "score in mltaln9.c = %f\n", score );
#endif
    return( (double)score );
}

void floatncpy( float *vec1, float *vec2, int len )
{
	while( len-- )
		*vec1++ = *vec2++;
}

float score_calc_a( char **seq, int s, double **eff )  /* algorithm A+ */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
                 ;
			score += 0.5 * (float)cob * penalty * efficient;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]] * (float)efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (float)score / nglen + 400.0 * !scoremtx );
}


float score_calc_s( char **seq, int s, double **eff )  /* algorithm S, not used */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
                 ;
			score += 0.5 * (float)cob * penalty * efficient;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]] * (float)efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (float)score / nglen + 400.0 );
}

double score_calc_for_score_s( int s, char **seq )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
                 ;
			score += 0.5 * (float)cob * penalty;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]];
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (double)score / nglen + 400.0 );
}

double SSPscore___( int s, char **seq, int ex )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	i=ex; for( j=0; j<s; j++ )
	{

		if( j == ex ) continue;

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2 * 2.0 

                 +  gb1  * !gc1
                 * !gb2  *  gc2 * 2.0 
      
                 ;
			score += 0.5 * (float)cob * penalty;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]];
			nglen += ( !gc1 * !gc2 ); /* tsukawanai */
		}
	}
	return( (double)score );
}

double SSPscore( int s, char **seq )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    float score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
                 ;
			score += 0.5 * (float)cob * penalty;
			score += (float)amino_dis[(int)seq[i][k]][(int)seq[j][k]];
			nglen += ( !gc1 * !gc2 ); /* tsukawanai */
		}
	}
	return( (double)score );
}

double DSPscore( int s, char **seq )  /* method 3 deha nai */
{
    int i, j, k;
    double c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;
#if DEBUG
	FILE *fp;
#endif

    score = 0.0;
    c = 0.0;

    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(int)mseq1[k]][(int)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore;
        }
    }

	return( score );
}


#define SEGMENTSIZE 150

int searchAnchors( int nseq, char **seq, Segment *seg )
{
	int i, j, k, kcyc;
	int status;
	double score;
	int value = 0;
	int len;
	int length;
	static double *stra = NULL;
	static int alloclen = 0;
	double cumscore;
	static double threshold;

	len = strlen( seq[0] );
	if( alloclen < len )
	{
		if( alloclen )
		{
			FreeDoubleVec( stra );
		}
		else
		{
			threshold = (int)divThreshold / 100.0 * 600.0 * divWinSize;
		}
		stra = AllocateDoubleVec( len );
		alloclen = len;
	}

	for( i=0; i<len; i++ )
	{
		stra[i] = 0.0;
		kcyc = nseq-1;
		for( k=0; k<kcyc; k++ ) for( j=k+1; j<nseq; j++ )
			stra[i] += n_dis[(int)amino_n[(int)seq[k][i]]][(int)amino_n[(int)seq[j][i]]];
		stra[i] /= (double)nseq * ( nseq-1 ) / 2;
	}

	(seg+0)->skipForeward = 0;
	(seg+1)->skipBackward = 0;
	status = 0;
	cumscore = 0.0;
	score = 0.0;
	length = 0; /* modified at 01/09/11 */
	for( j=0; j<divWinSize; j++ ) score += stra[j];
	for( i=1; i<len-divWinSize; i++ )
	{
		score = score - stra[i-1] + stra[i+divWinSize-1];
#if DEBUG
		reporterr(       "%d %f   ? %f", i, score, threshold );
		if( score > threshold ) reporterr(       "YES\n" );
		else                    reporterr(       "NO\n" );
#endif

		if( score > threshold )
		{
			if( !status )
			{
				status = 1;
				seg->start = i;
				length = 0;
				cumscore = 0.0;
			}
			length++;
			cumscore += score;
		}
		if( score <= threshold || length > SEGMENTSIZE )
		{
			if( status )
			{
				seg->end = i;
				seg->center = ( seg->start + seg->end + divWinSize ) / 2 ;
				seg->score = cumscore;
#if DEBUG
				reporterr(       "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
				if( length > SEGMENTSIZE )
				{
					(seg+0)->skipForeward = 1;
					(seg+1)->skipBackward = 1;
				}
				else
				{
					(seg+0)->skipForeward = 0;
					(seg+1)->skipBackward = 0;
				}
				length = 0;
				cumscore = 0.0;
				status = 0;
				value++;
				seg++;
				if( value > MAXSEG - 3 ) ErrorExit( "TOO MANY SEGMENTS!");
			}
		}
	}
	if( status )
	{
		seg->end = i;
		seg->center = ( seg->start + seg->end + divWinSize ) / 2 ;
		seg->score = cumscore;
#if DEBUG
reporterr(       "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
		value++;
	}
	return( value );
}

void dontcalcimportance( int nseq, double *eff, char **seq, LocalHom **localhom )
{
	int i, j;
	LocalHom *ptr;
	int *nogaplen;

	nogaplen = AllocateIntVec( nseq );

	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[%d] = %d\n", i, nogaplen[i] );
	}

	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<nseq; j++ )
		{
			for( ptr=localhom[i]+j; ptr; ptr=ptr->next )
			{
//				reporterr(       "i,j=%d,%d,ptr=%p\n", i, j, ptr );
#if 1
				ptr->importance = ptr->opt / ptr->overlapaa;
				ptr->fimportance = (float)ptr->importance;
#else
				ptr->importance = ptr->opt / MIN( nogaplen[i], nogaplen[j] );
#endif
			}
		}
	}
	free( nogaplen );
}

void dontcalcimportance_firstone( int nseq, double *eff, char **seq, LocalHom **localhom )
{
	int i, j, nseq1;
	LocalHom *ptr;
#if 1
#else
	int *nogaplen;
	nogaplen = AllocateIntVec( nseq );
	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[%d] = %d\n", i, nogaplen[i] );
	}
#endif

	nseq1 = nseq - 1;
	for( i=0; i<nseq1; i++ )
	{
		j=0;
		{
			for( ptr=localhom[i]+j; ptr; ptr=ptr->next )
			{
//				reporterr(       "i,j=%d,%d,ptr=%p\n", i, j, ptr );
#if 1
//				ptr->importance = ptr->opt / ptr->overlapaa;
				ptr->importance = ptr->opt * 0.5; // tekitou
				ptr->fimportance = (float)ptr->importance;
//				reporterr(       "i=%d, j=%d, importance = %f, opt=%f\n", i, j, ptr->fimportance, ptr->opt  );
#else
				ptr->importance = ptr->opt / MIN( nogaplen[i], nogaplen[j] );
#endif
			}
		}
	}
#if 1
#else
	free( nogaplen );
#endif
}

void calcimportance( int nseq, double *eff, char **seq, LocalHom **localhom )
{
	int i, j, pos, len;
	double *importance; // static -> local, 2012/02/25
	double tmpdouble;
	double *ieff, totaleff; // counteff_simple_float ni utsusu kamo
	int *nogaplen; // static -> local, 2012/02/25
	LocalHom *tmpptr;

	importance = AllocateDoubleVec( nlenmax );
	nogaplen = AllocateIntVec( nseq );
	ieff = AllocateDoubleVec( nseq );

	totaleff = 0.0;
	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[] = %d\n", nogaplen[i] );
		if( nogaplen[i] == 0 ) ieff[i] = 0.0;
		else ieff[i] = eff[i];
		totaleff += ieff[i];
	}
	for( i=0; i<nseq; i++ ) ieff[i] /= totaleff;
//	for( i=0; i<nseq; i++ ) reporterr(       "eff[%d] = %f\n", i, ieff[i] );


	for( i=0; i<nseq; i++ )
	{
//		reporterr(       "i = %d\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			importance[pos] = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( i == j ) continue;
			tmpptr = localhom[i]+j;
			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1 ) continue;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
#if 1
					importance[pos] += ieff[j];
#else
					importance[pos] += ieff[j] * tmpptr->opt / MIN( nogaplen[i], nogaplen[j] );
					importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
				}
			}
		}
		for( j=0; j<nseq; j++ )
		{
//			reporterr(       "i=%d, j=%d\n", i, j );
			if( i == j ) continue;
			if( localhom[i][j].opt == -1.0 ) continue;
#if 1
			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				tmpdouble = 0.0;
				len = 0;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
					tmpdouble += importance[pos];
					len++;
				}

				tmpdouble /= (double)len;

				tmpptr->importance = tmpdouble * tmpptr->opt;
				tmpptr->fimportance = (float)tmpptr->importance;
			}
#else
			tmpdouble = 0.0;
			len = 0;
			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
					tmpdouble += importance[pos];
					len++;
				}
			}

			tmpdouble /= (double)len;

			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				tmpptr->importance = tmpdouble * tmpptr->opt;
//				tmpptr->importance = tmpptr->opt / tmpptr->overlapaa; //なかったことにする
			}
#endif

//			reporterr(       "importance of match between %d - %d = %f\n", i, j, tmpdouble );
		}
	}


#if 1
//	reporterr(       "average?\n" );
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
	{
		double imp;
		LocalHom *tmpptr1, *tmpptr2;

//		reporterr(       "i=%d, j=%d\n", i, j );

		tmpptr1 = localhom[i]+j; tmpptr2 = localhom[j]+i;
		for( ; tmpptr1 && tmpptr2; tmpptr1 = tmpptr1->next, tmpptr2 = tmpptr2->next)
		{
			if( tmpptr1->opt == -1.0 || tmpptr2->opt == -1.0 ) 
			{
//				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
				continue;
			}
//			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
			imp = 0.5 * ( tmpptr1->importance + tmpptr2->importance );
			tmpptr1->importance = tmpptr2->importance = imp;
			tmpptr1->fimportance = tmpptr2->fimportance = (float)imp;

//			reporterr(       "## importance = %f\n", tmpptr1->importance );

		}
	}
#endif
	free( importance );
	free( nogaplen );
	free( ieff );
}




static void addlocalhom2_e( LocalHom *pt, LocalHom *lh, int sti, int stj, int eni, int enj, double opt, int overlp, int interm )
{
// dokka machigatteru
	if( pt != lh ) // susumeru
	{
		pt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
		pt = pt->next;
		pt->next = NULL;
		lh->last = pt;
	}
	else // sonomamatsukau
	{
		lh->last = pt;
	}
	lh->nokori++;
//	reporterr(       "in addlocalhom2_e, pt = %p, pt->next = %p, interm=%d, sti-eni-stj-enj=%d %d %d %d\n", pt, pt->next, interm, sti, eni, stj, enj );

	pt->start1 = sti;
	pt->start2 = stj;
	pt->end1 = eni;
	pt->end2 = enj;
	pt->opt = opt;
	pt->extended = interm;
	pt->overlapaa = overlp;
}

void extendlocalhom2( int nseq, LocalHom **localhom, double **dist )
{
	int overlp, plim;
	int i, j, k;
	int pi, pj, pk, len;
	int status, sti, stj;
	int *ipt;
	int co;
	static int *ini = NULL;
	static int *inj = NULL;
	LocalHom *pt;

	sti = 0; // by D.Mathog, a guess
	stj = 0; // by D.Mathog, a guess

	if( ini == NULL )
	{
		ini = AllocateIntVec( nlenmax+1 );
		inj = AllocateIntVec( nlenmax+1 );
	}


	for( i=0; i<nseq-1; i++ )
	{
		for( j=i+1; j<nseq; j++ )
		{
			for( k=0; k<nseq; k++ )
			{
//				reporterr(       "i=%d, j=%d, k=%d, dists = %f,%f,%f thrinter=%f\n", i, j, k, dist[i][j], dist[MIN(i,k)][MAX(i,k)], dist[MIN(j,k)][MAX(j,k)], thrinter );
				if( k == i || k == j ) continue; // mou yatta nomo habuita hoga ii 
				if( dist[MIN(i,k)][MAX(i,k)] > dist[i][j] * thrinter || dist[MIN(j,k)][MAX(j,k)] > dist[i][j] * thrinter ) continue;
				ipt = ini; co = nlenmax+1;
				while( co-- ) *ipt++ = -1;
				ipt = inj; co = nlenmax+1;
				while( co-- ) *ipt++ = -1;
				overlp = 0;

				{
					for( pt=localhom[i]+k; pt; pt=pt->next )
		        	{
//						reporterr(       "i=%d,k=%d,st1:st2=%d:%d,pt=%p,extended=%p\n", i, k, pt->start1, pt->start2, pt, pt->extended );
						if( pt->opt == -1 )
						{
							reporterr(       "opt kainaide tbfast.c = %f\n", pt->opt );
						}
						if( pt->extended > -1 ) break;
						pi = pt->start1;
						pk = pt->start2;
						len = pt->end1 - pt->start1 + 1;
						ipt = ini + pk;
						while( len-- ) *ipt++ = pi++;
					}
				}

				{
					for( pt=localhom[j]+k; pt; pt=pt->next )
		        	{
						if( pt->opt == -1 )
						{
							reporterr(       "opt kainaide tbfast.c = %f\n", pt->opt );
						}
						if( pt->extended > -1 ) break;
						pj = pt->start1;
						pk = pt->start2;
						len = pt->end1 - pt->start1 + 1;
						ipt = inj + pk;
						while( len-- ) *ipt++ = pj++;
					}
				}
				overlp = 0;
				plim = nlenmax+1;
				for( pk = 0; pk < plim; pk++ )
					if( ini[pk] != -1 && inj[pk] != -1 ) overlp++;


				status = 0;
				plim = nlenmax+1;
				for( pk=0; pk<plim; pk++ )
				{
//					reporterr(       "%d %d: %d-%d\n", i, j, ini[pk], inj[pk] );
					if( status )
					{
						if( ini[pk] == -1 || inj[pk] == -1 || ini[pk-1] != ini[pk] - 1 || inj[pk-1] != inj[pk] - 1 ) // saigonoshori
						{
							status = 0;
//							reporterr(       "end here!\n" );

							pt = localhom[i][j].last;
//							reporterr(       "in ex (ba), pt = %p, nokori=%d, i,j,k=%d,%d,%d\n", pt, localhom[i][j].nokori, i, j, k );
							addlocalhom2_e( pt, localhom[i]+j, sti, stj, ini[pk-1], inj[pk-1], MIN( localhom[i][k].opt, localhom[j][k].opt ) * 1.0, overlp, k );
//							reporterr(       "in ex, pt = %p, pt->next = %p, pt->next->next = %p\n", pt, pt->next, pt->next->next );

							pt = localhom[j][i].last;
//							reporterr(       "in ex (ba), pt = %p, pt->next = %p\n", pt, pt->next );
//							reporterr(       "in ex (ba), pt = %p, pt->next = %p, k=%d\n", pt, pt->next, k );
							addlocalhom2_e( pt, localhom[j]+i, stj, sti, inj[pk-1], ini[pk-1], MIN( localhom[i][k].opt, localhom[j][k].opt ) * 1.0, overlp, k );
//							reporterr(       "in ex, pt = %p, pt->next = %p, pt->next->next = %p\n", pt, pt->next, pt->next->next );
						}
					}
					if( !status ) // else deha arimasenn.
					{
						if( ini[pk] == -1 || inj[pk] == -1 ) continue;
						sti = ini[pk];
						stj = inj[pk];
//						reporterr(       "start here!\n" );
						status = 1;
					}
				}
//				if( status ) reporterr(       "end here\n" );

//				exit( 1 );
//					fprintf( hat3p, "%d %d %d %6.3f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->next ); 
			}
		}
	}
}

int makelocal( char *s1, char *s2, int thr )
{
	int start, maxstart, maxend;
	char *pt1, *pt2;
	double score;
	double maxscore;

	pt1 = s1;
	pt2 = s2;

	maxend = 0; // by D.Mathog, a guess

//	reporterr(       "thr = %d, \ns1 = %s\ns2 = %s\n", thr, s1, s2 );
	maxscore = 0.0;
	score = 0.0;
	start = 0;
	maxstart = 0;
	while( *pt1 )
	{
//		reporterr(       "*pt1 = %c*pt2 = %c\n", *pt1, *pt2 );
		if( *pt1 == '-' || *pt2 == '-' )
		{
//			reporterr(       "penalty = %d\n", penalty );
			score += penalty;
			while( *pt1 == '-' || *pt2 == '-' )
			{
				pt1++; pt2++;
			}
			continue;
		}

		score += ( amino_dis[(int)*pt1++][(int)*pt2++] - thr );
//		score += ( amino_dis[(int)*pt1++][(int)*pt2++] );
		if( score > maxscore ) 
		{
//			reporterr(       "score = %f\n", score );
			maxscore = score;
			maxstart = start;
//			reporterr(       "## max! maxstart = %d, start = %d\n", maxstart, start );
		}
		if( score < 0.0 )
		{
//			reporterr(       "## resetting, start = %d, maxstart = %d\n", start, maxstart );
			if( start == maxstart )
			{
				maxend = pt1 - s1;
//				reporterr(       "maxend = %d\n", maxend );
			}
			score = 0.0;
			start = pt1 - s1;
		}
	}
	if( start == maxstart )
		maxend = pt1 - s1 - 1;

//	reporterr(       "maxstart = %d, maxend = %d, maxscore = %f\n", maxstart, maxend, maxscore );
	s1[maxend+1] = 0;
	s2[maxend+1] = 0;
	return( maxstart );
}

void resetlocalhom( int nseq, LocalHom **lh )
{
	int i, j;
	LocalHom *pt;

	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
	{
		for( pt=lh[i]+j; pt; pt=pt->next )
			pt->opt = 1.0;
	}

}

void gapireru( char *res, char *ori, char *gt )
{
	char g;
	while( (g = *gt++) )
	{
		if( g == '-' )
		{
			*res++ = *newgapstr;
		}
		else
		{
			*res++ = *ori++;
		}
	}
	*res = 0;
}

void getkyokaigap( char *g, char **s, int pos, int n )
{
//	char *bk = g;
//	while( n-- ) *g++ = '-';
	while( n-- ) *g++ = (*s++)[pos];

//	reporterr(       "bk = %s\n", bk );
}

void new_OpeningGapCount( float *ogcp, int clus, char **seq, double *eff, int len, char *sgappat )
{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = ogcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		spt = seq[j];
		fpt = ogcp;
		gc = ( sgappat[j] == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( !gb *  gc ) *fpt += feff;
				fpt++;
			}
		}
	}
}

void new_OpeningGapCount_zure( float *ogcp, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )

{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = ogcp;
	i = len+2;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		spt = seq[j];
		fpt = ogcp;
		gc = ( sgappat[j] == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( !gb *  gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			if( !gb *  gc ) *fpt += feff;
		}
	}
}

void new_FinalGapCount_zure( float *fgcp, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )

{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len+2;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( sgappat[j] == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}

void new_FinalGapCount( float *fgcp, int clus, char **seq, double *eff, int len, char *egappat )
{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}

void st_OpeningGapCount( float *ogcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = ogcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		spt = seq[j];
		fpt = ogcp;
		gc = 0;
//		gc = 1;
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( !gb *  gc ) *fpt += feff;
				fpt++;
			}
		}
	}
	ogcp[len] = 0.0;
}

void st_FinalGapCount_zure( float *fgcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len+1;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		fpt = fgcp+1;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len;
//		for( i=1; i<len; i++ ) 
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = 0;
//			gc = 1;
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}

void st_FinalGapCount( float *fgcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len;
//		for( i=1; i<len; i++ ) 
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = 0;
//			gc = 1;
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}

void getGapPattern( float *fgcp, int clus, char **seq, double *eff, int len, char *xxx )
{
	int i, j, gc, gb; 
	float feff;
	float *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len+1;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (float)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len+1;
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
	}
	for( j=0; j<len; j++ )
	{
		reporterr(       "%d, %f\n", j, fgcp[j] );
	}
}

void getdigapfreq_st( float *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	float feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( 0 && seq[i][0] == '-' ) // machigai kamo
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] == '-' && seq[i][j-1] == '-' )
				freq[j] += feff;
		}
		if( 0 && seq[i][len-1] == '-' )
			freq[len] += feff;
	}
//	reporterr(       "\ndigapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdiaminofreq_x( float *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	float feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( seq[i][0] != '-' ) // tadashii
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] != '-' && seq[i][j-1] != '-' )
				freq[j] += feff;
		}
		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
			freq[len] += feff;
	}
//	reporterr(       "\ndiaaf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdiaminofreq_st( float *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	float feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( seq[i][0] != '-' )
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] != '-' && seq[i][j-1] != '-' )
				freq[j] += feff;
		}
//		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
			freq[len] += feff;
	}
//	reporterr(       "\ndiaaf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdigapfreq_part( float *freq, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )
{
	int i, j;
	float feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
//		if( seq[i][0] == '-' )
		if( seq[i][0] == '-' && sgappat[i] == '-' )
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] == '-' && seq[i][j-1] == '-' )
				freq[j] += feff;
		}
//		if( seq[i][len] == '-' && seq[i][len-1] == '-' ) // xxx wo tsukawanaitoki arienai
		if( egappat[i] == '-' && seq[i][len-1] == '-' )
			freq[len] += feff;
	}
//	reporterr(       "\ndigapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdiaminofreq_part( float *freq, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )
{
	int i, j;
	float feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( seq[i][0] != '-' && sgappat[i] != '-' )
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] != '-' && seq[i][j-1] != '-' )
				freq[j] += feff;
		}
//		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
		if( egappat[i] != '-'  && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
			freq[len] += feff;
	}
//	reporterr(       "\ndiaaf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getgapfreq_zure_part( float *freq, int clus, char **seq, double *eff, int len, char *sgap )
{
	int i, j;
	float feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( sgap[i] == '-' )
			freq[0] += feff;
		for( j=0; j<len; j++ ) 
		{
			if( seq[i][j] == '-' )
				freq[j+1] += feff;
		}
//		if( egap[i] == '-' )
//			freq[len+1] += feff;
	}
//	reporterr(       "\ngapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getgapfreq_zure( float *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	float feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		for( j=0; j<len; j++ ) 
		{
			if( seq[i][j] == '-' )
				freq[j+1] += feff;
		}
	}
	freq[len+1] = 0.0;
//	reporterr(       "\ngapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getgapfreq( float *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	float feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		for( j=0; j<len; j++ ) 
		{
			if( seq[i][j] == '-' )
				freq[j] += feff;
		}
	}
	freq[len] = 0.0;
//	reporterr(       "\ngapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void st_getGapPattern( Gappat **pat, int clus, char **seq, double *eff, int len )
{
	int i, j, k, gb, gc; 
	int known;
	float feff;
	Gappat **fpt;
	char *spt;
	int gaplen;

	fpt = pat;
	i = len+1;
	while( i-- ) 
	{
		if( *fpt ) free( *fpt );
		*fpt++ = NULL;
	}

	for( j=0; j<clus; j++ ) 
	{
//		reporterr(       "seq[%d] = %s\n", j, seq[j] );
		feff = (float)eff[j];

		fpt = pat;
		*fpt = NULL; // Falign.c kara yobareru tokiha chigau.
		spt = seq[j];
		gc = 0;
		gaplen = 0;

		for( i=0; i<len+1; i++ ) 
//		while( i-- )
		{
//			reporterr(       "i=%d, gaplen = %d\n", i, gaplen );
			gb = gc;
			gc = ( i != len && *spt++ == '-' );
			if( gc ) 
				gaplen++;
			else
			{
				if( gb && gaplen )
				{
					k = 1;
					known = 0;
					if( *fpt ) for( ; (*fpt)[k].len != -1; k++ )
					{
						if( (*fpt)[k].len == gaplen ) 
						{
//							reporterr(       "known\n" );
							known = 1;
							break;
						}
					}

					if( known == 0 )
					{
						*fpt = (Gappat *)realloc( *fpt, (k+3) *  sizeof( Gappat ) );  // mae1 (total), ato2 (len0), term
						if( !*fpt )
						{
							reporterr(       "Cannot allocate gappattern!'n" );
							reporterr(       "Use an approximate method, with the --mafft5 option.\n" );
							exit( 1 );
						}
						(*fpt)[k].freq = 0.0;
						(*fpt)[k].len = gaplen;
						(*fpt)[k+1].len = -1;
						(*fpt)[k+1].freq = 0.0; // iranai
//						reporterr(       "gaplen=%d, Unknown, %f\n", gaplen, (*fpt)[k].freq );
					}

//					reporterr(       "adding pos %d, len=%d, k=%d, freq=%f->", i, gaplen, k, (*fpt)[k].freq );
					(*fpt)[k].freq += feff;
//					reporterr(       "%f\n", (*fpt)[k].freq );
					gaplen = 0;
				}
			}
			fpt++;
		}
	}
#if 1
	for( j=0; j<len+1; j++ )
	{
		if( pat[j] )
		{
//			reporterr(       "j=%d\n", j );
//			for( i=1; pat[j][i].len!=-1; i++ )
//				reporterr(       "pos=%d, i=%d, len=%d, freq=%f\n", j, i, pat[j][i].len, pat[j][i].freq );

			pat[j][0].len = 0; // iminashi
			pat[j][0].freq = 0.0;
			for( i=1; pat[j][i].len!=-1;i++ )
			{
				pat[j][0].freq += pat[j][i].freq;
//				reporterr(       "totaling, i=%d, result = %f\n", i, pat[j][0].freq );
			}
//			reporterr(       "totaled, result = %f\n", pat[j][0].freq );

			pat[j][i].freq = 1.0 - pat[j][0].freq;
			pat[j][i].len = 0; // imiari
			pat[j][i+1].len = -1; 
		}
		else
		{
			pat[j] = (Gappat *)calloc( 3, sizeof( Gappat ) );
			pat[j][0].freq = 0.0;
			pat[j][0].len = 0; // iminashi

			pat[j][1].freq = 1.0 - pat[j][0].freq;
			pat[j][1].len = 0; // imiari
			pat[j][2].len = -1; 
		}
	}
#endif
}

static int minimum( int i1, int i2 )
{
	return MIN( i1, i2 );
}

static void commongappickpairfast( char *r1, char *r2, char *i1, char *i2, int *skip1, int *skip2 )
{
	CALLS && printf("called %s:commongappickpairfast()\n", __FILE__);
//	char *i1bk = i1;
	int skip, skipped1, skipped2;
//	int skip, skipped1, skipped2, scand1, scand2;
	skipped1 = skipped2 = 0;
//	reporterr("\n");
//	while( *i1 )
	while( 1 )
	{
//		fprintf( stderr, "i1 pos =%d\n", (int)(i1- i1bk) );
//		reporterr( "\nSkip cand %d-%d\n", *skip1-skipped1, *skip2-skipped2 );
		skip = minimum( *skip1-skipped1, *skip2-skipped2 );
//		reporterr( "Skip %d\n", skip );
		i1 += skip;
		i2 += skip;
		skipped1 += skip;
		skipped2 += skip;
//		fprintf( stderr, "i1 pos =%d, nlenmax=%d\n", (int)(i1- i1bk), nlenmax );
		if( !*i1 ) break;
//		reporterr( "%d, %c-%c\n", i1-i1bk, *i1, *i2 );
//		if( *i1 == '-' && *i2 == '-' )  // iranai?
//		{
//			reporterr( "Error in commongappickpairfast" );
//			exit( 1 );
//			i1++;
//			i2++;
//		}
		if( *i1 != '-' ) 
		{
			skipped1 = 0;
			skip1++;
		}
		else skipped1++;

		if( *i2 != '-' ) 
		{
			skipped2 = 0;
			skip2++;
		}
		else skipped2++;

		*r1++ = *i1++;
		*r2++ = *i2++;
	}
	*r1 = 0;
	*r2 = 0;
}

static void commongappickpair( char *r1, char *r2, char *i1, char *i2 )
{
	CALLS && printf("called %s:commongappickpair()\n", __FILE__);
//	strcpy( r1, i1 );
//	strcpy( r2, i2 );
//	return; // not SP
	while( *i1 )
	{
		if( *i1 == '-' && *i2 == '-' ) 
		{
			i1++;
			i2++;
		}
		else
		{
			*r1++ = *i1++;
			*r2++ = *i2++;
		}
	}
	*r1 = 0;
	*r2 = 0;
}

float naiveRpairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
//	return( 0 );
	int i, j;
	float val;
	float  valf;
	int  pv;
	double deff;
	char *p1, *p2, *p1p, *p2p;
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		deff = eff1[i] * eff2[j];
//		reporterr(       "feff %d-%d = %f\n", i, j, feff );
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		valf = 0;
		p1 = seq1[i]; p2 = seq2[j];
		pv = 0;
		if( *p1 == '-' && *p2 != '-' )
			pv = penal;
		if( *p1 != '-' && *p2 == '-' )
			pv = penal;
//		if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
		p1p = p1; p2p = p2;
		valf += (float)amino_dis[(int)*p1++][(int)*p2++] + 0.5 * pv;
		while( *p1p )
		{
			pv = 0;
			if( *p1p != '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p == '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
//					;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p != '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 == '-' )
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
					pv = penal;
//					;
			}
			if( *p1p == '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
					pv = penal;
//					;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
//			if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
			valf += amino_dis[(int)*p1++][(int)*p2++] + 0.5 * pv;
			p1p++; p2p++;
		}
//		reporterr(       "valf = %d\n", valf );
		val += deff * ( valf );
	}
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}
float naiveQpairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
	int i, j;
	float val;
	float  valf;
	int  pv;
	double deff;
	char *p1, *p2, *p1p, *p2p;
	return( 0 );
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		deff = eff1[i] * eff2[j];
//		reporterr(       "feff %d-%d = %f\n", i, j, feff );
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		valf = 0;
		p1 = seq1[i]; p2 = seq2[j];
		pv = 0;
		if( *p1 == '-' && *p2 != '-' )
			pv = penal;
		if( *p1 != '-' && *p2 == '-' )
			pv = penal;
//		if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
		p1p = p1; p2p = p2;
		valf += (float)amino_dis[(int)*p1++][(int)*p2++] + 0.5 * pv;
		while( *p1p )
		{
			pv = 0;
			if( *p1p != '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p == '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 == '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p != '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 == '-' )
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
			if( *p1p == '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
//			if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
			valf += amino_dis[(int)*p1++][(int)*p2++] + 0.5 * pv;
			p1p++; p2p++;
		}
//		reporterr(       "valf = %d\n", valf );
		val += deff * ( valf );
	}
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}
float naiveHpairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
	int i, j;
	float val;
	float  valf;
	int  pv;
//	float feff = 0.0; // by D.Mathog, a guess
	double deff;
	char *p1, *p2, *p1p, *p2p;
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		deff = eff1[i] * eff2[j];
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		valf = 0;
		p1 = seq1[i]; p2 = seq2[j];
		pv = 0;
		if( *p1 == '-' && *p2 != '-' )
			pv = penal;
		if( *p1 != '-' && *p2 == '-' )
			pv = penal;
		if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, (int)(p1-seq1[i]), (int)(p2-seq2[j]) );
		p1p = p1; p2p = p2;
		valf += (float)amino_dis[(int)*p1++][(int)*p2++] + 0.5 * pv;
		while( *p1p )
		{
			pv = 0;
			if( *p1p != '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p == '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 == '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p != '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 == '-' )
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
			if( *p1p == '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					;
				if( *p1 != '-' && *p2 == '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
//			if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
			valf += amino_dis[(int)*p1++][(int)*p2++] + 0.5 * pv;
			p1p++; p2p++;
		}
//		reporterr(       "valf = %d\n", valf );
		val += deff * ( valf );
	}
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}

float naivepairscorefast( char *seq1, char *seq2, int *skip1, int *skip2, int penal )
{
	float  vali;
	int len = strlen( seq1 );
	char *s1, *s2;
	char *p1, *p2;

	s1 = calloc( len+1, sizeof( char ) );
	s2 = calloc( len+1, sizeof( char ) );
	{
		vali = 0.0;
		commongappickpairfast( s1, s2, seq1, seq2, skip1, skip2 );
//		commongappickpair( s1, s2, seq1, seq2 );
//		reporterr(       "\n###s1 = %s\n", seq1 );
//		reporterr(       "###s2 = %s\n", seq2 );
//		reporterr(       "\n###i1 = %s\n", s1 );
//		reporterr(       "###i2 = %s\n", s2 );
//		reporterr( "allocated size, len+1 = %d\n", len+1 );
//		reporterr(       "###penal = %d\n", penal );

		p1 = s1; p2 = s2;
		while( *p1 )
		{
			if( *p1 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali += (float)penal;
//				while( *p1 == '-' || *p2 == '-' ) 
				while( *p1 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
			if( *p2 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali +=  (float)penal;
//				while( *p2 == '-' || *p1 == '-' ) 
				while( *p2 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
			vali += (float)amino_dis[(int)*p1++][(int)*p2++];
		}
	}
	free( s1 );
	free( s2 );
//	reporterr(       "###vali = %d\n", vali );
	return( vali );
}

float naivepairscore11( char *seq1, char *seq2, int penal )
{
	CALLS && printf("called %s:naivepairscore11()\n", __FILE__);
	float  vali;
	int len = strlen( seq1 );
	char *s1, *s2, *p1, *p2;

	s1 = calloc( len+1, sizeof( char ) );
	s2 = calloc( len+1, sizeof( char ) );
	{
		vali = 0.0;
		commongappickpair( s1, s2, seq1, seq2 );
//		reporterr(       "###i1 = %s\n", s1 );
//		reporterr(       "###i2 = %s\n", s2 );
//		reporterr(       "###penal = %d\n", penal );

		p1 = s1; p2 = s2;
		while( *p1 )
		{
			if( *p1 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali += (float)penal;
//				while( *p1 == '-' || *p2 == '-' ) 
				while( *p1 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
			if( *p2 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali +=  (float)penal;
//				while( *p2 == '-' || *p1 == '-' ) 
				while( *p2 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
			vali += (float)amino_dis[(int)*p1++][(int)*p2++];
		}
	}
	free( s1 );
	free( s2 );
//	reporterr(       "###vali = %d\n", vali );
	return( vali );
}

float naivepairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
//	return( 0.0 );
	int i, j;
	float val;
	int  vali;
	float feff;
	int len = strlen( seq1[0] );
	char *s1, *s2, *p1, *p2;
	s1 = calloc( len+1, sizeof( char ) );
	s2 = calloc( len+1, sizeof( char ) );
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		vali = 0;
		feff = eff1[i] * eff2[j];
//		reporterr(       "feff %d-%d = %f\n", i, j, feff );
		commongappickpair( s1, s2, seq1[i], seq2[j] );
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		p1 = s1; p2 = s2;
		while( *p1 )
		{
			if( *p1 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali += penal;
//				while( *p1 == '-' || *p2 == '-' ) 
				while( *p1 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
			if( *p2 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali +=  penal;
//				while( *p2 == '-' || *p1 == '-' ) 
				while( *p2 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
			vali += amino_dis[(int)*p1++][(int)*p2++];
		}
//		reporterr(       "vali = %d\n", vali );
		val += feff * vali;
	}
	free( s1 );
	free( s2 );
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}

double plainscore( int nseq, char **s )
{
	int i, j, ilim;
	double v = 0.0;
	
	ilim = nseq-1;
	for( i=0; i<ilim; i++ ) for( j=i+1; j<nseq; j++ )
	{
		v += (double)naivepairscore11( s[i], s[j], penalty );
	}

	reporterr(       "penalty = %d\n", penalty );

	return( v );
}

void intcat( int *s1, int *s2 )
{
	while( *s1 != -1 ) s1++;
	while( *s2 != -1 ) 
	{
//		reporterr(       "copying %d\n", *s2 );
		*s1++ = *s2++;
	}
	*s1 = -1;
}

void intcpy( int *s1, int *s2 )
{
	while( *s2 != -1 ) 
	{
//		reporterr(       "copying %d\n", *s2 );
		*s1++ = *s2++;
	}
	*s1 = -1;
}

void intncpy( int *s1, int *s2, int n )
{
	while( n-- ) *s1++ = *s2++;
}

void fltncpy( float *s1, float *s2, int n )
{
	while( n-- ) *s1++ = *s2++;
}

static int countmem( int *s )
{
	int v = 0;
	while( *s++ != -1 ) v++;
	return( v );
}

static int lastmem( int *s )
{
	while( *s++ != -1 ) 
		;
	return( *(s-2) );
}


int addonetip( int njobc, int ***topolc, float **lenc, float **iscorec, int ***topol, float **len, Treedep *dep, int treeout, Addtree *addtree, int iadd, char **name, int *alnleninnode, int *nogaplen, int noalign )
{
	int i, j, mem0, mem1, posinnew, m;
	int nstep;
	int norg;
	float minscore, minscoreo, eff0, eff1, addedlen, tmpmin;
	int nearest, nearesto;
	int repnorg;
	int *leaf2node;
	int *additionaltopol;
//	double (*clusterfuncpt[1])(double,double);
	Bchain *ac, *acpt, *acori, *acnext, *acprev;
	int neighbor;
	char *neighborlist;
	char *npt;
	int reflen, nearestnode, nogaplentoadd;
	int *topoldum0 = NULL;
	int *topoldum1 = NULL;
	int *topolo0;
	int *topolo1;
	int seqlengthcondition;
	double sueff1_double_local = 1.0 - sueff_global;
	double sueff05_double_local = sueff_global * 0.5;
//	char **tree; //static?
//	char *treetmp; //static?

//	for( i=0; i<njobc; i++ ) reporterr( "nogaplen of %d = %d\n", i+1, nogaplen[i] );
//exit( 1 );


//	treetmp = AllocateCharVec( njob*150 );
//	tree = AllocateCharMtx( njob, njob*150 );

//	sueff1_double = 1.0 - sueff_global;
//	sueff05_double = sueff_global * 0.5;
//	if ( treemethod == 'X' )
//		clusterfuncpt[0] = cluster_mix_double;
//	else if ( treemethod == 'E' )
//		clusterfuncpt[0] = cluster_average_double;
//	else if ( treemethod == 'q' )
//		clusterfuncpt[0] = cluster_minimum_double;
//	else
//	{
//		reporterr(       "Unknown treemethod, %c\n", treemethod );
//		exit( 1 );
//	}

	norg = njobc-1;
	nstep = njobc-2;

	additionaltopol = (int *)calloc( 2, sizeof( int ) );
	leaf2node= (int *)calloc( norg, sizeof( int ) );
	if( treeout )
	{
		neighborlist = calloc( norg * 30, sizeof( char ) );
	}
//	for( i=0; i<njobc; i++ ) sprintf( tree[i], "%d", i+1 );
	if( !leaf2node )
	{
		reporterr(       "Cannot allocate leaf2node.\n" );
		exit( 1 );
	}
	additionaltopol[0] = norg;
	additionaltopol[1] = -1;

	ac = (Bchain *)malloc( norg * sizeof( Bchain ) );
	for( i=0; i<norg; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[norg-1].next = NULL;


	acori = (Bchain *)malloc( 1 * sizeof( Bchain ) );
	acori->next = ac;
	acori->pos = -1;
	ac[0].prev = acori;


//	for( i=0; i<nstep; i++ )
//	{
//		reporterr(       "distfromtip = %f\n", dep[i].distfromtip );
//	}
//
//	for( i=0; i<norg; i++ )
//	{
//		reporterr(       "disttofrag(%d,%d) = %f\n", i, njobc-1, iscorec[i][norg-i] );
//	}


	minscore = 9999.9;
	nearest = -1;
	for( i=0; i<norg; i++ )
	{
		tmpmin = iscorec[i][norg-i];
		if( minscore > tmpmin )
		{
			minscore = tmpmin;
			nearest = i;
		}
	}
	nearesto = nearest;
	minscoreo = minscore;



//	for( i=0; i<njobc-1; i++ ) for( j=i+1; j<njobc; j++ )
//		reporterr(       "iscorec[%d][%d] = %f\n", i, j, iscorec[i][j-i] );
//	reporterr( "nearest = %d\n", nearest+1 );
//	reporterr( "nearesto = %d\n", nearesto+1 );

	posinnew = 0;
	repnorg = -1;
	nogaplentoadd = nogaplen[norg];



	for( i=0; i<norg; i++ ) leaf2node[i] = -1; 
	for( i=0; i<nstep; i++ )
	{
		mem0 = topol[i][0][0];
		mem1 = topol[i][1][0];
		nearestnode = leaf2node[nearest];
		if( nearestnode == -1 )
			reflen = nogaplen[nearest];
		else
			reflen = alnleninnode[nearestnode];
//			reflen = alnleninnode[i]; // BUG!!

		if( noalign ) seqlengthcondition =  1;
		else seqlengthcondition = ( nogaplentoadd <= reflen );

//seqlengthcondition = 1; // CHUUI
//seqlengthcondition = ( nogaplentoadd <= reflen ); // CHUUI

		if( repnorg == -1 && dep[i].distfromtip * 2 > minscore && seqlengthcondition )  // Keitouteki ichi ha fuseikaku.
//		if( repnorg == -1 && dep[i].distfromtip * 2 > minscore ) // Keitouteki ichi dake ga hitsuyouna baaiha kore wo tsukau.
		{
//			reporterr(       "INSERT HERE, %d-%d\n", nearest, norg );
//			reporterr(       "nearest = %d\n", nearest );
//			reporterr(       "\n\n\nminscore = %f\n", minscore );
//			reporterr(       "distfromtip *2 = %f\n", dep[i].distfromtip * 2 );
//			reporterr(       "nearest=%d, leaf2node[]=%d\n", nearest, leaf2node[nearest] );

			if( nearestnode == -1 )
			{
//				reporterr(       "INSERTING to 0!!!\n" );
//				reporterr(       "lastlength = %d\n", nogaplen[norg] );
//				reporterr(       "reflength = %d\n", nogaplen[nearest] );
				topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( 1 + 1 ) * sizeof( int ) );
				topolc[posinnew][0][0] = nearest;
				topolc[posinnew][0][1] = -1;

				addedlen = lenc[posinnew][0] = minscore / 2;

			}
			else
			{
//				reporterr(       "INSERTING to g, leaf2node = %d, cm=%d!!!\n", leaf2node[nearest], countmem(topol[leaf2node[nearest]][0] ) );
//				reporterr(       "alnleninnode[i] = %d\n", alnleninnode[i] );
//				reporterr(       "alnleninnode[leaf2node[nearest]] = %d\n", alnleninnode[leaf2node[nearest]] );

				topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( ( countmem( topol[nearestnode][0] ) + countmem( topol[nearestnode][1] ) + 1 ) * sizeof( int ) ) );
//				reporterr(       "leaf2node[%d] = %d\n", nearest, leaf2node[nearest] );
				intcpy( topolc[posinnew][0], topol[nearestnode][0] );
				intcat( topolc[posinnew][0], topol[nearestnode][1] );
//				addedlen = lenc[posinnew][0] = minscore / 2 - len[nearestnode][0]; // bug!!
				addedlen = lenc[posinnew][0] = dep[i].distfromtip - minscore / 2; // 2014/06/10
//				fprintf( stderr, "addedlen = %f, dep[i].distfromtip = %f, len[nearestnode][0] = %f, minscore/2 = %f, lenc[posinnew][0] = %f\n", addedlen, dep[i].distfromtip, len[nearestnode][0], minscore/2, lenc[posinnew][0] );

			}
			neighbor = lastmem( topolc[posinnew][0] );

			if( treeout )
			{
				addtree[iadd].nearest = nearesto;
				addtree[iadd].dist1 = minscoreo;
				addtree[iadd].dist2 = minscore;
				neighborlist[0] = 0;
				npt = neighborlist;
				for( j=0; topolc[posinnew][0][j]!=-1; j++ )
				{
					sprintf( npt, "%d ", topolc[posinnew][0][j]+1 );
					npt += strlen( npt );
				}
				addtree[iadd].neighbors = calloc( npt-neighborlist+1, sizeof( char ) );
				strcpy( addtree[iadd].neighbors, neighborlist );
			}

//			reporterr(       "INSERTING to 1!!!\n" );
			topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1], ( 1 + 1 ) * sizeof( int ) );
			topolc[posinnew][1][0] = norg;
			topolc[posinnew][1][1] = -1;
			lenc[posinnew][1] = minscore / 2;

//			reporterr(       "STEP %d (newnew)\n", posinnew );
//			for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j]+1 );
//			reporterr(       "\n len=%f\n", lenc[posinnew][0] );
//			for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j]+1 );
//			reporterr(       "\n len=%f\n", lenc[posinnew][1] );

			repnorg = nearest;

//			reporterr(       "STEP %d\n", posinnew );
//			for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j] );
//			reporterr(       "\n len=%f\n", lenc[i][0] );
//			for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j] );
//			reporterr(       "\n len=%f\n", lenc[i][1] );

//			im = topolc[posinnew][0][0];
//			jm = topolc[posinnew][1][0];
//			sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], lenc[posinnew][0], tree[jm], lenc[posinnew][1] );
//			strcpy( tree[im], treetmp );

			posinnew++;
		}

//		reporterr(       "minscore = %f\n", minscore );
//		reporterr(       "distfromtip = %f\n", dep[i].distfromtip );
//		reporterr(       "Modify matrix, %d-%d\n", nearest, norg );
		eff0 = iscorec[mem0][norg-mem0];
		eff1 = iscorec[mem1][norg-mem1];

//		iscorec[mem0][norg-mem0] = (clusterfuncpt[0])( eff0, eff1 );
		iscorec[mem0][norg-mem0] =  MIN( eff0, eff1 ) * sueff1_double_local + ( eff0 + eff1 ) * sueff05_double_local; 
		iscorec[mem1][norg-mem1] = 9999.9; // sukoshi muda

		acprev = ac[mem1].prev; 
		acnext = ac[mem1].next; 
		acprev->next = acnext;
		if( acnext != NULL ) acnext->prev = acprev;

		if( ( nearest == mem1 || nearest == mem0 ) )
		{
			minscore = 9999.9;
//			for( j=0; j<norg; j++ ) // sukoshi muda
//			{
//				if( minscore > iscorec[j][norg-j] )
//				{
//					minscore = iscorec[j][norg-j];
//					nearest = j;
//				}
//			}
//			reporterr(       "searching on modified ac " );
			for( acpt=acori->next; acpt!=NULL; acpt=acpt->next ) // sukoshi muda
			{
//				reporterr(       "." );
				j = acpt->pos;
				tmpmin = iscorec[j][norg-j];
				if( minscore > tmpmin )
				{
					minscore = tmpmin;
					nearest = j;
				}
			}
//			reporterr(       "done\n" );
		}

//		reporterr(       "posinnew = %d\n", posinnew );


		if( topol[i][0][0] == repnorg )
		{
			topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topol[i][0] ) + 2 ) * sizeof( int ) );
			intcpy( topolc[posinnew][0], topol[i][0] );
			intcat( topolc[posinnew][0], additionaltopol );
			lenc[posinnew][0] = len[i][0] - addedlen; // 2014/6/10
//			fprintf( stderr, "i=%d, dep[i].distfromtip=%f\n", i, dep[i].distfromtip );
//			fprintf( stderr, "addedlen=%f, len[i][0]=%f, lenc[][0]=%f\n", addedlen, len[i][0], lenc[posinnew][0] );
//			fprintf( stderr, "lenc[][1] = %f\n", lenc[posinnew][0] );
			addedlen = 0.0;
		}
		else
		{
			topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topol[i][0] ) + 1 ) * sizeof( int ) );
			intcpy( topolc[posinnew][0], topol[i][0] );
			lenc[posinnew][0] = len[i][0];
		}

		if( topol[i][1][0] == repnorg )
		{
			topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1], ( countmem( topol[i][1] ) + 2 ) * sizeof( int ) );
			intcpy( topolc[posinnew][1], topol[i][1] );
			intcat( topolc[posinnew][1], additionaltopol );
			lenc[posinnew][1] = len[i][1] - addedlen; // 2014/6/10
//			fprintf( stderr, "i=%d, dep[i].distfromtip=%f\n", i, dep[i].distfromtip );
//			fprintf( stderr, "addedlen=%f, len[i][1]=%f, lenc[][1]=%f\n", addedlen, len[i][1], lenc[posinnew][1] );
//			fprintf( stderr, "lenc[][1] = %f\n", lenc[posinnew][1] );
			addedlen = 0.0;

			repnorg = topolc[posinnew][0][0]; // juuyou
		}
		else
		{
			topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1], ( countmem( topol[i][1] ) + 1 ) * sizeof( int ) );
			intcpy( topolc[posinnew][1], topol[i][1] );
			lenc[posinnew][1] = len[i][1];
		}

//		reporterr(       "\nSTEP %d (new)\n", posinnew );
//		for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j]+1 );
//		reporterr(       "\n len=%f\n", lenc[posinnew][0] );
//		for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j]+1 );
//		reporterr(       "\n len=%f\n", lenc[posinnew][1] );

//		reporterr("\ni=%d\n####### leaf2node[nearest]= %d\n", i, leaf2node[nearest] );

		for( j=0; (m=topol[i][0][j])!=-1; j++ ) leaf2node[m] = i;
		for( j=0; (m=topol[i][1][j])!=-1; j++ ) leaf2node[m] = i;

//		reporterr("####### leaf2node[nearest]= %d\n", leaf2node[nearest] );

//		im = topolc[posinnew][0][0];
//		jm = topolc[posinnew][1][0];
//		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], lenc[posinnew][0], tree[jm], lenc[posinnew][1] );
//		strcpy( tree[im], treetmp );
//
//		reporterr(       "%s\n", treetmp );

		posinnew++;
	}

	if( nstep )
	{
		i--;
		topolo0 = topol[i][0];
		topolo1 = topol[i][1];
	}
	else
	{
//		i = 0;
//		free( topol[i][0] );//?
//		free( topol[i][1] );//?
//		topol[i][0] = calloc( 2, sizeof( int ) );
//		topol[i][1] = calloc( 1, sizeof( int ) );
//		topol[i][0][0] = 0;
//		topol[i][0][1] = -1;
//		topol[i][1][0] = -1;

		topoldum0 = calloc( 2, sizeof( int ) );
		topoldum1 = calloc( 1, sizeof( int ) );
		topoldum0[0] = 0;
		topoldum0[1] = -1;
		topoldum1[0] = -1;

		topolo0 = topoldum0;
		topolo1 = topoldum1;
	}
	if( repnorg == -1 )
	{
//		topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topol[i][0] ) + countmem( topol[i][1] ) + 1 ) * sizeof( int ) );
//		intcpy( topolc[posinnew][0], topol[i][0] );
//		intcat( topolc[posinnew][0], topol[i][1] );
		topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topolo0 ) + countmem( topolo1 ) + 1 ) * sizeof( int ) );
		intcpy( topolc[posinnew][0], topolo0 );
		intcat( topolc[posinnew][0], topolo1 );
//		lenc[posinnew][0] = len[i][0] + len[i][1] - minscore / 2; // BUG!! 2014/06/07 ni hakken
		if( nstep )
			lenc[posinnew][0] = minscore / 2 - dep[nstep-1].distfromtip; // only when nstep>0, 2014/11/21
		else
			lenc[posinnew][0] = minscore / 2;

//		reporterr( "\ndep[nstep-1].distfromtip = %f\n", dep[nstep-1].distfromtip );
//		reporterr( "lenc[][0] = %f\n", lenc[posinnew][0] );

		topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1],  2  * sizeof( int ) );
		intcpy( topolc[posinnew][1], additionaltopol );
		lenc[posinnew][1] = minscore / 2;

//		neighbor = lastmem( topolc[posinnew][0] );
		neighbor = norg-1;  // hakkirishita neighbor ga inai baai saigo ni hyouji

		if( treeout )
		{
			addtree[iadd].nearest = nearesto;
			addtree[iadd].dist1 = minscoreo;
			addtree[iadd].dist2 = minscore;
			neighborlist[0] = 0;
			npt = neighborlist;
			for( j=0; topolc[posinnew][0][j]!=-1; j++ )
			{
				sprintf( npt, "%d ", topolc[posinnew][0][j]+1 );
				npt += strlen( npt );
			}
			addtree[iadd].neighbors = calloc( npt-neighborlist+1, sizeof( char ) );
			strcpy( addtree[iadd].neighbors, neighborlist );
		}

//		reporterr(       "STEP %d\n", posinnew );
//		for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j] );
//		reporterr(       "\n len=%f", lenc[posinnew][0] );
//		for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j] );
//		reporterr(       "\n len=%f\n", lenc[posinnew][1] );
	}

	if( topoldum0 ) free( topoldum0 );
	if( topoldum1 ) free( topoldum1 );
	free( leaf2node );
	free( additionaltopol );
	free( ac );
	free( acori );
	if( treeout ) free( neighborlist );

	return( neighbor );
}

int samemember( int *mem, int *cand )
{
	CALLS && printf("called %s:samemember()\n", __FILE__);
	int i, j;
	int nm, nc;

	nm = 0; for( i=0; mem[i]>-1; i++ ) nm++;
	nc = 0; for( i=0; cand[i]>-1; i++ ) nc++;

	if( nm != nc ) return( 0 );

	for( i=0; mem[i]>-1; i++ )	
	{
		for( j=0; cand[j]>-1; j++ )
			if( mem[i] == cand[j] ) break;
		if( cand[j] == -1 ) return( 0 );
	}

	if( mem[i] == -1 )
	{
		return( 1 );
	}
	else
	{
		return( 0 );
	}
}

int samemembern( int *mem, int *cand, int nc )
{
	int i, j;
	int nm;

	nm = 0; 
	for( i=0; mem[i]>-1; i++ ) 
	{
		nm++;
		if( nm > nc ) return( 0 );
	}

	if( nm != nc ) return( 0 );

	for( i=0; mem[i]>-1; i++ )	
	{
		for( j=0; j<nc; j++ )
			if( mem[i] == cand[j] ) break;
		if( j == nc ) return( 0 );
	}

	if( mem[i] == -1 )
	{
		return( 1 );
	}
	else
	{
		return( 0 );
	}
}


int includemember( int *mem, int *cand ) // mem in cand 
{
	CALLS && printf("called %s:includemember()\n", __FILE__);
	int i, j;


	for( i=0; mem[i]>-1; i++ )
	{
		for( j=0; cand[j]>-1; j++ )
			if( mem[i] == cand[j] ) break;
		if( cand[j] == -1 ) return( 0 );
	}
//	reporterr(       "INCLUDED! mem[0]=%d\n", mem[0] );
	return( 1 );
}

int overlapmember( int *mem1, int *mem2 )
{
	int i, j;

	for( i=0; mem1[i]>-1; i++ )
		for( j=0; mem2[j]>-1; j++ )
			if( mem1[i] == mem2[j] ) return( 1 );
	return( 0 );
}
void gapcount( double *freq, char **seq, int nseq, double *eff, int lgth )
{
	int i, j;
	double fr;

//	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
//	return;

	for( i=0; i<lgth; i++ )
	{
		fr = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( seq[j][i] == '-' ) fr += eff[j];
		}
		freq[i] = fr;
//		reporterr(       "freq[%d] = %f\n", i, freq[i] );
	}
//	reporterr(       "\n" );
	return;
}

void gapcountf( float *freq, char **seq, int nseq, double *eff, int lgth )
{
	int i, j;
	double fr;

//	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
//	return;

	for( i=0; i<lgth; i++ )
	{
		fr = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( seq[j][i] == '-' ) fr += eff[j];
		}
		freq[i] = fr;
//		reporterr(       "freq[%d] = %f\n", i, freq[i] );
	}
//	reporterr(       "\n" );
	return;
}

void outgapcount( float *freq, int nseq, char *gappat, double *eff )
{
	int j;
	double fr;

	fr = 0.0;
	for( j=0; j<nseq; j++ )
	{
		if( gappat[j] == '-' ) fr += eff[j];
	}
	*freq = fr;
	return;
}

float dist2offset( float dist )
{
	float val = dist * 0.5 - specificityconsideration; // dist ha 0..2 dakara
//	float val = dist * 1.0 - specificityconsideration; // dist ha 0..2 dakara
	if( val > 0.0 ) val = 0.0;
	return val;
}

void makedynamicmtx( double **out, double **in, float offset )
{
	int i, j, ii, jj;
	double av;
 
	offset = dist2offset( offset * 2.0 ); // offset 0..1 -> 0..2

//	if( offset > 0.0 ) offset = 0.0;
//	reporterr(       "dynamic offset = %f\n", offset );

	for( i=0; i<nalphabets; i++ ) for( j=0; j<nalphabets; j++ )
	{
		out[i][j] = in[i][j];
	}
	if( offset == 0.0 ) return;

	for( i=0; i<nalphabets; i++ ) 
	{
		ii = (int)amino[i];
		if( ii == '-' ) continue; // text no toki arieru
		for( j=0; j<nalphabets; j++ )
		{
			jj = (int)amino[j];
			if( jj == '-' ) continue; // text no toki arieru
			out[i][j] = in[i][j] + offset * 600;
//			reporterr(       "%c-%c: %f\n", ii, jj, out[i][j] );
		}
	}

//	reporterr(       "offset = %f\n", offset );
//	reporterr(       "out[W][W] = %f\n", out[amino_n['W']][amino_n['W']] );
//	reporterr(       "out[A][A] = %f\n", out[amino_n['A']][amino_n['A']] );


	return;

// Taikaku youso no heikin ga 600 ni naruyouni re-scale.
// Hitaikaku youso ga ookiku narisugi.

	av = 0.0;
	for( i=0; i<nalphabets; i++ ) 
	{
		if( ii == '-' ) continue; // text no toki arieru
		av += out[i][i];
	}
	av /= (double)nalphabets;

	for( i=0; i<nalphabets; i++ ) 
	{
		if( amino[i] == '-' ) continue; // text no toki arieru
		for( j=0; j<nalphabets; j++ )
		{
			if( amino[j] == '-' ) continue; // text no toki arieru
			out[i][j] = out[i][j] * 600 / av;
			reporterr(       "%c-%c: %f\n", amino[i], amino[j], out[i][j] );
		}
	}
}
void FreeCommonIP()
{
	CALLS && printf("called %s:FreeCommonIP()\n", __FILE__);
	if( commonIP ) FreeIntMtx( commonIP ); 
	commonIP = NULL;
	commonAlloc1 = 0;
	commonAlloc2 = 0;
}

void makeskiptable( int n, int **skip, char **seq )
{
	CALLS && printf("called %s:makeskiptable()\n", __FILE__);
	char *nogapseq;
	int nogaplen, alnlen;
	int i, j, posinseq, gaplen;

	nogapseq = calloc( strlen( seq[0] )+1, sizeof( char ) );
	for( i=0; i<n; i++ )
	{
		gappick0( nogapseq, seq[i] );
		nogaplen = strlen( nogapseq );
		alnlen = strlen( seq[i] );
		skip[i] = calloc( nogaplen+1, sizeof( int ) );

//		reporterr( "%s\n", nogapseq );

		posinseq = 0;
		gaplen = 0;
		for( j=0; j<alnlen; j++ )
		{
			if( seq[i][j] == '-' )
			{
				skip[i][posinseq]++;
			}
			else
			{
				posinseq++;
			}
		}
//		for( j=0; j<nogaplen+1; j++ )
//			reporterr( "%d ", skip[i][j] );
//		reporterr( "\n" );
//		exit( 1 );
	}
	free( nogapseq );
}

int generatesubalignmentstable( int nseq, int ***tablept, int *nsubpt, int *maxmempt, int ***topol, double **len, double threshold )
{
	int i, j, rep0, rep1, nmem, mem;
	double distfromtip0, distfromtip1;
	double *distfromtip;
	reporterr( "\n\n\n" );

	*maxmempt = 0;
	*nsubpt = 0;

	distfromtip = calloc( nseq, sizeof( double ) );
	for( i=0; i<nseq-1; i++ )
	{

		rep0 = topol[i][0][0];
		distfromtip0 = distfromtip[rep0];
		distfromtip[rep0] += len[i][0];
//		reporterr( "distfromtip[%d] = %f->%f\n", rep0, distfromtip0, distfromtip[rep0] );

		rep1 = topol[i][1][0];
		distfromtip1 = distfromtip[rep1];
		distfromtip[rep1] += len[i][1];
//		reporterr( "distfromtip[%d] = %f->%f\n", rep1, distfromtip1, distfromtip[rep1] );

		if( topol[i][0][1] != -1 && distfromtip0 <= threshold && threshold < distfromtip[rep0] )
		{
//			reporterr( "HIT 0!\n" );
			*tablept = realloc( *tablept, sizeof( char * ) * (*nsubpt+2) );
			for( j=0, nmem=0; (mem=topol[i][0][j])!=-1; j++ ) 
				nmem++;
//			reporterr( "allocating %d\n", nmem+1 );
			(*tablept)[*nsubpt] = calloc( nmem+1, sizeof( int ) );
			(*tablept)[*nsubpt+1] = NULL;
			intcpy( (*tablept)[*nsubpt], topol[i][0] );
			if( *maxmempt < nmem ) *maxmempt = nmem;
			*nsubpt += 1;
		}

		if( topol[i][1][1] != -1 && distfromtip1 <= threshold && threshold < distfromtip[rep1] )
		{
//			reporterr( "HIT 1!\n" );
			*tablept = realloc( *tablept, sizeof( char * ) * (*nsubpt+2) );
			for( j=0, nmem=0; (mem=topol[i][1][j])!=-1; j++ )
				nmem++;
//			reporterr( "allocating %d\n", nmem+1 );
			(*tablept)[*nsubpt] = calloc( nmem+1, sizeof( int ) );
			(*tablept)[*nsubpt+1] = NULL;
			intcpy( (*tablept)[*nsubpt], topol[i][1] );
			if( *maxmempt < nmem ) *maxmempt = nmem;
			*nsubpt += 1;
		}

	}

	if( distfromtip[0] <= threshold ) 
	{
		free( distfromtip );
		return( 1 );
	}

	free( distfromtip );
	return( 0 );
}



float sumofpairsscore( int nseq, char **seq )
{
	CALLS && printf("called %s:sumofpairsscore()\n", __FILE__);
	float v = 0;
	int i, j;
	for( i=1; i<nseq; i++ )
	{
		for( j=0; j<i; j++ )
		{
			v += naivepairscore11( seq[i], seq[j], penalty ) / 600;
		}
	}
//	v /= ( (nseq-1) * nseq ) / 2;
	return( v );
}
