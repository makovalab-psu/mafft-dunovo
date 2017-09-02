#include "mltaln.h"

#ifdef PCALLS
#define CALLS 1
#else
#define CALLS 0
#endif
#ifdef PFILES
#define FILES 1
#else
#define FILES 0
#endif

static int upperCase = 0;

#define DEBUG   0
#define IODEBUG 0

char creverse( char f )
{
	CALLS && printf("called %s:creverse()\n", __FILE__);
	static char *table = NULL;
	if( table == NULL )
	{
		int i;
		table = AllocateCharVec(0x80);
		for( i=0; i<0x80; i++ ) table[i] = i;
		table['A'] = 'T';
		table['C'] = 'G';
		table['G'] = 'C';
		table['T'] = 'A';
		table['U'] = 'A';
		table['M'] = 'K';
		table['R'] = 'Y';
		table['W'] = 'W';
		table['S'] = 'S';
		table['Y'] = 'R';
		table['K'] = 'M';
		table['V'] = 'B';
		table['H'] = 'D';
		table['D'] = 'H';
		table['B'] = 'V';
		table['N'] = 'N';
		table['a'] = 't';
		table['c'] = 'g';
		table['g'] = 'c';
		table['t'] = 'a';
		table['u'] = 'a';
		table['m'] = 'k';
		table['r'] = 'y';
		table['w'] = 'w';
		table['s'] = 's';
		table['y'] = 'r';
		table['k'] = 'm';
		table['v'] = 'b';
		table['h'] = 'd';
		table['d'] = 'h';
		table['b'] = 'v';
		table['n'] = 'n';
//		table['-'] = '-';
//		table['.'] = '.';
//		table['*'] = '*';
	}
	return( table[(int)f] );
}

void sreverse( char *r, char *s )
{
	CALLS && printf("called %s:sreverse()\n", __FILE__);
	r += strlen( s );
	*r-- = 0;
	while( *s )
		*r-- = creverse( *s++ );
//		*r-- = ( *s++ );
}

void gappick_samestring( char *seq )
{
	CALLS && printf("called %s:gappick_samestring()\n", __FILE__);
	char *aseq = seq;

	for( ; *seq != 0; seq++ )
	{
		if( *seq != '-' )
			*aseq++ = *seq;
	}
	*aseq = 0;
}

#if 0

static int addlocalhom2( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, int skip )
{
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt;
	int st;
	int nlocalhom = 0;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	isumscore = 0;
	sumoverlap = 0;

#if 0
	fprintf( stderr, "nlocalhom = %d in addlocalhom\n", nlocalhom );
	fprintf( stderr, "al1 = %s, al2 = %s\n", al1, al2 );
	fprintf( stderr, "off1 = %d, off2 = %d\n", off1, off2 );
	fprintf( stderr, "localhopt = %p, skip = %d\n", localhompt, skip );
	fprintf( stderr, "pt1 = \n%s\n, pt2 = \n%s\n", pt1, pt2 );
#endif

	if( skip )
	{
		while( --skip > 0 ) localhompt = localhompt->next;
		localhompt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
		localhompt = localhompt->next;
//		fprintf( stderr, "tmppt = %p, localhompt = %p\n", tmppt, localhompt );
	}
	tmppt = localhompt;

	st = 0;
	iscore = 0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "In in while loop\n" );
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

#if 1
			isumscore += iscore;
			sumoverlap += end2-start2+1;
#else
			tmppt->overlapaa   = end2-start2+1;
			tmppt->opt = iscore * 5.8 / 600;
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "iscore (1)= %d\n", iscore );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			iscore = 0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			iscore += n_dis[(int)amino_n[(int)*pt1]][(int)amino_n[(int)*pt2]];
//			fprintf( stderr, "%c-%c, score(0) = %d\n", *pt1, *pt2, iscore );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}

	if( st )
	{
		if( nlocalhom++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}
		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;

#if 1
		isumscore += iscore;
		sumoverlap += end2-start2+1;
#else
		tmppt->overlapaa   = end2-start2+1;
		tmppt->opt = (double)iscore * 5.8 / 600;
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif
#if 0
		fprintf( stderr, "score (2)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
	{
		tmppt->overlapaa = sumoverlap;
		tmppt->opt = (double)sumscore * 5.8 / 600 / sumoverlap;
	}
	return( nlocalhom );
}

#endif



static int addlocalhom_r( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa, int skip )
{
	CALLS && printf("called %s:addlocalhom_r()\n", __FILE__);
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	double score;
	double sumscore;
	int sumoverlap;
	LocalHom *tmppt = NULL; // by D.Mathog, a guess
	int st;
	int nlocalhom = 0;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	sumscore = 0.0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

#if 0
	fprintf( stderr, "nlocalhom = %d in addlocalhom\n", nlocalhom );
	fprintf( stderr, "al1 = %s, al2 = %s\n", al1, al2 );
	fprintf( stderr, "off1 = %d, off2 = %d\n", off1, off2 );
	fprintf( stderr, "localhopt = %p, skip = %d\n", localhompt, skip );
#endif
	fprintf( stderr, "pt1 = \n%s\n, pt2 = \n%s\n", pt1, pt2 );

	if( skip )
	{
		while( --skip > 0 ) localhompt = localhompt->next;
		localhompt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
		localhompt = localhompt->next;
		fprintf( stderr, "tmppt = %p, localhompt = %p\n", (void *)tmppt, (void *)localhompt );
	}
	tmppt = localhompt;

	st = 0;
	score = 0.0;
	while( *pt1 != 0 )
	{
		fprintf( stderr, "In in while loop\n" );
		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

#if 1
			sumscore += score;
			sumoverlap += end2-start2+1;
#else
			tmppt->overlapaa   = end2-start2+1;
			tmppt->opt = score * 5.8 / 600;
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

			fprintf( stderr, "score (1)= %f\n", score );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
			score = 0.0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			score += (double)n_dis[(int)amino_n[(int)*pt1]][(int)amino_n[(int)*pt2]];
//			fprintf( stderr, "%c-%c, score(0) = %f\n", *pt1, *pt2, score );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( nlocalhom++ > 0 )
	{
//		fprintf( stderr, "reallocating ...\n" );
		tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//		fprintf( stderr, "done\n" );
		tmppt = tmppt->next;
		tmppt->next = NULL;
	}
	end1 = pos1 - 1;
	end2 = pos2 - 1;
	tmppt->start1 = start1;
	tmppt->start2 = start2;
	tmppt->end1   = end1  ;
	tmppt->end2   = end2  ;

#if 1
	sumscore += score;
	sumoverlap += end2-start2+1;
#else
	tmppt->overlapaa   = end2-start2+1;
	tmppt->opt = score * 5.8 / 600;
	tmppt->overlapaa   = overlapaa;
	tmppt->opt = (double)opt;
#endif

	fprintf( stderr, "score (2)= %f\n", score );
	fprintf( stderr, "al1: %d - %d\n", start1, end1 );
	fprintf( stderr, "al2: %d - %d\n", start2, end2 );

	for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
	{
		tmppt->overlapaa = sumoverlap;
		tmppt->opt = sumscore * 5.8 / 600 / sumoverlap;
	}
	return( nlocalhom );
}
void putlocalhom3( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa )
{
	CALLS && printf("called %s:putlocalhom3()\n", __FILE__);
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	double score;
	double sumscore;
	int sumoverlap;
	LocalHom *tmppt;
	LocalHom *subnosento;
	int st;
	int saisho;

	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	sumscore = 0.0;
	sumoverlap = 0;
	start1 = 0; // by Mathog, a guess
	start2 = 0; // by Mathog, a guess

	subnosento = localhompt;
	while( subnosento->next ) subnosento = subnosento->next;
	tmppt = subnosento;

	saisho = ( localhompt->nokori == 0 );

	fprintf( stderr, "localhompt = %p\n", (void *)localhompt );
	fprintf( stderr, "tmppt = %p\n", (void *)tmppt );
	fprintf( stderr, "subnosento = %p\n", (void *)subnosento );

	st = 0;
	score = 0.0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( localhompt->nokori++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
			}
			else
			{
				sumscore += score;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "score (1)= %f\n", score );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			score = 0.0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			score += (double)n_dis[(int)amino_n[(int)*pt1]][(int)amino_n[(int)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, score(0) = %f\n", *pt1, *pt2, score );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( *(pt1-1) != '-' && *(pt2-1) != '-'  )
	{
		if( localhompt->nokori++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}

		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;


#if 1
		if( divpairscore )
		{
			tmppt->overlapaa   = end2-start2+1;
			tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
		}
		else
		{
			sumscore += score;
			sumoverlap += end2-start2+1;
		}
#else
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif

#if 0
		fprintf( stderr, "score (2)= %f\n", score );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	fprintf( stderr, "sumscore = %f\n", sumscore );
	if( !divpairscore )
	{

		if( !saisho ) subnosento = subnosento->next;
		for( tmppt=subnosento; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
			tmppt->opt = sumscore * 5.8 / 600 / sumoverlap;
			fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
		}
	}
}
void putlocalhom_ext( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa )
{
	CALLS && printf("called %s:putlocalhom_ext()\n", __FILE__);
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;


	isumscore = 0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	st = 0;
	iscore = 0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			}
			else
			{
				isumscore += iscore;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "iscore (1)= %d\n", iscore );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			iscore = 0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			iscore += n_dis[(int)amino_n[(int)*pt1]][(int)amino_n[(int)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, iscore(0) = %d\n", *pt1, *pt2, iscore );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( *(pt1-1) != '-' && *(pt2-1) != '-'  )
	{
		if( nlocalhom++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}
		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;
	
#if 1
		if( divpairscore )
		{
			tmppt->overlapaa   = end2-start2+1;
			tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
		}
		else
		{
			isumscore += iscore;
			sumoverlap += end2-start2+1;
		}
#else
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif

#if 0
		fprintf( stderr, "iscore (2)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	if( !divpairscore )
	{
		for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
//			tmppt->opt = (double)isumscore * 5.8 / ( 600 * sumoverlap );
			tmppt->opt = (double)600 * 5.8 / 600;
//			fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
		}
	}
}

void putlocalhom_str( char *al1, char *al2, double *equiv, double scale, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa )
{
	CALLS && printf("called %s:putlocalhom_str()\n", __FILE__);
	int posinaln, pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
//	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;

	isumscore = 0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	posinaln = 0;
	while( *pt1 != 0 )
	{
		if( *pt1 != '-' && *pt2 != '-' && equiv[posinaln] > 0.0 )
		{
			start1 = end1 = pos1; start2 = end2 = pos2;
			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ... (posinaln=%d)\n", posinaln );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

			tmppt->overlapaa   = 1;
//			tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			tmppt->opt = equiv[posinaln] * scale;
//			fprintf( stdout, "*pt1=%c, *pt2=%c, equiv=%f\n", *pt1, *pt2, equiv[posinaln] );

		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
		posinaln++;
	}
}


void putlocalhom2( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa )
{
	CALLS && printf("called %s:putlocalhom2()\n", __FILE__);
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	int iscore;
	int isumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;


	isumscore = 0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	st = 0;
	iscore = 0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
			}
			else
			{
				isumscore += iscore;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "iscore (1)= %d\n", iscore );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			iscore = 0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			iscore += n_dis[(int)amino_n[(int)*pt1]][(int)amino_n[(int)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, iscore(0) = %d\n", *pt1, *pt2, iscore );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( *(pt1-1) != '-' && *(pt2-1) != '-'  )
	{
		if( nlocalhom++ > 0 )
		{
//			fprintf( stderr, "reallocating ...\n" );
			tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//			fprintf( stderr, "done\n" );
			tmppt = tmppt->next;
			tmppt->next = NULL;
		}
		end1 = pos1 - 1;
		end2 = pos2 - 1;
		tmppt->start1 = start1;
		tmppt->start2 = start2;
		tmppt->end1   = end1  ;
		tmppt->end2   = end2  ;
	
#if 1
		if( divpairscore )
		{
			tmppt->overlapaa   = end2-start2+1;
			tmppt->opt = (double)iscore / tmppt->overlapaa * 5.8 / 600;
		}
		else
		{
			isumscore += iscore;
			sumoverlap += end2-start2+1;
		}
#else
		tmppt->overlapaa   = overlapaa;
		tmppt->opt = (double)opt;
#endif

#if 0
		fprintf( stderr, "iscore (2)= %d\n", iscore );
		fprintf( stderr, "al1: %d - %d\n", start1, end1 );
		fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
	}

	if( !divpairscore )
	{
		for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
			tmppt->opt = (double)isumscore * 5.8 / ( 600 * sumoverlap );
//			fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
		}
	}
}
void putlocalhom( char *al1, char *al2, LocalHom *localhompt, int off1, int off2, int opt, int overlapaa )
{
	CALLS && printf("called %s:putlocalhom()\n", __FILE__);
	int pos1, pos2, start1, start2, end1, end2;
	char *pt1, *pt2;
	double score;
	double sumscore;
	int sumoverlap;
	LocalHom *tmppt = localhompt;
	int nlocalhom = 0;
	int st;
	pt1 = al1; pt2 = al2;
	pos1 = off1; pos2 = off2;


	sumscore = 0.0;
	sumoverlap = 0;
	start1 = 0; // by D.Mathog, a guess
	start2 = 0; // by D.Mathog, a guess

	st = 0;
	score = 0.0;
	while( *pt1 != 0 )
	{
//		fprintf( stderr, "pt = %c, %c, st=%d\n", *pt1, *pt2, st );
		if( st == 1 && ( *pt1 == '-' || *pt2 == '-' ) )
		{
			end1 = pos1 - 1;
			end2 = pos2 - 1;

			if( nlocalhom++ > 0 )
			{
//				fprintf( stderr, "reallocating ...\n" );
				tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//				fprintf( stderr, "done\n" );
				tmppt = tmppt->next;
				tmppt->next = NULL;
			}
			tmppt->start1 = start1;
			tmppt->start2 = start2;
			tmppt->end1   = end1  ;
			tmppt->end2   = end2  ;

#if 1
			if( divpairscore )
			{
				tmppt->overlapaa   = end2-start2+1;
				tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
			}
			else
			{
				sumscore += score;
				sumoverlap += end2-start2+1;
			}
#else
			tmppt->overlapaa   = overlapaa;
			tmppt->opt = (double)opt;
#endif

#if 0
			fprintf( stderr, "score (1)= %f\n", score );
			fprintf( stderr, "al1: %d - %d\n", start1, end1 );
			fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif
			score = 0.0;
			st = 0;
		}
		else if( *pt1 != '-' && *pt2 != '-' )
		{
			if( st == 0 )
			{
				start1 = pos1; start2 = pos2;
				st = 1;
			}
			score += (double)n_dis[(int)amino_n[(int)*pt1]][(int)amino_n[(int)*pt2]]; // - offset はいらないかも
//			fprintf( stderr, "%c-%c, score(0) = %f\n", *pt1, *pt2, score );
		}
		if( *pt1++ != '-' ) pos1++;
		if( *pt2++ != '-' ) pos2++;
	}
	if( nlocalhom++ > 0 )
	{
//		fprintf( stderr, "reallocating ...\n" );
		tmppt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
//		fprintf( stderr, "done\n" );
		tmppt = tmppt->next;
		tmppt->next = NULL;
	}
	end1 = pos1 - 1;
	end2 = pos2 - 1;
	tmppt->start1 = start1;
	tmppt->start2 = start2;
	tmppt->end1   = end1  ;
	tmppt->end2   = end2  ;

#if 1
	if( divpairscore )
	{
		tmppt->overlapaa   = end2-start2+1;
		tmppt->opt = score / tmppt->overlapaa * 5.8 / 600;
	}
	else
	{
		sumscore += score;
		sumoverlap += end2-start2+1;
	}
#else
	tmppt->overlapaa   = overlapaa;
	tmppt->opt = (double)opt;
#endif

#if 0
	fprintf( stderr, "score (2)= %f\n", score );
	fprintf( stderr, "al1: %d - %d\n", start1, end1 );
	fprintf( stderr, "al2: %d - %d\n", start2, end2 );
#endif

	if( !divpairscore )
	{
		for( tmppt=localhompt; tmppt; tmppt=tmppt->next )
		{
			tmppt->overlapaa = sumoverlap;
			tmppt->opt = sumscore * 5.8 / 600 / sumoverlap;
//			fprintf( stderr, "tmpptr->opt = %f\n", tmppt->opt );
		}
	}
}

char *cutal( char *al, int al_display_start, int start, int end )
{
	CALLS && printf("called %s:cutal()\n", __FILE__);
	int pos;
	char *pt = al;
	char *val = NULL;

	pos = al_display_start;
	do
	{
		if( start == pos ) val = pt;
		if( end == pos ) break;
//		fprintf( stderr, "pos=%d, *pt=%c, val=%p\n", pos, *pt, val );
		if( *pt != '-' ) pos++;
	} while( *pt++ != 0 );
	*(pt+1) = 0;
	return( val );
}

void ErrorExit( char *message )
{
	CALLS && printf("called %s:ErrorExit()\n", __FILE__);
	FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
	fprintf( stderr, "%s\n", message );
	exit( 1 );
}

void strncpy_caseC( char *str1, char *str2, int len )
{
	CALLS && printf("called %s:strncpy_caseC()\n", __FILE__);
	if( dorp == 'd' && upperCase > 0 ) 
	{
		while( len-- )
			*str1++ = toupper( *str2++ );
	}
	else strncpy( str1, str2, len );
}
	
void seqUpper( int nseq, char **seq )
{
	CALLS && printf("called %s:seqUpper()\n", __FILE__);
	int i, j, len;
	for( i=0; i<nseq; i++ ) 
	{
		len = strlen( seq[i] );
		for( j=0; j<len; j++ ) 
			seq[i][j] = toupper( seq[i][j] );
	}
}

void seqLower( int nseq, char **seq ) // yes
{
	CALLS && printf("called %s:seqLower()\n", __FILE__);
	int i, j, len;
	for( i=0; i<nseq; i++ ) 
	{
		len = strlen( seq[i] );
		for( j=0; j<len; j++ ) 
			seq[i][j] = tolower( seq[i][j] );
	}
}

int getaline_fp_eof( char *s, int l, FILE *fp )  /* end of file -> return 1 */
{
	CALLS && printf("called %s:getaline_fp_eof()\n", __FILE__);
    int c, i = 0 ;
    int noteofflag = 0;
    for( i=0; i<l && ( noteofflag = ( (c=getc(fp)) != EOF ) ) && c != '\n'; i++ ) 
    	*s++ = c;
    *s = '\0' ;
     return( !noteofflag );
}

int getaline_fp_eof_new(s, l, fp)  /* end of file -> return 1 */
char    s[] ; int l ; FILE *fp ;
{
		CALLS && printf("called %s:getaline_fp_eof_new()\n", __FILE__);
        int c = 0, i = 0 ;
		int noteofflag = 0;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( noteofflag = ( (c=getc(fp)) != EOF ) ) && c != '\n'; i++ ) 
        { *s++ = c ; }
        *s = '\0' ;
		if( c != '\n' && c != EOF ) while( getc(fp) != '\n' )
			;
		return( !noteofflag );
}

// fgets() reads at most l characters into buffer s from file fp, stopping at newlines and EOFs.
int myfgets(s, l, fp) // yes
char    s[] ; int l ; FILE *fp ;
{
		CALLS && printf("called %s:myfgets()\n", __FILE__);
        int     c = 0, i = 0 ;

		if( feof( fp ) ) return( 1 );

		for( i=0; i<l && ( c=getc( fp ) ) != '\n'; i++ ) 
        	*s++ = c;
        *s = '\0' ;
		if( c != '\n' ) 
			while( getc(fp) != '\n' )
				;
		return( 0 );
}

float input_new( FILE *fp, int d )
{
	CALLS && printf("called %s:input_new()\n", __FILE__);
	char mojiretsu[10];
	int i, c;

	c = getc( fp );
	if( c != '\n' )
		ungetc( c, fp );

	for( i=0; i<d; i++ )
		mojiretsu[i] = getc( fp );
	mojiretsu[i] = 0;

	return( atof( mojiretsu ) );
}


void PreRead( FILE *fp, int *locnjob, int *locnlenmax )
{
	CALLS && printf("called %s:PreRead()\n", __FILE__);
	int i, nleni;
	char b[B];

	fgets( b, B-1, fp ); *locnjob = atoi( b );
	*locnlenmax = 0;
	i=0; 
	while( i<*locnjob )
	{
		fgets( b, B-1, fp );
		if( !strncmp( b, "=", 1 ) )
		{
			i++;
			fgets( b, B-1, fp ); nleni = atoi( b );
			if( nleni > *locnlenmax ) *locnlenmax = nleni;
		}
	}
	if( *locnlenmax > N )
	{
		fprintf( stderr, "TOO LONG SEQUENCE!\n" );
		exit( 1 );
	}
	if( njob > M  ) 
	{
		fprintf( stderr, "TOO MANY SEQUENCE!\n" );
		fprintf( stderr, "%d > %d\n", njob, M );
		exit( 1 );
	}
}	

int allSpace( char *str )
{
	CALLS && printf("called %s:allSpace()\n", __FILE__);
	int value = 1;
	while( *str ) value *= ( !isdigit( *str++ ) );
	return( value );
}
	
void Read( char name[M][B], int nlen[M], char **seq )
{
	CALLS && printf("called %s:Read()\n", __FILE__);
	extern void FRead( FILE *x, char y[M][B], int z[M], char **w );
	FRead( stdin, name, nlen, seq );
}


void FRead( FILE *fp, char name[][B], int nlen[], char **seq )
{
	CALLS && printf("called %s:FRead()\n", __FILE__);
	int i, j; 
	char b[B];

	fgets( b, B-1, fp );
#if DEBUG
	fprintf( stderr, "b = %s\n", b );
#endif

    if( strstr( b, "onnet" ) ) scoremtx = 1;
    else if( strstr( b, "DnA" ) ) 
	{
		scoremtx = -1; 
		upperCase = -1;
	}
    else if( strstr( b, "dna" ) ) 
	{
		scoremtx = -1; 
		upperCase = 0;
	}
	else if( strstr( b, "DNA" ) )
	{
		scoremtx = -1; 
		upperCase = 1;
	}
    else if( strstr( b, "M-Y" ) || strstr( b, "iyata" ) ) scoremtx = 2; 
    else scoremtx = 0;
#if DEBUG
	fprintf( stderr, " %s->scoremtx = %d\n", b, scoremtx );
#endif

	geta2 = GETA2;

#if 0
	if( strlen( b ) >=25 )
	{
		b[25] = 0;
	#if DEBUG
		fprintf( stderr, "kimuraR = %s\n", b+20 );
	#endif
		kimuraR = atoi( b+20 );

		if( kimuraR < 0 || 20 < kimuraR ) ErrorExit( "Illeagal kimuraR value.\n" );
		if( allSpace( b+20 ) ) kimuraR = NOTSPECIFIED;
	}
	else kimuraR = NOTSPECIFIED;
	#if DEBUG
	fprintf( stderr, "kimuraR = %d\n", kimuraR );
	#endif

	if( strlen( b ) >=20 )
	{
		b[20] = 0;
	#if DEBUG
		fprintf( stderr, "pamN = %s\n", b+15 );
	#endif
		pamN = atoi( b+15 );
		if( pamN < 0 || 400 < pamN ) ErrorExit( "Illeagal pam value.\n" );
		if( allSpace( b+15 ) ) pamN = NOTSPECIFIED;
	}
	else pamN = NOTSPECIFIED;

	if( strlen( b ) >= 15 )
	{
		b[15] = 0;
	#if DEBUG
		fprintf( stderr, "poffset = %s\n", b+10 );
	#endif
		poffset = atoi( b+10 );
		if( poffset > 500 ) ErrorExit( "Illegal extending gap ppenalty\n" );
		if( allSpace( b+10 ) ) poffset = NOTSPECIFIED;
	}
	else poffset = NOTSPECIFIED;

	if( strlen( b ) >= 10 )
	{
		b[10] = 0;
	#if DEBUG
		fprintf( stderr, "ppenalty = %s\n", b+5 );
	#endif
		ppenalty = atoi( b+5 );
		if( ppenalty > 0 ) ErrorExit( "Illegal opening gap ppenalty\n" );
		if( allSpace( b+5 ) ) ppenalty = NOTSPECIFIED;
	}
	else ppenalty = NOTSPECIFIED;
#endif

	for( i=0; i<njob; i++ )
	{
		getaline_fp_eof_new( b, B-1, fp );
		strcpy( name[i], b );
#if DEBUG
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		fgets( b, B-1, fp ); nlen[i] = atoi( b );      /* seq i no nagasa */
		seq[i][0] = 0;
		if( nlen[i] ) for( j=0; j <= (nlen[i]-1)/C; j++ )
		{
			getaline_fp_eof_new( b, B-1, fp );
			/*	b[C] = 0;  */
			strcat( seq[i], b );
		} 
		seq[i][nlen[i]] = 0;
	}
	if( scoremtx == -1 && upperCase != -1 ) seqLower( njob, seq );
}


static int countKUorWA( FILE *fp ) // yes
{
	CALLS && printf("called %s:countKUorWA()\n", __FILE__);
	int value;
	int c, b;

	value= 0;
	b = '\n';
	FILES && printf("file read %d %s:%d\n", fp, __FILE__, __LINE__);
	while( ( c = getc( fp ) ) != EOF )
	{
		if( b == '\n' && ( c == '>' ) )
			value++;
		b = c;
	}
	rewind( fp );
	return( value );
}

void searchKUorWA( FILE *fp ) // yes
{
	CALLS && printf("called %s:searchKUorWA()\n", __FILE__);
	int c, b;
	b = '\n';
	FILES && printf("file read %d %s:%d\n", fp, __FILE__, __LINE__);
	while( !( ( ( c = getc( fp ) ) == '>' || c == EOF ) && b == '\n' ) )
		b = c;
	ungetc( c, fp );
}

static int onlyGraph( char *str )
{
	CALLS && printf("called %s:onlyGraph()\n", __FILE__);
	char tmp;
	char *res = str;
	char *bk = str;

//	while( (tmp=*str++) ) if( isgraph( tmp ) ) *res++ = tmp;
	while( (tmp=*str++) ) 
	{
		if( 0x20 < tmp && tmp < 0x7f ) *res++ = tmp;
		if( tmp == '>' || tmp == '(' )
		{
			fprintf( stderr, "========================================================\n" );
			fprintf( stderr, "========================================================\n" );
			fprintf( stderr, "=== \n" );
			fprintf( stderr, "=== ERROR!! \n" );
			fprintf( stderr, "=== In the '--anysymbol' and '--preservecase' modes, \n" );
			fprintf( stderr, "=== '>' and '(' are unacceptable.\n" );
			fprintf( stderr, "=== \n" );
			fprintf( stderr, "========================================================\n" );
			fprintf( stderr, "========================================================\n" );
			exit( 1 );
		}
	}
	*res = 0;
	return( res - bk );
}

static int charfilter( char *str )
{
	CALLS && printf("called %s:charfilter()\n", __FILE__);
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
	{
		if( tmp == '=' || tmp == '*' || tmp == '<' || tmp == '>' || tmp == '(' || tmp == ')' )
		{
			fprintf( stderr, "\n" );
			fprintf( stderr, "Characters '= * < > ( )' are not accepted in the --text mode, \nalthough most printable characters are ok.\n" );
			fprintf( stderr, "\n" );
			exit( 1 );
		}
		if( 0x20 < tmp && tmp < 0x7f )
//		if( tmp != '\n' && tmp != ' ' && tmp != '\t' ) // unprintable characters mo ok.
			*res++ = tmp;
	}
	*res = 0;
	return( res - bk );
}
static int onlyAlpha_lower( char *str ) // yes
{
	CALLS && printf("called %s:onlyAlpha_lower()\n", __FILE__);
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
		if( isalpha( tmp ) || tmp == '-' || tmp == '*' || tmp == '.' )
			*res++ = tolower( tmp );
	*res = 0;
	return( res - bk );
}
static int onlyAlpha_upper( char *str )
{
	CALLS && printf("called %s:onlyAlpha_upper()\n", __FILE__);
	char tmp;
	char *res = str;
	char *bk = str;

	while( (tmp=*str++) )
		if( isalpha( tmp ) || tmp == '-' || tmp == '*' || tmp == '.' )
			*res++ = toupper( tmp );
	*res = 0;
	return( res - bk );
}

void kake2hiku( char *str ) // yes
{
	CALLS && printf("called %s:kake2hiku()\n", __FILE__);
	do
		if( *str == '*' ) *str = '-';
	while( *str++ );
}

char *load1SeqWithoutName_realloc_casepreserve( FILE *fpp )
{
	CALLS && printf("called %s:load1SeqWithoutName_realloc_casepreserve()\n", __FILE__);
	int c, b;
	char *cbuf;
	int size = N;
	char *val;

	val = malloc( (size+1) * sizeof( char ) );
	cbuf = val;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&           
          !( ( c == '>' || c == '(' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;  /* 長すぎてもしらない */
		if( cbuf - val == size )
		{
			size += N;
			fprintf( stderr, "reallocating...\n" );
			val = (char *)realloc( val, (size+1) * sizeof( char ) );
			if( !val )
			{
				fprintf( stderr, "Allocation error in load1SeqWithoutName_realloc \n" );
				exit( 1 );
			}
			fprintf( stderr, "done.\n" );
			cbuf = val + size-N;
		}
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;
	onlyGraph( val );
//	kake2hiku( val );
	return( val );
}

char *load1SeqWithoutName_realloc( FILE *fpp ) // yes
{
	CALLS && printf("called %s:load1SeqWithoutName_realloc()\n", __FILE__);
	int c, b;
	char *cbuf;
	int size = N;
	char *val;

	val = malloc( (size+1) * sizeof( char ) );
	cbuf = val;

	b = '\n';
	FILES && printf("file read %d %s:%d\n", fpp, __FILE__, __LINE__);
	while( ( c = getc( fpp ) ) != EOF &&           
          !( ( c == '>' || c == '(' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;  /* [mojibake] */
		if( cbuf - val == size )
		{
			size += N;
			FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
			fprintf( stderr, "reallocating...\n" );
			val = (char *)realloc( val, (size+1) * sizeof( char ) );
			if( !val )
			{
				fprintf( stderr, "Allocation error in load1SeqWithoutName_realloc \n" );
				exit( 1 );
			}
			fprintf( stderr, "done.\n" );
			cbuf = val + size-N;
		}
		b = c;
	}
	FILES && printf("file seek %d %s:%d\n", fpp, __FILE__, __LINE__);
	ungetc( c, fpp );
	*cbuf = 0;

	if( nblosum == -2 )
	{
		charfilter( val );
	}
	else
	{
		if( dorp == 'd' )
			onlyAlpha_lower( val );
		else
			onlyAlpha_upper( val );
		kake2hiku( val );
	}
	return( val );
}

int load1SeqWithoutName_new( FILE *fpp, char *cbuf )
{
	CALLS && printf("called %s:load1SeqWithoutName_new()\n", __FILE__);
	int c, b;
	char *bk = cbuf;

	b = '\n';
	while( ( c = getc( fpp ) ) != EOF &&                    /* by T. Nishiyama */
          !( ( c == '>' || c == '(' || c == EOF ) && b == '\n' ) )
	{
		*cbuf++ = (char)c;  /* 長すぎてもしらない */
		b = c;
	}
	ungetc( c, fpp );
	*cbuf = 0;
	if( dorp == 'd' )
		onlyAlpha_lower( bk );
	else
		onlyAlpha_upper( bk );
	kake2hiku( bk );
	return( 0 );
}


void readDataforgaln( FILE *fp, char **name, int *nlen, char **seq )
{
	CALLS && printf("called %s:readDataforgaln()\n", __FILE__);
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		nlen[i] = strlen( seq[i] );
		free( tmpseq );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}

void readData_varlen( FILE *fp, char **name, int *nlen, char **seq )
{
	CALLS && printf("called %s:readData_varlen()\n", __FILE__);
	int i; 
	static char *tmpseq = NULL;

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		nlen[i] = strlen( tmpseq );
//		fprintf( stderr, "nlen[%d] = %d\n", i+1, nlen[i] );
		seq[i] = calloc( nlen[i]+1, sizeof( char ) );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}

void readData_pointer2( FILE *fp, int nseq, char **name, int *nlen, char **seq )
{
	CALLS && printf("called %s:readData_pointer2()\n", __FILE__);
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<nseq; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
		nlen[i] = strlen( seq[i] );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( nseq, seq );
#if 0
	free( tmpseq );
#endif
	if( outnumber )
	{
		char *namebuf;
		char *cptr;
		namebuf = calloc( B+100, sizeof( char ) );
		for( i=0; i<nseq; i++ )
		{
			namebuf[0] = '=';
			cptr = strstr( name[i], "_numo_e_" );
			if( cptr )
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, cptr+8 );
			else
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, name[i]+1 );
			strncpy( name[i], namebuf, B );
			name[i][B-1] = 0;
		}
		free( namebuf );
//		exit( 1 );
	}
}


void readData_pointer_casepreserve( FILE *fp, char **name, int *nlen, char **seq )
{
	CALLS && printf("called %s:readData_pointer_casepreserve()\n", __FILE__);
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
		nlen[i] = strlen( seq[i] );
	}
}


int copydatafromgui( char **namegui, char **seqgui, char **name, int *nlen, char **seq )
{
	CALLS && printf("called %s:copydatafromgui()\n", __FILE__);
	int i; 


	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; 
		strncpy( name[i]+1, namegui[i], B-2 );
		name[i][B-1] = 0;

		strcpy( seq[i], seqgui[i] );
		nlen[i] = strlen( seq[i] );
	}
	if( dorp == 'd' ) 
		seqLower( njob, seq );
	else if( dorp == 'p' ) 
		seqUpper( njob, seq );
	else
	{
		reporterr( "DNA or Protein?\n" );
		return( 1 );
	}
#if 0
	free( tmpseq );
#endif
	if( outnumber )
	{
		char *namebuf;
		char *cptr;
		namebuf = calloc( B+100, sizeof( char ) );
		for( i=0; i<njob; i++ )
		{
			namebuf[0] = '=';
			cptr = strstr( name[i], "_numo_e_" );
			if( cptr )
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, cptr+8 );
			else
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, name[i]+1 );
			strncpy( name[i], namebuf, B );
			name[i][B-1] = 0;
		}
		free( namebuf );
	}
	return( 0 );
}

void readData_pointer( FILE *fp, char **name, int *nlen, char **seq ) // yes
{
	CALLS && printf("called %s:readData_pointer()\n", __FILE__);
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	FILES && printf("file seek %d %s:%d\n", fp, __FILE__, __LINE__);
	rewind( fp );
	searchKUorWA( fp );

	FILES && printf("file read %d %s:%d\n", fp, __FILE__, __LINE__);
	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		FILES && printf("file read %d %s:%d\n", fp, __FILE__, __LINE__);
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		free( tmpseq );
		nlen[i] = strlen( seq[i] );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
	if( outnumber )
	{
		char *namebuf;
		char *cptr;
		namebuf = calloc( B+100, sizeof( char ) );
		for( i=0; i<njob; i++ )
		{
			namebuf[0] = '=';
			cptr = strstr( name[i], "_numo_e_" );
			if( cptr )
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, cptr+8 );
			else
				sprintf( namebuf+1, "_numo_s_%08d_numo_e_%s", i+1, name[i]+1 );
			strncpy( name[i], namebuf, B );
			name[i][B-1] = 0;
		}
		free( namebuf );
//		exit( 1 );
	}
}

void readData( FILE *fp, char name[][B], int nlen[], char **seq )
{
	CALLS && printf("called %s:readData()\n", __FILE__);
	int i; 
	static char *tmpseq = NULL;

#if 0
	if( !tmpseq )
	{
		tmpseq = AllocateCharVec( N );
	}
#endif

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		name[i][0] = '='; getc( fp ); 
#if 0
		fgets( name[i]+1, B-2, fp ); 
		j = strlen( name[i] );
		if( name[i][j-1] != '\n' )
			ErrorExit( "Too long name\n" );
		name[i][j-1] = 0;
#else
		myfgets( name[i]+1, B-2, fp ); 
#endif
#if 0
		fprintf( stderr, "name[%d] = %s\n", i, name[i] );
#endif
		tmpseq = load1SeqWithoutName_realloc( fp );
		strcpy( seq[i], tmpseq );
		nlen[i] = strlen( seq[i] );
		free( tmpseq );
	}
	if( dorp == 'd' && upperCase != -1 ) seqLower( njob, seq );
#if 0
	free( tmpseq );
#endif
}

void cutAlignment( FILE *fp, int **regtable, char **revtable, int *outtable, char **name, char **outseq )
{
	CALLS && printf("called %s:cutAlignment()\n", __FILE__);
	int i, j; 
	int outlen;
	static char *tmpseq = NULL;
	static char *dumname = NULL;
	char *fs, *rs;
	int npos, lpos;
	int startpos, endpos, seqlen;

	if( dumname == NULL )
	{
		dumname = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );


	npos = 0;
	for( i=0; i<njob; i++ )
	{
		dumname[0] = '>'; getc( fp ); 
		myfgets( dumname+1, B-1, fp ); 
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );

		if( outtable[i] )
		{
//			putc( '>', stdout );
//			puts( dumname+1 );

	
			strncat( name[npos], dumname, B-1 );
			name[npos][B-1] = 0;
	
			if( dorp == 'd' && upperCase != -1 ) seqLower( 1, &tmpseq );
			seqlen = strlen( tmpseq );
			lpos = 0;
			for( j=0; j<5; j++ )
			{
				if( regtable[0][j*2] == -1 && regtable[0][j*2+1] == -1 ) continue;

				startpos = regtable[0][j*2];
				endpos   = regtable[0][j*2+1];
				if( startpos > endpos )
				{
					endpos   = regtable[0][j*2];
					startpos = regtable[0][j*2+1];
				}

				if( startpos < 0 ) startpos = 0;
				if( endpos   < 0 ) endpos   = 0;
				if( endpos   >= seqlen ) endpos   = seqlen-1;
				if( startpos >= seqlen ) startpos = seqlen-1;

//				fprintf( stderr, "startpos = %d, endpos = %d\n", startpos, endpos );

				outlen = endpos - startpos+1;
				if( revtable[0][j] == 'f' )
				{
//					fprintf( stderr, "regtable[%d][st] = %d\n", i, regtable[0][j*2+0] );
//					fprintf( stderr, "regtable[%d][en] = %d\n", i, regtable[0][j*2+1] );
//					fprintf( stderr, "outlen = %d\n", outlen );
//					fprintf( stdout, "%.*s\n", outlen, tmpseq+regtable[0][j*2] );
					strncpy( outseq[npos] + lpos, tmpseq+startpos, outlen );
					lpos += outlen;
				}
				else
				{
					fs = AllocateCharVec( outlen+1 );
					rs = AllocateCharVec( outlen+1 );

					fs[outlen] = 0;
					strncpy( fs, tmpseq+startpos, outlen );
					sreverse( rs, fs );
//					fprintf( stdout, "%s\n", rs );
					strncpy( outseq[npos] + lpos, rs, outlen );
					lpos += outlen;
					free( fs );
					free( rs );
				}
				outseq[npos][lpos] = 0;
			}
			npos++;
		}
		free( tmpseq );
	}
}

void cutData( FILE *fp, int **regtable, char **revtable, int *outtable )
{
	CALLS && printf("called %s:cutData()\n", __FILE__);
	int i, j; 
	int outlen, seqlen, startpos, endpos;
	static char *tmpseq = NULL;
	static char *dumname = NULL;
	char *fs, *rs;

	if( dumname == NULL )
	{
		dumname = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		dumname[0] = '='; getc( fp ); 
		myfgets( dumname+1, B-2, fp ); 
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );

		if( outtable[i] )
		{
			gappick_samestring( tmpseq );
			putc( '>', stdout );
			puts( dumname+1 );

			seqlen = strlen( tmpseq );

			if( dorp == 'd' && upperCase != -1 ) seqLower( 1, &tmpseq );
			if( outtable[i] == 2 )
			{
				startpos = 0;
				endpos   = seqlen-1;
				outlen = endpos - startpos + 1;
				fprintf( stdout, "%.*s\n", outlen, tmpseq+startpos );
			}
			else
			{
				for( j=0; j<5; j++ )
				{
					if( regtable[i][j*2] == -1 && regtable[i][j*2+1] == -1 ) continue;
	
					startpos = regtable[i][j*2];
					endpos   = regtable[i][j*2+1];
	
					if( startpos > endpos )
					{
						endpos   = regtable[i][j*2];
						startpos = regtable[i][j*2+1];
					}
	
					if( startpos < 0 ) startpos = 0;
					if( endpos   < 0 ) endpos   = 0;
					if( endpos   >= seqlen ) endpos   = seqlen-1;
					if( startpos >= seqlen ) startpos = seqlen-1;
	
					outlen = endpos - startpos + 1;
					if( revtable[i][j] == 'f' )
					{
						fprintf( stderr, "startpos = %d\n", startpos );
						fprintf( stderr, "endpos   = %d\n", endpos );
						fprintf( stderr, "outlen = %d\n", outlen );
						fprintf( stdout, "%.*s\n", outlen, tmpseq+startpos );
					}
					else
					{
						fs = AllocateCharVec( outlen+1 );
						rs = AllocateCharVec( outlen+1 );
	
						fs[outlen] = 0;
						strncpy( fs, tmpseq+startpos, outlen );
						sreverse( rs, fs );
						fprintf( stdout, "%s\n", rs );
						free( fs );
						free( rs );
					}
				}
			}
		}
		free( tmpseq );
	}
}

void catData( FILE *fp )
{
	CALLS && printf("called %s:catData()\n", __FILE__);
	int i; 
	static char *tmpseq = NULL;
	static char *dumname = NULL;
//	char *cptr;

	if( dumname == NULL )
	{
		dumname = AllocateCharVec( N );
	}

	rewind( fp );
	searchKUorWA( fp );

	for( i=0; i<njob; i++ )
	{
		dumname[0] = '='; getc( fp ); 
		myfgets( dumname+1, B-2, fp ); 
		if( outnumber )
		{
			fprintf( stdout, ">_numo_s_%08d_numo_e_", i+1 );
		}
		else
		{
			putc( '>', stdout );
		}
		puts( dumname+1 );
		tmpseq = load1SeqWithoutName_realloc( fp );
		if( dorp == 'd' && upperCase != -1 ) seqLower( 1, &tmpseq );
		puts( tmpseq );
		free( tmpseq );
	}
}

int countATGC( char *s, int *total ) // yes
{
	CALLS && printf("called %s:countATGC()\n", __FILE__);
	int nATGC;
	int nChar;
	char c;
	nATGC = nChar = 0;

	if( *s == 0 ) 
	{
		*total = 0;
		return( 0 );
	}

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
		}
	}
	while( *++s );

	*total = nChar;
	return( nATGC );
}

double countATGCbk( char *s )
{
	CALLS && printf("called %s:countATGCbk()\n", __FILE__);
	int nATGC;
	int nChar;
	char c;
	nATGC = nChar = 0;

	do
	{
		c = tolower( *s );
		if( isalpha( c ) )
		{
			nChar++;
			if( c == 'a' || c == 't' || c == 'g' || c == 'c' || c == 'u' || c == 'n' )
				nATGC++;
		}
	}
	while( *++s );
	return( (double)nATGC / nChar );
}


int countnogaplen( char *seq )
{
	CALLS && printf("called %s:countnogaplen()\n", __FILE__);
	int val = 0;
	while( *seq )
		if( *seq++ != '-' ) val++;
	return( val );
}

int countnormalletters( char *seq, char *ref )
{
	CALLS && printf("called %s:countnormalletters()\n", __FILE__);
	int val = 0;
	while( *seq )
		if( strchr( ref, *seq++ ) ) val++;
	return( val );
}

void getnumlen_casepreserve( FILE *fp, int *nlenminpt )
{
	CALLS && printf("called %s:getnumlen_casepreserve()\n", __FILE__);
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq, *tmpname;
	double atgcfreq;
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		tmp = strlen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		atgcnum += countATGC( tmpseq, &nsite );
		total += nsite;
		free( tmpseq );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
//	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
}

void getnumlen_nogap( FILE *fp, int *nlenminpt )
{
	CALLS && printf("called %s:getnumlen_nogap()\n", __FILE__);
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq, *tmpname;
	double atgcfreq;
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpseq = load1SeqWithoutName_realloc( fp );
		tmp = countnogaplen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		atgcnum += countATGC( tmpseq, &nsite );
		total += nsite;
		free( tmpseq );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
}


void getnumlen_nogap_outallreg( FILE *fp, int *nlenminpt )
{
	CALLS && printf("called %s:getnumlen_nogap_outallreg()\n", __FILE__);
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq, *tmpname;
	double atgcfreq;
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		fprintf( stdout, "%s\n", tmpname );
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		tmp = countnogaplen( tmpseq );
		fprintf( stdout, "%d\n", tmp );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		atgcnum += countATGC( tmpseq, &nsite );
		total += nsite;
		free( tmpseq );
	}
	free( tmpname );
	atgcfreq = (double)atgcnum / total;
	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
}

static void escapehtml( char *res, char *ori, int maxlen )
{
	CALLS && printf("called %s:escapehtml()\n", __FILE__);
	char *res0 = res;
	while( *ori )
	{
		if( *ori == '<' ) 
		{
			strcpy( res, "&lt;" );
			res += 3;
		}
		else if( *ori == '>' ) 
		{
			strcpy( res, "&gt;" );
			res += 3;
		}
		else if( *ori == '&' ) 
		{
			strcpy( res, "&amp;" );
			res += 4;
		}
		else if( *ori == '"' ) 
		{
			strcpy( res, "&quot;" );
			res += 5;
		}
		else if( *ori == ' ' ) 
		{
			strcpy( res, "&nbsp;" );
			res += 5;
		}
		else
		{
			*res = *ori;
		}
		res++;
		ori++;

		if( res - res0 -10 > N ) break;
	}
	*res = 0;
}

void getnumlen_nogap_outallreg_web( FILE *fp, FILE *ofp, int *nlenminpt, int *isalignedpt )
{
	CALLS && printf("called %s:getnumlen_nogap_outallreg_web()\n", __FILE__);
	int total;
	int nsite = 0;
	int atgcnum;
	int alnlen = 0, alnlen_prev;
	int i, tmp, lennormalchar;
	char *tmpseq, *tmpname, *tmpname2;
	double atgcfreq;
	tmpname = AllocateCharVec( N );
	tmpname2 = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	*nlenminpt = 99999999;
	atgcnum = 0;
	total = 0;
	alnlen_prev = -1;
	*isalignedpt = 1;
	for( i=0; i<njob; i++ )
	{
		myfgets( tmpname, N-1, fp );
		tmpname2[0] = tmpname[0];
		escapehtml( tmpname2+1, tmpname+1, N );
//		fprintf( stdout, "%s\n", tmpname );
//		fprintf( stdout, "%s\n", tmpname2 );
//		exit(1);
		tmpseq = load1SeqWithoutName_realloc_casepreserve( fp );
		tmp = countnogaplen( tmpseq );
//		fprintf( stdout, "%d\n", tmp );
		if( tmp > nlenmax ) nlenmax  = tmp;
		if( tmp < *nlenminpt ) *nlenminpt  = tmp;
		atgcnum += countATGC( tmpseq, &nsite );
		total += nsite;

		alnlen = strlen( tmpseq );
//		fprintf( stdout, "##### alnlen, alnlen_prev = %d, %d\n", alnlen, alnlen_prev );
		if( i>0 && alnlen_prev != alnlen ) *isalignedpt = 0;
		alnlen_prev = alnlen;

		atgcfreq = (double)atgcnum / total;
//		fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
//		if( dorp == NOTSPECIFIED ) // you kentou
		{
			if( atgcfreq > 0.75 ) 	
			{
				dorp = 'd';
				upperCase = -1;
			}
			else                  
			{
				dorp = 'p';
				upperCase = 0;
			}
		}
	
		if( dorp == 'd' ) lennormalchar = countnormalletters( tmpseq, "atgcuATGCU" );
		else              lennormalchar = countnormalletters( tmpseq, "ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv" );
		free( tmpseq );

		fprintf( ofp, " <label for='s%d'><span id='ss%d'><input type='checkbox' id='s%d' name='s%d' checked></span> <input type='text' class='ll' id='ll%d' style='display:none' size='6' value='%d' readonly='readonly'>%s</label>\n", i, i, i, i, i, lennormalchar, tmpname2 );
		fprintf( ofp, "<span id='t%d-0' style='display:none'>", i );
		fprintf( ofp, " <a href='javascript:void(0)' onclick='ddcycle(this.form,\"t%d\")'>+reg</a>", i );
		fprintf( ofp, " Begin:<input type='text' name='b%d-0' size='8' value='1' class='ie'> End:<input type='text' name='e%d-0' size='8' value='%d' class='ie'>", i, i, tmp );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-0'><input type='checkbox' name='r%d-0' id='r%d-0'>Reverse</label>", i, i, i );
//		fprintf( ofp, "  Sequence Length:<input type='text' name='l%d' size='8' value='%d' readonly='readonly'>", i, tmp );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-1' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-1' size='8' value='' class='ie'> End:<input type='text' name='e%d-1' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-1'><input type='checkbox' name='r%d-1' id='r%d-1'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-2' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-2' size='8' value='' class='ie'> End:<input type='text' name='e%d-2' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-2'><input type='checkbox' name='r%d-2' id='r%d-2'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-3' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-3' size='8' value='' class='ie'> End:<input type='text' name='e%d-3' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-3'><input type='checkbox' name='r%d-3' id='r%d-3'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='t%d-4' style='display:none'>", i );
		fprintf( ofp, "      Begin:<input type='text' name='b%d-4' size='8' value='' class='ie'> End:<input type='text' name='e%d-4' size='8' value='' class='ie'>", i, i );
		if( dorp == 'd' ) fprintf( ofp, " <label for='r%d-4'><input type='checkbox' name='r%d-4' id='r%d-4'>Reverse</label>", i, i, i );
		fprintf( ofp, "\n</span>" );
	}
	free( tmpname );
	free( tmpname2 );
	atgcfreq = (double)atgcnum / total;
	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
//	if( dorp == NOTSPECIFIED ) // you kentou
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
	fprintf( ofp, "\n" );
	if( *isalignedpt )
	{
		fprintf( ofp, "<span id='tall-0' style='display:none'>" );
		fprintf( ofp, "Cut the alignment\n" );
		fprintf( ofp, " <a href='javascript:void(0)' onclick='ddcycle(this.form,\"tall\")'>+reg</a>" );
		fprintf( ofp, " Begin:<input type='text' name='ball-0' size='8' value='1'> End:<input type='text' name='eall-0' size='8' value='%d'>", alnlen );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-0'><input type='checkbox' name='rall-0' id='rall-0'>Reverse</label>" );
		fprintf( ofp, "  Alignment length:<input type='text' name='lall' size='8' value='%d' readonly='readonly'>", alnlen );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-1' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-1' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-1'><input type='checkbox' name='rall-1' id='rall-1'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-2' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-2' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-2'><input type='checkbox' name='rall-2' id='rall-2'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-3' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-3' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-3'><input type='checkbox' name='rall-3' id='rall-3'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
		fprintf( ofp, "<span id='tall-4' style='display:none'>" );
		fprintf( ofp, "      Begin:<input type='text' name='ball-4' size='8' value=''> End:<input type='text' name='eall-1' size='8' value=''>" );
		if( dorp == 'd' ) fprintf( ofp, " <label for='rall-4'><input type='checkbox' name='rall-4' id='rall-4'>Reverse</label>" );
		fprintf( ofp, "\n</span>" );
	}

}

void getnumlen( FILE *fp ) // yes
{
	CALLS && printf("called %s:getnumlen()\n", __FILE__);
	int total;
	int nsite = 0;
	int atgcnum;
	int i, tmp;
	char *tmpseq;
	char *tmpname;
	double atgcfreq;
	tmpname = AllocateCharVec( N );
	njob = countKUorWA( fp );
	searchKUorWA( fp );
	nlenmax = 0;
	atgcnum = 0;
	total = 0;
	for( i=0; i<njob; i++ )
	{
		FILES && printf("file read %d %s:%d\n", fp, __FILE__, __LINE__);
		myfgets( tmpname, N-1, fp );
		tmpseq = load1SeqWithoutName_realloc( fp );
		tmp = strlen( tmpseq );
		if( tmp > nlenmax ) nlenmax  = tmp;
		atgcnum += countATGC( tmpseq, &nsite );
		total += nsite;
//		fprintf( stderr, "##### total = %d\n", total );
		free( tmpseq );
	}


	atgcfreq = (double)atgcnum / total;
//	fprintf( stderr, "##### atgcfreq = %f\n", atgcfreq );
	if( dorp == NOTSPECIFIED )
	{
		if( atgcfreq > 0.75 ) 	
		{
			dorp = 'd';
			upperCase = -1;
		}
		else                  
		{
			dorp = 'p';
			upperCase = 0;
		}
	}
	free( tmpname );
}
	


void WriteGapFill( FILE *fp, int locnjob, char name[][B], int nlen[M], char **aseq )
{
	CALLS && printf("called %s:WriteGapFill()\n", __FILE__);
	static char b[N];
	int i, j;
	int nalen[M];
	static char gap[N];
	static char buff[N];

#if IODEBUG
	fprintf( stderr, "IMAKARA KAKU\n" );
#endif
	nlenmax = 0;
	for( i=0; i<locnjob; i++ )
	{
		int len = strlen( aseq[i] );
		if( nlenmax < len ) nlenmax = len;
	}

	for( i=0; i<nlenmax; i++ ) gap[i] = '-';
	gap[nlenmax] = 0;

	fprintf( fp, "%5d", locnjob );
	fprintf( fp, "\n" );

	for( i=0; i<locnjob; i++ )
	{
		strcpy( buff, aseq[i] );
		strncat( buff, gap, nlenmax-strlen( aseq[i] ) );
		buff[nlenmax] = 0;
		nalen[i] = strlen( buff );
		fprintf( fp, "%s\n", name[i] );
		fprintf( fp, "%5d\n", nalen[i] );
		for( j=0; j<nalen[i]; j=j+C )
		{
			strncpy_caseC( b, buff+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
		}
	}
#if DEBUG
	fprintf( stderr, "nalen[0] = %d\n", nalen[0] );
#endif
#if IODEBUG
	fprintf( stderr, "KAKIOWATTA\n" );
#endif
}

void writeDataforgaln( FILE *fp, int locnjob, char **name, int *nlen, char **aseq )
{
	CALLS && printf("called %s:writeDataforgaln()\n", __FILE__);
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		nalen = strlen( aseq[i] );
		fprintf( fp, ">%s\n", name[i]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[i]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[i]+j );
#endif
		}
	}
}

void writeData_pointer( FILE *fp, int locnjob, char **name, int *nlen, char **aseq ) // yes
{
	CALLS && printf("called %s:writeData_pointer()\n", __FILE__);
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
#if DEBUG
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		nalen = strlen( aseq[i] );
		FILES && printf("file write %d %s:%d\n", fp, __FILE__, __LINE__);
		fprintf( fp, ">%s\n", name[i]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[i]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			FILES && printf("file write %d %s:%d\n", fp, __FILE__, __LINE__);
			fprintf( fp, "%.*s\n", C, aseq[i]+j );
#endif
		}
	}
}

void writeData( FILE *fp, int locnjob, char name[][B], int nlen[], char **aseq )
{
	CALLS && printf("called %s:writeData()\n", __FILE__);
	int i, j;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
#if DEBUG
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		nalen = strlen( aseq[i] );
		fprintf( fp, ">%s\n", name[i]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[i]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[i]+j );
#endif
		}
	}
}


void write1seq( FILE *fp, char *aseq )
{
	CALLS && printf("called %s:write1seq()\n", __FILE__);
	int j;
	int nalen;

	nalen = strlen( aseq );
	for( j=0; j<nalen; j=j+C )
		fprintf( fp, "%.*s\n", C, aseq+j );
}

void readhat2_floathalf_part_pointer( FILE *fp, int nseq, int nadd, char **name, float **mtx )
{
	CALLS && printf("called %s:readhat2_floathalf_part_pointer()\n", __FILE__);
    int i, j, nseq0, norg;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) 
	{
		fprintf( stderr, "%d != %d\n", nseq, nseq0 );
		ErrorExit( "hat2 is wrong." );
	}
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
	norg = nseq-nadd;
    for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ )
    {
        mtx[i][j] = ( input_new( fp, D ) );
    }
}

void readhat2_floathalf_pointer( FILE *fp, int nseq, char **name, float **mtx )
{
	CALLS && printf("called %s:readhat2_floathalf_pointer()\n", __FILE__);
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) 
	{
		fprintf( stderr, "%d != %d\n", nseq, nseq0 );
		ErrorExit( "hat2 is wrong." );
	}
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j-i] = ( input_new( fp, D ) );
    }
}
void readhat2_floathalf( FILE *fp, int nseq, char name[M][B], float **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j-i] = ( input_new( fp, D ) );
    }
}
void readhat2_float( FILE *fp, int nseq, char name[M][B], float **mtx )
{
	CALLS && printf("called %s:readhat2_float()\n", __FILE__);
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = ( input_new( fp, D ) );
    }
}
void readhat2_int( FILE *fp, int nseq, char name[M][B], int **mtx )
{
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = (int)( input_new( fp, D ) * INTMTXSCALE + 0.5 );
    }
}

void readhat2_pointer( FILE *fp, int nseq, char **name, double **mtx )
{
	CALLS && printf("called %s:readhat2_pointer()\n", __FILE__);
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = (double)input_new( fp, D);
    }
}
void readhat2( FILE *fp, int nseq, char name[M][B], double **mtx )
{
	CALLS && printf("called %s:readhat2()\n", __FILE__);
    int i, j, nseq0;
    char b[B];

    fgets( b, B, fp );
    fgets( b, B, fp ); b[5] = 0; nseq0 = atoi( b ); if( nseq != nseq0 ) ErrorExit( "hat2 is wrong." );
    fgets( b, B, fp );
    for( i=0; i<nseq; i++ )
    {
#if 0
        getaline_fp_eof( b, B, fp ); 
#else
		myfgets( b, B-2, fp );
#endif
#if 0
		j = MIN( strlen( b+6 ), 10 );
        if( strncmp( name[i], b+6 , j ) ) 
		{
			fprintf( stderr, "Error in hat2\n" );
			fprintf( stderr, "%s != %s\n", b, name[i] );
			exit( 1 );
		}
#endif
    }
    for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
    {
        mtx[i][j] = (double)input_new( fp, D);
    }
}

void WriteFloatHat2_pointer_halfmtx( FILE *hat2p, int locnjob, char **name, float **mtx )
{
	CALLS && printf("called %s:WriteFloatHat2_pointer_halfmtx()\n", __FILE__);
	int i, j, ijsa;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=1; j<locnjob-i; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob; i++ )
	{
		for( j=i+1; j<njob; j++ )
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j-i] );
			ijsa = j-i;
			if( ijsa % 12 == 0 || ijsa == locnjob-i-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteFloatHat2_pointer( FILE *hat2p, int locnjob, char **name, float **mtx )
{
	CALLS && printf("called %s:WriteFloatHat2_pointer()\n", __FILE__);
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=1; j<locnjob-i; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob; i++ )
	{
		for( j=1; j<locnjob-i; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j] );
			if( j % 12 == 0 || j == locnjob-i-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteFloatHat2( FILE *hat2p, int locnjob, char name[M][B], float **mtx )
{
	CALLS && printf("called %s:WriteFloatHat2()\n", __FILE__);
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=1; j<locnjob-i; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob; i++ )
	{
		for( j=1; j<locnjob-i; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j] );
			if( j % 12 == 0 || j == locnjob-i-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2_int( FILE *hat2p, int locnjob, char name[M][B], int **mtx )
{
	CALLS && printf("called %s:WriteHat2_int()\n", __FILE__);
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];
	max /= INTMTXSCALE;

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", (float)mtx[i][j] / INTMTXSCALE );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2_part_pointer( FILE *hat2p, int locnjob, int nadd, char **name, double **mtx )
{
	CALLS && printf("called %s:WriteHat2_part_pointer()\n", __FILE__);
	int i, j;
	int norg = locnjob-nadd;
	double max = 0.0;
//	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<norg; i++ )
	{
		for( j=0; j<nadd; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j] );
			if( (j+1) % 12 == 0 || j == nadd-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2_pointer( FILE *hat2p, int locnjob, char **name, double **mtx )
{
	CALLS && printf("called %s:WriteHat2_pointer()\n", __FILE__);
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j] );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

void WriteHat2( FILE *hat2p, int locnjob, char name[M][B], double **mtx )
{
	CALLS && printf("called %s:WriteHat2()\n", __FILE__);
	int i, j;
	double max = 0.0;
	for( i=0; i<locnjob-1; i++ ) for( j=i+1; j<locnjob; j++ ) if( mtx[i][j] > max ) max = mtx[i][j];

	fprintf( hat2p, "%5d\n", 1 );
	fprintf( hat2p, "%5d\n", locnjob );
	fprintf( hat2p, " %#6.3f\n", max * 2.5 );

	for( i=0; i<locnjob; i++ ) fprintf( hat2p, "%4d. %s\n", i+1, name[i] );
	for( i=0; i<locnjob-1; i++ )
	{
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, "%#6.3f", mtx[i][j] );
			if( (j-i) % 12 == 0 || j == locnjob-1 ) fprintf( hat2p, "\n" );
		}
	}
}

#if 0
void WriteHat2plain( FILE *hat2p, int locnjob, double **mtx )
{
	int i, j, ilim;

	ilim = locnjob-1;
	for( i=0; i<ilim; i++ )
	{
		fprintf( hat2p, "%d-%d d=%.3f\n", i+1, i+1, 0.0 );
		for( j=i+1; j<locnjob; j++ ) 
		{
			fprintf( hat2p, "%d-%d d=%.3f\n", i+1, j+1, mtx[i][j] );
		}
	}
}
#endif

int ReadFasta_sub( FILE *fp, double *dis, int nseq, char name[M][B] )
{
	CALLS && printf("called %s:ReadFasta_sub()\n", __FILE__);
    int i, count=0;
    char b[B];
    int junban[M];

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
            count++;
        }
    }

	for( i=0; i<nseq; i++ ) dis[i] = 0.0;
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
		if( fgets( b, B-1, fp ) ) break;
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            fgets( b, B-1, fp );
            dis[junban[count]] = atof( b );
            count++;
        }
    }
    return 0;
}


int ReadSsearch( FILE *fp, double *dis, int nseq, char name[M][B] )
{
	CALLS && printf("called %s:ReadSsearch()\n", __FILE__);
    int i, count=0;
    char b[B];
    int junban[M];
	int opt;

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
			sscanf( b+75, "%d", &opt ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }

/*
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            dis[junban[count]] = atof( b+65 );
            count++;
        }
    }
*/
    return 0;
}

int ReadBlastm7_avscore( FILE *fp, double *dis, int nin )
{
	CALLS && printf("called %s:ReadBlastm7_avscore()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    int *junban;
	double score, sumscore;
	double len, sumlen;
	int qstart, qend, tstart, tend;
	double scorepersite;
	static char qal[N], tal[N], al[N];
	int nlocalhom;

	junban = calloc( nin, sizeof( int ) );

	count = 0;
	sumscore = 0.0;
	sumlen = 0.0;
	score = 0.0;
	len = 0.0;
	scorepersite = 0.0; // by D.Mathog, a guess
    while( 1 )
	{

		if( feof( fp ) ) break;

		while( fgets( b, B-1, fp ) )
		{
			if( !strncmp( "          <Hit_def>", b, 19 ) || !strncmp( "              <Hsp_num>", b, 23 ) ) break;
		}

		if( !strncmp( "          <Hit_def>", b, 19 ) )
		{
			junban[count] = atoi( b+31 );
			nlocalhom = 0;
		}


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_score>", b, 25 ) ) break;
		pt = b + 25;
		score = atof( pt );
		sumscore += score;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-from>", b, 30 ) ) break;
		pt = b + 30;
		qstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-to>", b, 28 ) ) break;
		pt = b + 28;
		qend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-from>", b, 28 ) ) break;
		pt = b + 28;
		tstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-to>", b, 26 ) ) break;
		pt = b + 26;
		tend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_align-len>", b, 29 ) ) break;
		pt = b + 29;
		len = atoi( pt );
		sumlen += len;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_qseq>", al, 24 ) ) break;

		strcpy( qal, al+24 );
		pt = qal;
		while( *++pt != '<' )
			;
		*pt = 0;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_hseq>", al, 24 ) ) break;

		strcpy( tal, al+24 );
		pt = tal;
		while( *++pt != '<' )
			;
		*pt = 0;


//		fprintf( stderr, "t=%d, score = %f, qstart=%d, qend=%d, tstart=%d, tend=%d, overlapaa=%d\n", junban[count], score, qstart, qend, tstart, tend, overlapaa );


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "            </Hsp>:", b, 18 ) ) break;


		fgets( b, B-1, fp );


		if( !strncmp( "          </Hit_hsps>", b, 21 ) )
		{
			dis[junban[count++]] = sumscore;
			sumscore = 0.0;
			fgets( b, B-1, fp );
			fgets( b, B-1, fp );
			scorepersite = sumscore / sumlen;
			if( scorepersite != (int)scorepersite )
			{
				fprintf( stderr, "ERROR! sumscore=%f, sumlen=%f, and scorepersite=%f\n", sumscore, sumlen, scorepersite );
				exit( 1 );
			}

			if( !strncmp( "      </Iteration_hits>", b, 23 ) ) break;
		}
	}

	free( junban );

    return (int)scorepersite;
}
int ReadBlastm7_scoreonly( FILE *fp, double *dis, int nin )
{
	CALLS && printf("called %s:ReadBlastm7_scoreonly()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    int *junban;
	int overlapaa;
	double score, sumscore;
	int qstart, qend, tstart, tend;
	static char qal[N], tal[N], al[N];
	int nlocalhom;

	junban = calloc( nin, sizeof( int ) );

	count = 0;
	sumscore = 0.0;
	score = 0.0;
    while( 1 )
	{

		if( feof( fp ) ) break;

		while( fgets( b, B-1, fp ) )
		{
			if( !strncmp( "          <Hit_def>", b, 19 ) || !strncmp( "              <Hsp_num>", b, 23 ) ) break;
		}

		if( !strncmp( "          <Hit_def>", b, 19 ) )
		{
			junban[count] = atoi( b+31 );
			nlocalhom = 0;
		}


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_score>", b, 25 ) ) break;
		pt = b + 25;
		score = atof( pt );
		sumscore += score;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-from>", b, 30 ) ) break;
		pt = b + 30;
		qstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-to>", b, 28 ) ) break;
		pt = b + 28;
		qend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-from>", b, 28 ) ) break;
		pt = b + 28;
		tstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-to>", b, 26 ) ) break;
		pt = b + 26;
		tend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_align-len>", b, 29 ) ) break;
		pt = b + 29;
		overlapaa = atoi( pt );


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_qseq>", al, 24 ) ) break;

		strcpy( qal, al+24 );
		pt = qal;
		while( *++pt != '<' )
			;
		*pt = 0;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_hseq>", al, 24 ) ) break;

		strcpy( tal, al+24 );
		pt = tal;
		while( *++pt != '<' )
			;
		*pt = 0;


//		fprintf( stderr, "t=%d, score = %f, qstart=%d, qend=%d, tstart=%d, tend=%d, overlapaa=%d\n", junban[count], score, qstart, qend, tstart, tend, overlapaa );

//		nlocalhom += addlocalhom_r( qal, tal, localhomlist+junban[count], qstart, tstart, score, overlapaa, nlocalhom );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "            </Hsp>:", b, 18 ) ) break;


		fgets( b, B-1, fp );


		if( !strncmp( "          </Hit_hsps>", b, 21 ) )
		{
			dis[junban[count++]] = sumscore;
			sumscore = 0.0;
			fgets( b, B-1, fp );
			fgets( b, B-1, fp );
			if( !strncmp( "      </Iteration_hits>", b, 23 ) ) break;
		}
	}

	free( junban );

    return count;
}

int ReadBlastm7( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
	CALLS && printf("called %s:ReadBlastm7()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	double score, sumscore;
	int qstart, qend, tstart, tend;
	static char qal[N], tal[N], al[N];
	int nlocalhom;



	count = 0;
	sumscore = 0.0;
	score = 0.0;
	nlocalhom = 0;
    while( 1 )
	{

		if( feof( fp ) ) break;

		while( fgets( b, B-1, fp ) )
		{
			if( !strncmp( "          <Hit_def>", b, 19 ) || !strncmp( "              <Hsp_num>", b, 23 ) ) break;
		}

		if( !strncmp( "          <Hit_def>", b, 19 ) )
		{
			junban[count] = atoi( b+31 );
			nlocalhom = 0;
		}


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_score>", b, 25 ) ) break;
		pt = b + 25;
		score = atof( pt );
		sumscore += score;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-from>", b, 30 ) ) break;
		pt = b + 30;
		qstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_query-to>", b, 28 ) ) break;
		pt = b + 28;
		qend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-from>", b, 28 ) ) break;
		pt = b + 28;
		tstart = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_hit-to>", b, 26 ) ) break;
		pt = b + 26;
		tend = atoi( pt ) - 1;


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "              <Hsp_align-len>", b, 29 ) ) break;
		pt = b + 29;
		overlapaa = atoi( pt );


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_qseq>", al, 24 ) ) break;

		strcpy( qal, al+24 );
		pt = qal;
		while( *++pt != '<' )
			;
		*pt = 0;


		while( fgets( al, N-100, fp ) )
			if( !strncmp( "              <Hsp_hseq>", al, 24 ) ) break;

		strcpy( tal, al+24 );
		pt = tal;
		while( *++pt != '<' )
			;
		*pt = 0;


//		fprintf( stderr, "t=%d, score = %f, qstart=%d, qend=%d, tstart=%d, tend=%d, overlapaa=%d\n", junban[count], score, qstart, qend, tstart, tend, overlapaa );

		nlocalhom += addlocalhom_r( qal, tal, localhomlist+junban[count], qstart, tstart, score, overlapaa, nlocalhom );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "            </Hsp>:", b, 18 ) ) break;


		fgets( b, B-1, fp );


		if( !strncmp( "          </Hit_hsps>", b, 21 ) )
		{
			dis[junban[count++]] = sumscore;
			sumscore = 0.0;
			fgets( b, B-1, fp );
			fgets( b, B-1, fp );
			if( !strncmp( "      </Iteration_hits>", b, 23 ) ) break;
		}
	}
    return count;
}

int ReadFasta34noalign( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
	CALLS && printf("called %s:ReadFasta34noalign()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int opt;
	double z, bits;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
    }

    return count;
}
int ReadFasta34m10_nuc( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
	CALLS && printf("called %s:ReadFasta34m10_nuc()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	int opt, qstart, qend, tstart, tend;
	double z, bits;
	int qal_display_start, tal_display_start;
	static char qal[N], tal[N];
	char *qal2, *tal2;
	int c;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			if( strchr( b, 'r' ) ) continue;

			pt = strchr( b, ']' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
		else if( 0 == strncmp( ">>+==========+", b, 14 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

	count = 0;
    while( 1 )
	{
		if( strncmp( ">>+==========+", b, 14 ) )
		{
			fgets( b, B-1, fp );
			if( feof( fp ) ) break;
			continue;
		}
		junban[count++] = atoi( b+14 );
//		fprintf( stderr, "t = %d\n", atoi( b+14 ) );
		while( fgets( b, B-1, fp ) )
			if( !strncmp( "; fa_opt:", b, 9 ) || !strncmp( "; sw_s-w opt:", b, 13 ) ) break;
		pt = strstr( b, ":" ) +1;
		opt = atoi( pt );


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_overlap:", b+4, 9 ) ) break;
		pt = strstr( b, ":" ) +1;
		overlapaa = atoi( pt );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) +1;
		qstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) +1;
		qend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) +1;
		qal_display_start = atoi( pt ) - 1;

		pt = qal;
		while( (c = fgetc( fp )) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tal_display_start = atoi( pt ) - 1;

		pt = tal;
		while( ( c = fgetc( fp ) ) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

//		fprintf( stderr, "(%d-%d:%d-%d)\n", qstart, qend, tstart, tend );
//		fprintf( stderr, "qal_display_start = %d, tal_display_start = %d\n", qal_display_start, tal_display_start );

//		fprintf( stderr, "qal = %s\n", qal );
//		fprintf( stderr, "tal = %s\n", tal );

		qal2 = cutal( qal, qal_display_start, qstart, qend );
		tal2 = cutal( tal, tal_display_start, tstart, tend );

//		fprintf( stderr, "qal2 = %s\n", qal2 );
//		fprintf( stderr, "tal2 = %s\n", tal2 );

//		fprintf( stderr, "putting   %d - %d, opt = %d\n", qmem, junban[count-1], opt );
		putlocalhom( qal2, tal2, localhomlist+junban[count-1], qstart, tstart, opt, overlapaa );
	}
//	fprintf( stderr, "count = %d\n", count );
    return count;
}
int ReadFasta34m10( FILE *fp, double *dis, int qmem, char **name, LocalHom *localhomlist )
{
	CALLS && printf("called %s:ReadFasta34m10()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	int opt, qstart, qend, tstart, tend;
	double z, bits;
	int qal_display_start, tal_display_start;
	static char qal[N], tal[N];
	char *qal2, *tal2;
	int c;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
		else if( 0 == strncmp( ">>+==========+", b, 14 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

	count = 0;
    while( 1 )
	{
		if( strncmp( ">>+==========+", b, 14 ) )
		{
			fgets( b, B-1, fp );
			if( feof( fp ) ) break;
			continue;
		}
		junban[count++] = atoi( b+14 );
//		fprintf( stderr, "t = %d\n", atoi( b+14 ) );
		while( fgets( b, B-1, fp ) )
			if( !strncmp( "; fa_opt:", b, 9 ) || !strncmp( "; sw_s-w opt:", b, 13 ) ) break;
		pt = strstr( b, ":" ) +1;
		opt = atoi( pt );


		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_overlap:", b+4, 9 ) ) break;
		pt = strstr( b, ":" ) +1;
		overlapaa = atoi( pt );

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) +1;
		qstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) +1;
		qend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) +1;
		qal_display_start = atoi( pt ) - 1;

		pt = qal;
		while( (c = fgetc( fp )) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_start:", b+4, 7 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tstart = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_stop:", b+4, 6 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tend = atoi( pt ) - 1;

		while( fgets( b, B-1, fp ) )
			if( !strncmp( "_display_start:", b+4, 15 ) ) break;
		pt = strstr( b, ":" ) + 1;
		tal_display_start = atoi( pt ) - 1;

		pt = tal;
		while( ( c = fgetc( fp ) ) )
		{
			if( c == '>' ) 
			{
				ungetc( c, fp );
				break;
			}
			if( isalpha( c ) || c == '-' ) 
			*pt++ = c;
		}
		*pt = 0;

//		fprintf( stderr, "(%d-%d:%d-%d)\n", qstart, qend, tstart, tend );
//		fprintf( stderr, "qal_display_start = %d, tal_display_start = %d\n", qal_display_start, tal_display_start );

//		fprintf( stderr, "qal = %s\n", qal );
//		fprintf( stderr, "tal = %s\n", tal );

		qal2 = cutal( qal, qal_display_start, qstart, qend );
		tal2 = cutal( tal, tal_display_start, tstart, tend );

//		fprintf( stderr, "qal2 = %s\n", qal2 );
//		fprintf( stderr, "tal2 = %s\n", tal2 );

//		fprintf( stderr, "putting   %d - %d, opt = %d\n", qmem, junban[count-1], opt );
		putlocalhom( qal2, tal2, localhomlist+junban[count-1], qstart, tstart, opt, overlapaa );
	}
//	fprintf( stderr, "count = %d\n", count );
    return count;
}
int ReadFasta34m10_scoreonly_nucbk( FILE *fp, double *dis, int nin )
{
	CALLS && printf("called %s:ReadFasta34m10_scoreonly_nucbk()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    int pos;
	int opt;
	double z, bits;

    count = 0;
    while( !feof( fp ) )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+===========+", b, 13 ) )
        {
            pos = atoi( b+13 );

			if( strchr( b, 'r' ) ) continue;

//			pt = strchr( b, ')' ) + 1;
			pt = strchr( b, ']' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[pos] += (double)opt;
            count++;
#if 0
			fprintf( stderr, "b=%s\n", b );
			fprintf( stderr, "opt=%d\n", opt );
			fprintf( stderr, "pos=%d\n", pos );
			fprintf( stderr, "dis[pos]=%f\n", dis[pos] );
#endif

        }
		else if( 0 == strncmp( ">>><<<", b, 6 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

    return count;
}

int ReadFasta34m10_scoreonly_nuc( FILE *fp, double *dis, int nin )
{
	CALLS && printf("called %s:ReadFasta34m10_scoreonly_nuc()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    int pos;
	int opt;
	double z, bits;
	int c;
	int *yonda;


	yonda = AllocateIntVec( nin );
	for( c=0; c<nin; c++ ) yonda[c] = 0;
	for( c=0; c<nin; c++ ) dis[c] = 0.0;

    count = 0;
    while( !feof( fp ) )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+===========+", b, 13 ) )
        {
            pos = atoi( b+13 );

			if( strchr( b, 'r' ) ) continue;

//			pt = strchr( b, ')' ) + 1;
			pt = strchr( b, ']' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
			if( yonda[pos] == 0 )
			{
	            dis[pos] += (double)opt;
				yonda[pos] = 1;
			}
            count++;
#if 0
			fprintf( stderr, "b=%s\n", b );
			fprintf( stderr, "opt=%d\n", opt );
			fprintf( stderr, "pos=%d\n", pos );
			fprintf( stderr, "dis[pos]=%f\n", dis[pos] );
#endif

        }
        else if( !strncmp( ">>>", b, 3 ) )
		{
			for( c=0; c<nin; c++ ) yonda[c] = 0;
		}
		else if( 0 == strncmp( ">>><<<", b, 6 ) )
		{
			break;
		}

    }

	free( yonda );

	if( !count ) return -1;

    return count;
}

int ReadFasta34m10_scoreonly( FILE *fp, double *dis, int nin )
{
	CALLS && printf("called %s:ReadFasta34m10_scoreonly()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    int pos;
	int opt;
	double z, bits;
	int c;
	int *yonda;


	yonda = AllocateIntVec( nin );
	for( c=0; c<nin; c++ ) yonda[c] = 0;
	for( c=0; c<nin; c++ ) dis[c] = 0.0;

    count = 0;
    while( !feof( fp ) )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+===========+", b, 13 ) )
        {
            pos = atoi( b+13 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
			if( yonda[pos] == 0 )
			{
	            dis[pos] += (double)opt;
				yonda[pos] = 1;
			}
            count++;
#if 0
			fprintf( stderr, "b=%s\n", b );
			fprintf( stderr, "opt=%d\n", opt );
			fprintf( stderr, "pos=%d\n", pos );
			fprintf( stderr, "dis[pos]=%f\n", dis[pos] );
#endif

        }
        else if( !strncmp( ">>>", b, 3 ) )
		{
			for( c=0; c<nin; c++ ) yonda[c] = 0;
		}
		else if( 0 == strncmp( ">>><<<", b, 6 ) )
		{
			break;
		}

    }

	free( yonda );

	if( !count ) return -1;

    return count;
}
int ReadFasta34( FILE *fp, double *dis, int nseq, char name[M][B], LocalHom *localhomlist )
{
	CALLS && printf("called %s:ReadFasta34()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    static int junban[M];
	int overlapaa;
	int opt, qstart, qend, tstart, tend;
	double z, bits;


    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %lf %lf",  &opt, &bits, &z ); 
            dis[junban[count]] = (double)opt;
            count++;

        }
		else if( 0 == strncmp( ">>+==========+", b, 14 ) )
		{
			break;
		}

    }
	if( !count ) return -1;

	count = 0;
    while( !feof( fp ) )
	{
		if( !strncmp(">>+==========+", b, 14 ) )
		{
            junban[count] = atoi( b+14 );
            count++;
        	fgets( b, B-1, fp ); // initn:
			pt = strstr( b, "opt: " ) + 5;
			localhomlist[junban[count-1]].opt = atof( pt );
        	fgets( b, B-1, fp ); // Smith-Waterman score
			pt = strstr( b, "ungapped) in " ) + 13;
			sscanf( pt, "%d", &overlapaa ); 
			fprintf( stderr, "pt = %s, overlapaa = %d\n", pt, overlapaa );
			pt = strstr( b, "overlap (" ) + 8;
			sscanf( pt, "(%d-%d:%d-%d)", &qstart, &qend, &tstart, &tend ); 
			localhomlist[junban[count-1]].overlapaa = overlapaa;
			localhomlist[junban[count-1]].start1 = qstart-1;
			localhomlist[junban[count-1]].end1   = qend-1;
			localhomlist[junban[count-1]].start2 = tstart-1;
			localhomlist[junban[count-1]].end2   = tend-1;
		}
        fgets( b, B-1, fp );
	}
	fprintf( stderr, "count = %d\n", count );
    return count;
}

int ReadFasta3( FILE *fp, double *dis, int nseq, char name[M][B] )
{
	CALLS && printf("called %s:ReadFasta3()\n", __FILE__);
    int count=0;
    char b[B];
	char *pt;
    int junban[M];
	int initn, init1, opt;
	double z;

    count = 0;
#if 0
    for( i=0; i<10000000 && count<nseq; i++ )
#else
    while( !feof( fp ) )
#endif
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			pt = strchr( b, ')' ) + 1;
			sscanf( pt, "%d %d %d %lf", &initn, &init1, &opt, &z ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }
    return 0;
}

int ReadFasta( FILE *fp, double *dis, int nseq, char name[M][B] )
{
	CALLS && printf("called %s:ReadFasta()\n", __FILE__);
    int i, count=0;
    char b[B];
    int junban[M];
	int initn, init1, opt;

    count = 0;
	for( i=0; i<nseq; i++ ) dis[i] = 0.0;
    for( i=0; !feof( fp ) && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );

			sscanf( b+50, "%d %d %d", &initn, &init1, &opt ); 
            dis[junban[count]] = (double)opt;
            count++;
        }
    }

/*
    count = 0;
    for( i=0; i<100000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( name[junban[count]], b, 20  ) )
        {
            dis[junban[count]] = atof( b+65 );
            count++;
        }
    }
*/
    return 0;
}


int ReadOpt( FILE *fp, int opt[M], int nseq, char name[M][B] )
{
	CALLS && printf("called %s:ReadOpt()\n", __FILE__);
    int i, count=0;
    char b[B];
    int junban[M];
	int optt, initn, init1;

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
			sscanf( b+50, "%d %d %d", &initn, &init1, &optt ); 
            opt[junban[count]] = (double)optt;
            count++;
        }
    }
    return 0;
}

int ReadOpt2( FILE *fp, int opt[M], int nseq, char name[M][B] )
{
	CALLS && printf("called %s:ReadOpt2()\n", __FILE__);
    int i, count=0;
    char b[B];
    int junban[M];

    count = 0;
    for( i=0; i<10000000 && count<nseq; i++ )
    {
        fgets( b, B-1, fp );
        if( !strncmp( "+==========+", b, 12 ) )
        {
            junban[count] = atoi( b+12 );
            opt[junban[count]] = atoi( b+65 );
            count++;
        }
    }
    return 0;
}



int writePre( int nseq, char **name, int nlen[M], char **aseq, int force )
{
	CALLS && printf("called %s:writePre()\n", __FILE__);
#if USE_XCED
	int i, value;
	if( !signalSM )
	{
		if( force ) 
		{
			rewind( prep_g );
			if( devide ) dvWrite( prep_g, nseq, name, nlen, aseq );
#if 0
			else    WriteGapFill( prep_g, nseq, name, nlen, aseq );
#else
			else    writeData( prep_g, nseq, name, nlen, aseq );
#endif
		}
		return( 0 );
	}
	for( i=0; i<10; i++ )
	{
#if IODEBUG
		fprintf( stderr, "SEMAPHORE = %d\n", signalSM[SEMAPHORE] );
#endif
		if( signalSM[SEMAPHORE]-- > 0 )
		{
#if 0 /* /tmp/pre の関係ではずした */
			if( ferror( prep_g ) ) prep_g = fopen( "pre", "w" );
			if( !prep_g ) ErrorExit( "Cannot re-open pre." ); 
#endif
			rewind( prep_g );
			signalSM[STATUS] = IMA_KAITERU;
#if IODEBUG
			if( force ) fprintf( stderr, "FINAL " );
#endif
			if( devide ) dvWrite( prep_g, nseq, name, nlen, aseq );
			else    WriteGapFill( prep_g, nseq, name, nlen, aseq );
			/*
			fprintf( prep_g, '\EOF' );
			*/
			fflush( prep_g );
			if( force ) signalSM[STATUS] = OSHIMAI;
			else        signalSM[STATUS] = KAKIOWATTA;
			value = 1;
			signalSM[SEMAPHORE]++;
#if IODEBUG
			fprintf( stderr, "signalSM[STATUS] = %c\n", signalSM[STATUS] );
#endif
			break;
		}
		else
		{
#if IODEBUG
			fprintf( stderr, "YONDERUKARA_AKIRAMERU\n" );
#endif
			value = 0;
			signalSM[SEMAPHORE]++;
			if( !force ) break;
#if IODEBUG
			fprintf( stderr, "MATSU\n" );
#endif
			sleep( 1 );
		}
	}
	if( force && !value ) ErrorExit( "xced ga pre wo hanasanai \n" );
	return( value );
#else
	if( force ) 
	{
		rewind( prep_g );
		writeData_pointer( prep_g, nseq, name, nlen, aseq );
	}
#endif
	return( 0 );
}


void readOtherOptions( int *ppidptr, int *fftThresholdptr, int *fftWinSizeptr )
{
	CALLS && printf("called %s:readOtherOptions()\n", __FILE__);
	if( calledByXced )
	{
		FILE *fp = fopen( "pre", "r" );
		char b[B];
		if( !fp ) ErrorExit( "Cannot open pre.\n" );
		fgets( b, B-1, fp );
		sscanf( b, "%d %d %d", ppidptr, fftThresholdptr, fftWinSizeptr );
		fclose( fp );
#if IODEBUG
	fprintf( stderr, "b = %s\n", b );
	fprintf( stderr, "ppid = %d\n", ppid );
	fprintf( stderr, "fftThreshold = %d\n", fftThreshold );
	fprintf( stderr, "fftWinSize = %d\n", fftWinSize );
#endif
	}
	else
	{
		*ppidptr = 0;
		*fftThresholdptr = FFT_THRESHOLD;
		if( dorp == 'd' )
			*fftWinSizeptr = FFT_WINSIZE_D;
		else
			*fftWinSizeptr = FFT_WINSIZE_P;
	}
#if 0
	fprintf( stderr, "fftThresholdptr=%d\n", *fftThresholdptr );
	fprintf( stderr, "fftWinSizeptr=%d\n", *fftWinSizeptr );
#endif
}

void initSignalSM( void ) // yes
{
	CALLS && printf("called %s:initSignalSM()\n", __FILE__);
//	int signalsmid;

#if IODEBUG
	FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
	if( ppid ) fprintf( stderr, "PID of xced = %d\n", ppid );
#endif
	if( !ppid )
	{
		signalSM = NULL;
		return;
	}

#if 0
	signalsmid = shmget( (key_t)ppid, 3, IPC_ALLOC | 0666 );
	if( signalsmid == -1 ) ErrorExit( "Cannot get Shared memory for signal.\n" );
	signalSM = shmat( signalsmid, 0, 0 );
	if( (int)signalSM == -1 ) ErrorExit( "Cannot attatch Shared Memory for signal!\n" );
	signalSM[STATUS] = IMA_KAITERU;
	signalSM[SEMAPHORE] = 1;
#endif
}

void initFiles( void ) // yes
{
	CALLS && printf("called %s:initFiles()\n", __FILE__);
	char pname[100];
	if( ppid )
		sprintf( pname, "/tmp/pre.%d", ppid );
	else
		sprintf( pname, "pre" );
	prep_g = fopen( pname, "w" );
	FILES && printf("file open w \"%s\" (%d) %s:%d\n", pname, prep_g, __FILE__, __LINE__);
	if( !prep_g ) ErrorExit( "Cannot open pre" );

	trap_g = fopen( "trace", "w" );
	FILES && printf("file open w \"trace\" (%d) %s:%d\n", trap_g, __FILE__, __LINE__);
	if( !trap_g ) ErrorExit( "cannot open trace" );
	FILES && printf("file write \"trace\" (%d) %s:%d\n", trap_g, __FILE__, __LINE__);
	fprintf( trap_g, "PID = %d\n", getpid() );
	fflush( trap_g );
}

void closeFiles( void ) // yes
{
	CALLS && printf("called %s:closeFiles()\n", __FILE__);
	FILES && printf("file close %d %s:%d\n", prep_g, __FILE__, __LINE__);
	fclose( prep_g );
	FILES && printf("file close %d %s:%d\n", trap_g, __FILE__, __LINE__);
	fclose( trap_g );
}


void WriteForFasta( FILE *fp, int locnjob, char **name, int nlen[M], char **aseq )
{
	CALLS && printf("called %s:WriteForFasta()\n", __FILE__);
    static char b[N];
    int i, j;
    int nalen[M];

    for( i=0; i<locnjob; i++ )
    {
        nalen[i] = strlen( aseq[i] );
        FILES && printf("file write %d %s:%d\n", fp, __FILE__, __LINE__);
        fprintf( fp, ">%s\n", name[i] );
        for( j=0; j<nalen[i]; j=j+C ) 
        {
            strncpy( b, aseq[i]+j, C ); b[C] = 0;
            fprintf( fp, "%s\n",b );
        }
    }
}

void readlocalhomtable2( FILE*fp, int njob, LocalHom **localhomtable, char *kozoarivec )
{
	CALLS && printf("called %s:readlocalhomtable2()\n", __FILE__);
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	LocalHom *tmpptr1, *tmpptr2;

//	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif

//		if( i < j )
		{
			if( localhomtable[i][j].nokori++ > 0 )
			{
				tmpptr1 = localhomtable[i][j].last;
//				fprintf( stderr, "reallocating, localhomtable[%d][%d].nokori = %d\n", i, j, localhomtable[i][j].nokori );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->extended = -1;
				tmpptr1->next = NULL;
				localhomtable[i][j].last = tmpptr1;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j].nokori );
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", i, j, localhomtable[i][j].nokori );
			}
	
			tmpptr1->start1 = start1;
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1;
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
		}
//		else
		{
			if( localhomtable[j][i].nokori++ > 0 )
			{
				tmpptr2 = localhomtable[j][i].last;
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->extended = -1;
				tmpptr2->next = NULL;
				localhomtable[j][i].last = tmpptr2;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", j, i, localhomtable[j][i].nokori );
			}
			else
			{
				tmpptr2 = localhomtable[j]+i;
//				fprintf( stderr, "### i,j = %d,%d, nokori=%d\n", j, i, localhomtable[j][i].nokori );
			}
	
			tmpptr2->start2 = start1;
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1;
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, st1=%d, en1=%d, opt = %f\n", i, j, tmpptr1->start1, tmpptr1->end1, opt );
		}

	}
}
void readlocalhomtable( FILE*fp, int njob, LocalHom **localhomtable, char *kozoarivec )
{
	CALLS && printf("called %s:readlocalhomtable()\n", __FILE__);
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	int **nlocalhom = NULL;
	LocalHom *tmpptr1=NULL, *tmpptr2=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( njob, njob );
	for( i=0; i<njob; i++ ) for( j=0; j<njob; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) kozoarivec[i] = kozoarivec[j] = 1;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif


//		if( i < j )
		{
			if( nlocalhom[i][j]++ > 0 )
			{
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );
		}
//		else
		{
			if( nlocalhom[j][i]++ > 0 )
			{
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->next = NULL;
			}
			else
				tmpptr2 = localhomtable[j]+i;
	
			tmpptr2->start2 = start1; // CHUUI!!!!
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1; // CHUUI!!!!
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;

//			fprintf( stderr, "j=%d, i=%d, opt = %f\n", j, i, opt );
		}

	}
	FreeIntMtx( nlocalhom );
}


void readlocalhomtable_two( FILE*fp, int norg, int nadd, LocalHom **localhomtable, LocalHom **localhomtablex, char *kozoarivec ) // for test only
{
	CALLS && printf("called %s:readlocalhomtable_two()\n", __FILE__);
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	int **nlocalhom = NULL;
	int **nlocalhomx = NULL;
	LocalHom *tmpptr1=NULL, *tmpptr2=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( norg, nadd );
	for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ ) nlocalhom[i][j] = 0;
	nlocalhomx = AllocateIntMtx( nadd, norg );
	for( i=0; i<nadd; i++ ) for( j=0; j<norg; j++ ) nlocalhomx[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) 
		{
			fprintf( stderr, "Not supported!\n" );
			exit( 1 );
		}
		j -= norg;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif


//		if( i < j )
		{
			if( nlocalhom[i][j]++ > 0 )
			{
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );
		}

		{
			if( nlocalhomx[j][i]++ > 0 )
			{
				tmpptr2->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr2 = tmpptr2->next;
				tmpptr2->next = NULL;
			}
			else
				tmpptr2 = localhomtablex[j]+i;
	
			tmpptr2->start2 = start1+1; // CHUUI!!!!
			tmpptr2->start1 = start2;
			tmpptr2->end2 = end1+1; // CHUUI!!!!
			tmpptr2->end1 = end2;
//			tmpptr2->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr2->opt = opt;
			tmpptr2->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr2->overlapaa = overlapaa;
			tmpptr2->korh = *infor;

//			fprintf( stderr, "j=%d, i=%d, opt = %f\n", j, i, opt );
		}

	}
	FreeIntMtx( nlocalhom );
	FreeIntMtx( nlocalhomx );
}

void readlocalhomtable_one( FILE*fp, int norg, int nadd, LocalHom **localhomtable, char *kozoarivec ) // for test only
{
	CALLS && printf("called %s:readlocalhomtable_one()\n", __FILE__);
	double opt;
	static char buff[B];
	char infor[100];
	int i, j, overlapaa, start1, end1, start2, end2;
	int **nlocalhom = NULL;
	LocalHom *tmpptr1=NULL; // by D.Mathog, a guess

	nlocalhom = AllocateIntMtx( norg, nadd );
	for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ ) nlocalhom[i][j] = 0;

	while ( NULL != fgets( buff, B-1, fp ) )
	{
//		fprintf( stderr, "\n" );
		sscanf( buff, "%d %d %d %lf %d %d %d %d %s",  &i, &j, &overlapaa, &opt, &start1, &end1, &start2, &end2, infor );
		if( *infor == 'k' ) 
		{
			fprintf( stderr, "Not supported!\n" );
			exit( 1 );
		}
		j -= norg;

#if 0
		if( start1 == end1 || start2 == end2 ) continue; //mondai ari
#endif


//		if( i < j )
		{
			if( nlocalhom[i][j]++ > 0 )
			{
//				fprintf( stderr, "reallocating, nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
				tmpptr1->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
				tmpptr1 = tmpptr1->next;
				tmpptr1->next = NULL;
			}
			else
			{
				tmpptr1 = localhomtable[i]+j;
//				fprintf( stderr, "nlocalhom[%d][%d] = %d\n", i, j, nlocalhom[i][j] );
			}
	
			tmpptr1->start1 = start1; // CHUUI!!!!
			tmpptr1->start2 = start2;
			tmpptr1->end1 = end1; // CHUUI!!!!
			tmpptr1->end2 = end2;
//			tmpptr1->opt = ( opt / overlapaa + 0.00 ) / 5.8  * 600;
//			tmpptr1->opt = opt;
			tmpptr1->opt = ( opt + 0.00 ) / 5.8  * 600;
			tmpptr1->overlapaa = overlapaa;
			tmpptr1->korh = *infor;
	
//			fprintf( stderr, "i=%d, j=%d, opt = %f\n", i, j, opt );
		}

	}
	FreeIntMtx( nlocalhom );
}

void outlocalhom_part( LocalHom **localhom, int norg, int nadd )
{
	CALLS && printf("called %s:outlocalhom_part()\n", __FILE__);
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<norg; i++ ) for( j=0; j<nadd; j++ )
	{
		tmpptr = localhom[i]+j;
		fprintf( stderr, "%d-%d\n", i, j+norg );
		do
		{
			fprintf( stderr, "reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt );
		}
		while( (tmpptr=tmpptr->next) );
	}
}

void outlocalhom( LocalHom **localhom, int nseq )
{
	CALLS && printf("called %s:outlocalhom()\n", __FILE__);
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j;
		fprintf( stderr, "%d-%d\n", i, j );
		do
		{
			fprintf( stderr, "reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt );
		}
		while( (tmpptr=tmpptr->next) );
	}
}

void outlocalhompt( LocalHom ***localhom, int n1, int n2 )
{
	CALLS && printf("called %s:outlocalhompt()\n", __FILE__);
	int i, j;
	LocalHom *tmpptr;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		tmpptr = localhom[i][j];
//		fprintf( stdout, "%d-%d\n", i, j );
		do
		{
			fprintf( stdout, "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f, opt=%f, wimp=%f\n", i, j, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->importance, tmpptr->opt, tmpptr->wimportance );
		}
		while( (tmpptr=tmpptr->next) );
	}
}

void FreeLocalHomTable_part( LocalHom **localhomtable, int n, int m ) 
{
	CALLS && printf("called %s:FreeLocalHomTable_part()\n", __FILE__);
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<m; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}
#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

void FreeLocalHomTable_two( LocalHom **localhomtable, int n, int m ) 
{
	CALLS && printf("called %s:FreeLocalHomTable_two()\n", __FILE__);
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<m; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}

	for( i=n; i<n+m; i++ ) 
	{
		for( j=0; j<n; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}
#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

void FreeLocalHomTable_one( LocalHom **localhomtable, int n, int m ) 
{
	CALLS && printf("called %s:FreeLocalHomTable_one()\n", __FILE__);
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<m; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}

#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

void FreeLocalHomTable( LocalHom **localhomtable, int n ) 
{
	CALLS && printf("called %s:FreeLocalHomTable()\n", __FILE__);
	int i, j;
	LocalHom *ppp, *tmpptr;
	for( i=0; i<n; i++ ) 
	{
		for( j=0; j<n; j++ )
		{
			tmpptr=localhomtable[i]+j;
			ppp = tmpptr->next;
			for( ; tmpptr; tmpptr=ppp )
			{
#if DEBUG
				fprintf( stderr, "i=%d, j=%d\n", i, j ); 
#endif
				ppp = tmpptr->next;
				if( tmpptr!=localhomtable[i]+j ) 
				{
#if DEBUG
					fprintf( stderr, "freeing %p\n", tmpptr );
#endif
					free( tmpptr );
				}
			}
		}
#if DEBUG
		fprintf( stderr, "freeing localhomtable[%d]\n", i );
#endif
		free( localhomtable[i] );
	}
#if DEBUG
	fprintf( stderr, "freeing localhomtable\n" );
#endif
	free( localhomtable );
#if DEBUG
	fprintf( stderr, "freed\n" );
#endif
}

char *progName( char *str ) // yes
{
	CALLS && printf("called %s:progName()\n", __FILE__);
    char *value; 
    if( ( value = strrchr( str, '/' ) ) != NULL )
        return( value+1 );
    else    
        return( str );
}

static void tabtospace( char *str )
{
	CALLS && printf("called %s:tabtospace()\n", __FILE__);
	char *p;
//	fprintf( stderr, "before = %s\n", str );
	while( NULL != ( p = strchr( str , '\t' ) ) )
	{
		*p = ' ';
	}
//	fprintf( stderr, "after = %s\n", str );
}

static char *extractfirstword( char *str )
{
	CALLS && printf("called %s:extractfirstword()\n", __FILE__);
	char *val = str;

	tabtospace( str );
	while( *str )
	{
		if( val == str && *str == ' ' )
		{
			val++; str++;
		}
		else if( *str != ' ' )
		{
			str++;
		}
		else if( *str == ' ' )
		{
			*str = 0;
		}
	}
	return( val );
}

void phylipout_pointer( FILE *fp, int nseq, int maxlen, char **seq, char **name, int *order, int namelen )
{
	CALLS && printf("called %s:phylipout_pointer()\n", __FILE__);
	int pos, pos2, j;
	if( namelen == -1 ) namelen = 10;
	pos = 0;

	fprintf( fp, " %d %d\n", nseq, maxlen );
	
	while( pos < maxlen )
	{
		for( j=0; j<nseq; j++ )
		{
			if( pos == 0 )
				fprintf( fp, "%-*.*s", namelen, namelen, extractfirstword( name[order[j]]+1 ) );
			else
				fprintf( fp, "%-*.*s", namelen, namelen, "" );

			pos2 = pos;
			while( pos2 < pos+41 && pos2 < maxlen )
			{
				fprintf( fp, " %.10s", seq[order[j]]+pos2 );
				pos2 += 10;
			}
			fprintf( fp, "\n" );
		}
		fprintf( fp, "\n" );
		pos += 50;
	}
}

void clustalout_pointer( FILE *fp, int nseq, int maxlen, char **seq, char **name, char *mark, char *comment, int *order, int namelen )
{
	CALLS && printf("called %s:clustalout_pointer()\n", __FILE__);
	int pos, j;
	if( namelen == -1 ) namelen = 15;
	pos = 0;
	if( comment == NULL )
		fprintf( fp, "CLUSTAL format alignment by MAFFT (v%s)\n\n", VERSION );
	else
		fprintf( fp, "CLUSTAL format alignment by MAFFT %s (v%s)\n\n", comment, VERSION );
	
	while( pos < maxlen )
	{
		fprintf( fp, "\n" );
		for( j=0; j<nseq; j++ )
		{
			fprintf( fp, "%-*.*s ", namelen, namelen, extractfirstword( name[order[j]]+1 ) );
			fprintf( fp, "%.60s\n", seq[order[j]]+pos ); // 長さが違うとだめ。
		}
		if( mark )
		{
			fprintf( fp, "%-*.*s ", namelen, namelen, "" );
			fprintf( fp, "%.60s\n", mark + pos ); // 長さが違うとだめ。
		}
		pos += 60;
	}
}


void writeData_reorder_pointer( FILE *fp, int locnjob, char **name, int *nlen, char **aseq, int *order )
{
	CALLS && printf("called %s:writeData_reorder_pointer()\n", __FILE__);
	int i, j, k;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		k = order[i];
#if DEBUG
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		nalen = strlen( aseq[k] );
		fprintf( fp, ">%s\n", name[k]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[k]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[k]+j );
#endif
		}
	}
}
void writeData_reorder( FILE *fp, int locnjob, char name[][B], int nlen[], char **aseq, int *order )
{
	CALLS && printf("called %s:writeData_reorder()\n", __FILE__);
	int i, j, k;
	int nalen;

	for( i=0; i<locnjob; i++ )
	{
		k = order[i];
#if DEBUG
		fprintf( stderr, "i = %d in writeData\n", i );
#endif
		nalen = strlen( aseq[k] );
		fprintf( fp, ">%s\n", name[k]+1 );
		for( j=0; j<nalen; j=j+C )
		{
#if 0
			strncpy( b, aseq[k]+j, C ); b[C] = 0;
			fprintf( fp, "%s\n",b );
#else
			fprintf( fp, "%.*s\n", C, aseq[k]+j );
#endif
		}
	}
}
static void showaamtxexample()
{
	CALLS && printf("called %s:showaamtxexample()\n", __FILE__);
	fprintf( stderr, "Format error in aa matrix\n" );
	fprintf( stderr, "# Example:\n" );
	fprintf( stderr, "# comment\n" );
	fprintf( stderr, "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V\n" );
	fprintf( stderr, "A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0\n" );
	fprintf( stderr, "R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3\n" );
	fprintf( stderr, "...\n" );
	fprintf( stderr, "V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4\n" );
	fprintf( stderr, "frequency 0.07 0.05 0.04 0.05 0.02 .. \n" );
	fprintf( stderr, "# Example end\n" );
	fprintf( stderr, "Only the lower half is loaded\n" );
	fprintf( stderr, "The last line (frequency) is optional.\n" );
	exit( 1 );
}

double *loadaamtx( void )
{
	CALLS && printf("called %s:loadaamtx()\n", __FILE__);
	int i, j, k, ii, jj;
	double *val;
	double **raw;
	int *map;
	char *aaorder = "ARNDCQEGHILKMFPSTWYV";
	char *inorder;
	char *line;
	char *ptr1;
	char *ptr2;
	char *mtxfname = "_aamtx";
	FILE *mf;

	raw = AllocateDoubleMtx( 21, 20 );
	val = AllocateDoubleVec( 420 );
	map = AllocateIntVec( 20 );

	if( dorp != 'p' )
	{
		fprintf( stderr, "User-defined matrix is not supported for DNA\n" );
		exit( 1 );
	}

	mf = fopen( mtxfname, "r" );
	if( mf == NULL )
	{
		fprintf( stderr, "Cannot open the _aamtx file\n" );
		exit( 1 );
	}

	inorder = calloc( 1000, sizeof( char ) );
	line = calloc( 1000, sizeof( char ) );
	

	while( !feof( mf ) )
	{
		fgets( inorder, 999, mf );
		if( inorder[0] != '#' ) break;
	}
	ptr1 = ptr2 = inorder;
	while( *ptr2 )
	{
		if( isalpha( *ptr2 ) )
		{
			*ptr1 = toupper( *ptr2 );
			ptr1++;
		}
		ptr2++;
	}
	inorder[20] = 0;

	for( i=0; i<20; i++ )
	{
		ptr2 = strchr( inorder, aaorder[i] );
		if( ptr2 == NULL )
		{
			fprintf( stderr, "%c: not found in the first 20 letters.\n", aaorder[i] );
			showaamtxexample();
		}
		else
		{
			map[i] = ptr2 - inorder;
		}
	}

	i = 0;
	while( !feof( mf ) )
	{
		fgets( line, 999, mf );
//		fprintf( stderr, "line = %s\n", line );
		if( line[0] == '#' ) continue;
		ptr1 = line;
//		fprintf( stderr, "line = %s\n", line );
		for( j=0; j<=i; j++ )
		{
			while( !isdigit( *ptr1 ) && *ptr1 != '-' && *ptr1 != '.' )
				ptr1++;

			raw[i][j] = atof( ptr1 );
//			fprintf( stderr, "raw[][]=%f, %c-%c %d-%d\n", raw[i][j], inorder[i], inorder[j], i, j );
			ptr1 = strchr( ptr1, ' ' );
			if( ptr1 == NULL && j<i) showaamtxexample();
		}
		i++;
		if( i > 19 ) break;
	}

	for( i=0; i<20; i++ ) raw[20][i] = -1.0;
	while( !feof( mf ) )
	{
		fgets( line, 999, mf );
		if( line[0] == 'f' )
		{
//			fprintf( stderr, "line = %s\n", line );
			ptr1 = line;
			for( j=0; j<20; j++ )
			{
				while( !isdigit( *ptr1 ) && *ptr1 != '-' && *ptr1 != '.' )
					ptr1++;
	
				raw[20][j] = atof( ptr1 );
//				fprintf( stderr, "raw[20][]=%f, %c %d\n", raw[20][j], inorder[i], j );
				ptr1 = strchr( ptr1, ' ' );
				if( ptr1 == NULL && j<19) showaamtxexample();
			}
			break;
		}
	}

	k = 0;
	for( i=0; i<20; i++ )
	{
		for( j=0; j<=i; j++ )
		{
			if( i != j )
			{
				ii = MAX( map[i], map[j] );
				jj = MIN( map[i], map[j] );
			}
			else ii = jj = map[i];
			val[k++] = raw[ii][jj];
//			fprintf( stderr, "%c-%c, %f\n", aaorder[i], aaorder[j], val[k-1] );
		}
	}
	for( i=0; i<20; i++ ) val[400+i] = raw[20][map[i]];

	fprintf( stderr, "inorder = %s\n", inorder );
	fclose( mf );
	free( inorder );
	free( line );
	FreeDoubleMtx( raw );
	free( map );
	return( val );
}

static void tab2space( char *s ) // nen no tame
{
	CALLS && printf("called %s:tab2space()\n", __FILE__);
	while( *s )
	{
		if( *s == '\t' ) *s = ' ';
		s++;
	}
}

static int readasubalignment( char *s, int *t, int *preservegaps )
{
	CALLS && printf("called %s:readasubalignment()\n", __FILE__);
	int v = 0;
	char status = 's';
	char *pt = s;
	*preservegaps = 0;
	tab2space( s );
	while( *pt )
	{
		if( *pt == ' ' )
		{
			status = 's';
		}
		else
		{
			if( status == 's' ) 
			{
				if( *pt == '\n' || *pt == '#' ) break;
				status = 'n';
				t[v] = atoi( pt );
				if( t[v] == 0 )
				{
					fprintf( stderr, "Format error? Sequences must be specified as 1, 2, 3...\n" );
					exit( 1 );
				}
				if( t[v] < 0 ) *preservegaps = 1;
				t[v] = abs( t[v] );
				t[v] -= 1;
				v++;
			}
		}
		pt++;
	}
	t[v] = -1;
	return( v );
}

static int countspace( char *s )
{
	CALLS && printf("called %s:countspace()\n", __FILE__);
	int v = 0;
	char status = 's';
	char *pt = s;
	tab2space( s );
	while( *pt )
	{
		if( *pt == ' ' )
		{
			status = 's';
		}
		else
		{
			if( status == 's' ) 
			{
				if( *pt == '\n' || *pt == '#' ) break;
				v++;
				status = 'n';
				if( atoi( pt ) == 0 )
				{
					fprintf( stderr, "Format error? Sequences should be specified as 1, 2, 3...\n" );
					exit( 1 );
				}
			}
		}
		pt++;
	}
	return( v );
}


void readsubalignmentstable( int nseq, int **table, int *preservegaps, int *nsubpt, int *maxmempt )
{
	CALLS && printf("called %s:readsubalignmentstable()\n", __FILE__);
	FILE *fp;
	char *line;
	int linelen = 1000000;
	int nmem;
	int lpos;
	int i, p;
	int *tab01;

	line = calloc( linelen, sizeof( char ) );
	fp = fopen( "_subalignmentstable", "r" );
	if( !fp )
	{
		fprintf( stderr, "Cannot open _subalignmentstable\n" );
		exit( 1 );
	}
	if( table == NULL )
	{
		*nsubpt = 0;
		*maxmempt = 0;
		while( 1 )
		{
			fgets( line, linelen-1, fp );
			if( feof( fp ) ) break;
			if( line[strlen(line)-1] != '\n' )
			{
				fprintf( stderr, "too long line? \n" );
				exit( 1 );
			}
			if( line[0] == '#' ) continue;
			if( atoi( line ) == 0 ) continue;
			nmem = countspace( line );
			if( nmem > *maxmempt ) *maxmempt = nmem;
			(*nsubpt)++;
		}
	}
	else
	{
		tab01 = calloc( nseq, sizeof( int ) );
		for( i=0; i<nseq; i++ ) tab01[i] = 0;
		lpos = 0;
		while( 1 )
		{
			fgets( line, linelen-1, fp );
			if( feof( fp ) ) break;
			if( line[strlen(line)-1] != '\n' )
			{
				fprintf( stderr, "too long line? \n" );
				exit( 1 );
			}
			if( line[0] == '#' ) continue;
			if( atoi( line ) == 0 ) continue;
			readasubalignment( line, table[lpos], preservegaps+lpos );
			for( i=0; (p=table[lpos][i])!=-1; i++ )
			{
				if( tab01[p] )
				{
					fprintf( stderr, "\nSequence %d appears in different groups.\n", p+1 );
					fprintf( stderr, "Hierarchical grouping is not supported.\n\n" );
					exit( 1 );
				}
				tab01[p] = 1;
				if( p > nseq-1 )
				{
					fprintf( stderr, "Sequence %d does not exist in the input sequence file.\n", p+1 );
					exit( 1 );
				}
			}
			lpos++;
		}
		free( tab01 );
	}
	fclose( fp );
	free( line );
}


void readmccaskill( FILE *fp, RNApair **pairprob, int length )
{
	CALLS && printf("called %s:readmccaskill()\n", __FILE__);
	char gett[1000];
	int *pairnum;
	int i;
	int left, right;
	float prob;
	int c;

	pairnum = (int *)calloc( length, sizeof( int ) );
	for( i=0; i<length; i++ ) pairnum[i] = 0;

	c = getc( fp );
	{
		if( c != '>' )
		{
			fprintf( stderr, "format error in hat4 - 1\n" );
			exit( 1 );
		}
	}
	fgets( gett, 999, fp );
	while( 1 )
	{
		if( feof( fp ) ) break;
		c = getc( fp );
		ungetc( c, fp );
		if( c == '>' || c == EOF )
		{
			break;
		}
		fgets( gett, 999, fp );
//		fprintf( stderr, "gett = %s\n", gett );
		sscanf( gett, "%d %d %f", &left, &right, &prob );

		if( left >= length || right >= length )
		{
			fprintf( stderr, "format error in hat4 - 2\n" );
			exit( 1 );
		}

		if( prob < 0.01 ) continue; // 080607, mafft ni dake eikyou

		if( left != right && prob > 0.0 )
		{
			pairprob[left] = (RNApair *)realloc( pairprob[left], (pairnum[left]+2) * sizeof( RNApair ) );
			pairprob[left][pairnum[left]].bestscore = prob;
			pairprob[left][pairnum[left]].bestpos = right;
			pairnum[left]++;
			pairprob[left][pairnum[left]].bestscore = -1.0;
			pairprob[left][pairnum[left]].bestpos = -1;
//			fprintf( stderr, "%d-%d, %f\n", left, right, prob );

			pairprob[right] = (RNApair *)realloc( pairprob[right], (pairnum[right]+2) * sizeof( RNApair ) );
			pairprob[right][pairnum[right]].bestscore = prob;
			pairprob[right][pairnum[right]].bestpos = left;
			pairnum[right]++;
			pairprob[right][pairnum[right]].bestscore = -1.0;
			pairprob[right][pairnum[right]].bestpos = -1;
//			fprintf( stderr, "%d-%d, %f\n", right, left, prob );
		}
	}
	free( pairnum );
}

void readpairfoldalign( FILE *fp, char *s1, char *s2, char *aln1, char *aln2, int q1, int q2, int *of1, int *of2, int sumlen )
{
	CALLS && printf("called %s:readpairfoldalign()\n", __FILE__);
	char gett[1000];
	int *maptoseq1;
	int *maptoseq2;
	char dumc;
	int dumi;
	char sinseq[100], sinaln[100];
	int posinseq, posinaln;
	int alnlen;
	int i;
	int pos1, pos2;
	char *pa1, *pa2;
	char qstr[1000];

	*of1 = -1;
	*of2 = -1;

	maptoseq1 = AllocateIntVec( sumlen+1 );
	maptoseq2 = AllocateIntVec( sumlen+1 );

	posinaln = 0; // foldalign ga alingment wo kaesanaitok no tame.

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; ALIGNING", 10 ) ) break;
	}
	sprintf( qstr, "; ALIGNING            %d against %d\n", q1+1, q2+1 );
	if( strcmp( gett, qstr ) )
	{
		fprintf( stderr, "Error in FOLDALIGN\n" );
		fprintf( stderr, "qstr = %s, but gett = %s\n", qstr, gett );
		exit( 1 );
	}

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; --------", 10 ) ) break;
	}


	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; ********", 10 ) ) break;
//		fprintf( stderr, "gett = %s\n", gett );
		sscanf( gett, "%c %c %s %s %d %d", &dumc, &dumc, sinseq, sinaln, &dumi, &dumi );
		posinaln = atoi( sinaln );
		posinseq = atoi( sinseq );
//		fprintf( stderr, "posinseq = %d\n", posinseq );
//		fprintf( stderr, "posinaln = %d\n", posinaln );
		maptoseq1[posinaln-1] = posinseq-1;
	}
	alnlen = posinaln;

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; --------", 10 ) ) break;
	}

	while( !feof( fp ) )
	{
		fgets( gett, 999, fp );
		if( !strncmp( gett, "; ********", 10 ) ) break;
//		fprintf( stderr, "gett = %s\n", gett );
		sscanf( gett, "%c %c %s %s %d %d", &dumc, &dumc, sinseq, sinaln, &dumi, &dumi );
		posinaln = atof( sinaln );
		posinseq = atof( sinseq );
//		fprintf( stderr, "posinseq = %d\n", posinseq );
//		fprintf( stderr, "posinaln = %d\n", posinaln );
		maptoseq2[posinaln-1] = posinseq-1;
	}
	if( alnlen != posinaln )
	{
		fprintf( stderr, "Error in foldalign?\n" );
		exit( 1 );
	}

	pa1 = aln1;
	pa2 = aln2;
	for( i=0; i<alnlen; i++ )
	{
		pos1 = maptoseq1[i];
		pos2 = maptoseq2[i];

		if( pos1 > -1 )
			*pa1++ = s1[pos1];
		else
			*pa1++ = '-';

		if( pos2 > -1 )
			*pa2++ = s2[pos2];
		else
			*pa2++ = '-';
	}
	*pa1 = 0;
	*pa2 = 0;

	*of1 = 0;
	for( i=0; i<alnlen; i++ )
	{
		*of1 = maptoseq1[i];
		if( *of1 > -1 ) break;
	}
	*of2 = 0;
	for( i=0; i<alnlen; i++ )
	{
		*of2 = maptoseq2[i];
		if( *of2 > -1 ) break;
	}

//	fprintf( stderr, "*of1=%d, aln1 = :%s:\n", *of1, aln1 );
//	fprintf( stderr, "*of2=%d, aln2 = :%s:\n", *of2, aln2 );

	free( maptoseq1 );
	free( maptoseq2 );
}

int myatoi( char *in ) // yes
{
	CALLS && printf("called %s:myatoi()\n", __FILE__);
	if( in == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Error in myatoi()\n" );
		exit( 1 );
	}
	return( atoi( in ) );
}

float myatof( char *in ) // yes
{
	CALLS && printf("called %s:myatof()\n", __FILE__);
	if( in == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Error in myatof()\n" );
		exit( 1 );
	}
	return( atof( in ) );
}

void reporterr( const char *str, ... ) // yes
{
	CALLS && printf("called %s:reporterr()\n", __FILE__);
//	static int loglen = 0;
	va_list args;

	if( gmsg )
	{
# if 1  // ato de sakujo
		static FILE *errtmpfp = NULL;
		if( errtmpfp == NULL ) {
			errtmpfp = fopen( "maffterr", "w" );
			FILES && printf("file open w \"maffterr\" (%d) %s:%d\n", errtmpfp, __FILE__, __LINE__);
		} else {
			errtmpfp = fopen( "maffterr", "a" );
			FILES && printf("file open a \"maffterr\" (%d) %s:%d\n", errtmpfp, __FILE__, __LINE__);
		}
		va_start( args, str );
		FILES && printf("file write %d %s:%d\n", errtmpfp, __FILE__, __LINE__);
		vfprintf( errtmpfp, str, args );
		va_end( args );
		FILES && printf("file close %d %s:%d\n", errtmpfp, __FILE__, __LINE__);
		fclose( errtmpfp );
#endif

#if 0
		char *tmpptr;
		tmpptr = (char *)realloc( *gmsg, (loglen+10000) * sizeof( char ) );
		if( tmpptr == NULL )
		{
			fprintf( stderr, "Cannot relloc *gmsg\n" );
			exit( 1 );
		}
		*gmsg = tmpptr;
		va_start( args, str );
		loglen += vsprintf( *gmsg + loglen, str, args );
		va_end( args );


		va_start( args, str );
		loglen += vsprintf( *gmsg + loglen, str, args );
		va_end( args );
		*(*gmsg + loglen) = 0;
		if( loglen > gmsglen - 100 ) loglen = 0; // tekitou
#endif

	}
	else
	{
		va_start( args, str );
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		vfprintf( stderr, str, args );
		va_end( args );
//		fflush( stderr ); // iru?
	}
	return;
}
