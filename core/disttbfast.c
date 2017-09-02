#include "mltaln.h"

#ifdef PCALLS
#define CALLS 1
#else
#define CALLS 0
#endif
#ifdef PBRANCHES
#define BRANCHES 1
#else
#define BRANCHES 0
#endif
#ifdef PFILES
#define FILES 1
#else
#define FILES 0
#endif

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SKIP 1

#define END_OF_VEC -1

static int nadd;
static int treein;
static int topin;
static int treeout;
static int noalign;
static int distout;
static float lenfaca, lenfacb, lenfacc, lenfacd;
static int tuplesize;
static int subalignment;
static int subalignmentoffset;
static int nguidetree;
static int sparsepickup;
#if 0
#define PLENFACA 0.0123
#define PLENFACB 10252
#define PLENFACC 10822
#define PLENFACD 0.5
#define DLENFACA 0.01
#define DLENFACB 2445
#define DLENFACC 2412
#define DLENFACD 0.1
#else
#define PLENFACA 0.01
#define PLENFACB 10000
#define PLENFACC 10000
#define PLENFACD 0.1
#define D6LENFACA 0.01
#define D6LENFACB 2500
#define D6LENFACC 2500
#define D6LENFACD 0.1
#define D10LENFACA 0.01
#define D10LENFACB 1000000
#define D10LENFACC 1000000
#define D10LENFACD 0.0
#endif

typedef struct _jobtable
{
    int i;  
    int j;  
} Jobtable;

typedef struct _msadistmtxthread_arg
{
	int njob;
	int thread_no;
	float *selfscore;
	float **iscore;
	char **seq;
	int **skiptable;
	Jobtable *jobpospt;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msadistmtxthread_arg_t;

#ifdef enablemultithread
// ue futatsu ha singlethread demo tsukau
typedef struct _treebasethread_arg
{
	int thread_no;
	int njob;
	int *nrunpt;
	int *nlen;
	int *jobpospt;
	int ***topol;
	Treedep *dep;
	char **aseq;
	double *effarr;
	int *alloclenpt;
	int *fftlog;
	char *mergeoralign;
	float **newdistmtx;
	float *selfscore;
	pthread_mutex_t *mutex;
	pthread_cond_t *treecond;
} treebasethread_arg_t;

typedef struct _distancematrixthread_arg
{
	int thread_no;
	int njob;
	int *jobpospt;
	int **pointt;
	float **mtx;
	pthread_mutex_t *mutex;
} distancematrixthread_arg_t;
#endif


void arguments( int argc, char *argv[] )
{
	CALLS && printf("called %s:arguments()\n", __FILE__);
    int c;

	nthread = 1;
	outnumber = 0;
	topin = 0;
	treein = 0;
	treeout = 0;
	distout = 0;
	noalign = 0;
	nevermemsave = 0;
	inputfile = NULL;
	nadd = 0;
	addprofile = 1;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0;
	force_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'X';
	sueff_global = 0.1;
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
	dorp = NOTSPECIFIED;
	ppenalty_dist = NOTSPECIFIED;
	ppenalty = -1530;
	ppenalty_ex = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	poffset = -123;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	TMorJTT = JTT;
	scoreout = 0;
	spscoreout = 0;
	tuplesize = 6;
	subalignment = 0;
	subalignmentoffset = 0;
	legacygapcost = 0;
	specificityconsideration = 0.0;
	nguidetree = 1;
	sparsepickup = 0;

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
        while ( (c = *++argv[0]) )
		{
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
					reporterr(       "inputfile = %s\n", inputfile );
					--argc;
					goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
					reporterr(       "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'V':
					ppenalty_dist = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
					reporterr(       "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					reporterr(       "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
					reporterr(       "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					reporterr(       "blosum %d / kimura 200 \n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
					reporterr(       "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
					reporterr(       "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
					reporterr(       "nthread = %d\n", nthread );
					--argc; 
					goto nextoption;
				case 's':
					specificityconsideration = (double)myatof( *++argv );
//					reporterr(       "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
				case 'y':
					distout = 1;
					break;
				case 't':
					treeout = 1;
					break;
				case 'T':
					noalign = 1;
					break;
#if 0
				case 'r':
					fmodel = -1;
					break;
#endif
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'L':
					legacygapcost = 1;
					break;
				case 'e':
					fftscore = 0;
					break;
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'R':
					fftRepeatStop = 1;
					break;
#endif
				case 'n' :
					outnumber = 1;
					break;
#if 0
				case 's':
					treemethod = 's';
					break;
				case 'q':
					treemethod = 'q'; // minimum
					break;
#endif
				case 'q':
					sparsepickup = myatoi( *++argv );
					reporterr(       "sparsepickup = %d\n", sparsepickup );
					--argc; 
					goto nextoption;
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
					FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
				case 'E':
					nguidetree = myatoi( *++argv );
					reporterr(       "nguidetree = %d\n", nguidetree );
					--argc; 
					goto nextoption;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'H':
					alg = 'H';
					break;
#endif
				case 'R':
					alg = 'R';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'M':
					alg = 'M';
					break;
#if 0
				case 'S' :
					scoreout = 1; // for checking parallel calculation
					break;
#else
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
#endif
				case 'B': // hitsuyou! memopt -M -B no tame
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					use_fft = 1;
					force_fft = 1;
					break;
#if 0
				case 'V':
					topin = 1;
					break;
#endif
				case 'U':
					treein = 1;
					break;
				case 'u':
					weight = 0;
					tbrweight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					disp = 1;
					break;
#if 1
				case 'O':
					outgap = 0;
					break;
#else
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
				case 'J':
					tbutree = 0;
					break;
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'W':
					tuplesize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'Z':
					checkC = 1;
					break;
                default:
                    reporterr(       "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}
    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        reporterr(       "options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		reporterr(       "conflicting options : o, m or u\n" );
		exit( 1 );
	}
}


static void pickup( int n, int *seqlen, int ***topol, char **name, char **seq )
{
	int i, j, k, m;
	int **longestseq;
	int **longestlen;
	int *select;
	char **nameout, **seqout;
	int *nlenout;
	char **namenotused, **seqnotused;
	int *nlennotused;
	FILE *notusedfp;

	longestseq = AllocateIntMtx( n-1, 2 );
	longestlen = AllocateIntMtx( n-1, 2 );
	select = AllocateIntVec( n );
	for( i=0; i<n; i++ ) select[i] = 0;
	nameout = AllocateCharMtx( n, 0 );
	seqout = AllocateCharMtx( n, 0 );
	nlenout = AllocateIntVec( n );
	namenotused = AllocateCharMtx( n, 0 );
	seqnotused = AllocateCharMtx( n, 0 );
	nlennotused = AllocateIntVec( n );

	for( i=0; i<n-1; i++ )
	{
//		reporterr( "STEP %d\n", i );
		longestlen[i][0] = -1;
		longestseq[i][0] = -1;
		for( j=0; (m=topol[i][0][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][0] )
			{
				longestlen[i][0] = seqlen[m];
				longestseq[i][0] = m;
			}
//			reporterr( "%d ", topol[i][0][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][0], longestseq[i][0] );


		longestlen[i][1] = -1;
		longestseq[i][1] = -1;
		for( j=0; (m=topol[i][1][j])!=-1; j++ ) // sukoshi muda
		{
			if( seqlen[m] > longestlen[i][1] )
			{
				longestlen[i][1] = seqlen[m];
				longestseq[i][1] = m;
			}
//			reporterr( "%d ", topol[i][1][j] );
		}
//		reporterr( "longest = %d (%d)\n", longestlen[i][1], longestseq[i][1] );
	}

	m = 1;
	for( i=n-2; i>-1; i-- )
	{
//		reporterr( "longest[%d][0] = %d (%d)\n", i, longestlen[i][0], longestseq[i][0] );
//		reporterr( "longest[%d][1] = %d (%d)\n", i, longestlen[i][1], longestseq[i][1] );
		select[longestseq[i][0]] = 1;
		select[longestseq[i][1]] = 1;
		m += 1;
		if( m >= sparsepickup ) break;
	}
	for( i=0, k=0, j=0; i<n; i++ ) 
	{
		if( select[i] )
		{
			nameout[k] = name[i];
			seqout[k] = seq[i];
			nlenout[k] = strlen( seqout[k] );
			k++;
		}
		else
		{
			namenotused[j] = name[i];
			seqnotused[j] = seq[i];
			nlennotused[j] = strlen( seqnotused[j] );
			j++;
		}
	}
	writeData_pointer( stdout, m, nameout, nlenout, seqout );

	notusedfp = fopen( "notused", "w" );
	FILES && printf("file open w \"notused\" (%d) %s:%d\n", notusedfp, __FILE__, __LINE__);
	writeData_pointer( notusedfp, n-m, namenotused, nlennotused, seqnotused );
	FILES && printf("file close \"notused\" %s:%d\n", __FILE__, __LINE__);
	fclose( notusedfp );


	free( nameout );
	free( nlenout );
	free( seqout );
	free( namenotused );
	free( nlennotused );
	free( seqnotused );
	FreeIntMtx( longestseq );
	FreeIntMtx( longestlen );
	free( select );
}


static int maxl;
static int tsize;
static int nunknown = 0;

void seq_grp_nuc( int *grp, char *seq )
{
	CALLS && printf("called %s:seq_grp_nuc()\n", __FILE__);
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 4 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < tuplesize )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void seq_grp( int *grp, char *seq )
{
	CALLS && printf("called %s:seq_grp()\n", __FILE__);
	int tmp;
	int *grpbk = grp;
	while( *seq )
	{
		tmp = amino_grp[(int)*seq++];
		if( tmp < 6 )
			*grp++ = tmp;
		else
			nunknown++;
	}
	*grp = END_OF_VEC;
	if( grp - grpbk < 6 )
	{
//		reporterr(       "\n\nWARNING: Too short.\nPlease also consider use mafft-ginsi, mafft-linsi or mafft-ginsi.\n\n\n" );
//		exit( 1 );
		*grpbk = -1;
	}
}

void makecompositiontable_p( int *table, int *pointt )
{
	CALLS && printf("called %s:makecompositiontable_p()\n", __FILE__);
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
		table[point]++;
}

int commonsextet_p( int *table, int *pointt )
{
	CALLS && printf("called %s:commonsextet_p()\n", __FILE__);
	int value = 0;
	int tmp;
	int point;
	static TLS int *memo = NULL;
	static TLS int *ct = NULL;
	static TLS int *cp;

	if( table == NULL )
	{
		if( memo ) free( memo );
		if( ct ) free( ct );
		memo = NULL;
		ct = NULL;
		return( 0 );
	}

	if( *pointt == -1 )
		return( 0 );

	if( !memo )
	{
		memo = (int *)calloc( tsize, sizeof( int ) );
		if( !memo ) ErrorExit( "Cannot allocate memo\n" );
		ct = (int *)calloc( MIN( maxl, tsize )+1, sizeof( int ) ); // chuui!!
		if( !ct ) ErrorExit( "Cannot allocate ct\n" );
	}

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
	}
	*cp = END_OF_VEC;
	
	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	return( value );
}

void makepointtable_nuc_dectet( int *pointt, int *n )
{
	CALLS && printf("called %s:makepointtable_nuc_dectet()\n", __FILE__);
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *262144;
	point += *n++ * 65536;
	point += *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ *262144;
		point *= 4;
		point += *n++;
		*pointt++ = point;

	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc_octet( int *pointt, int *n )
{
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ * 16384;
	point += *n++ *  4096;
	point += *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 16384;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable_nuc( int *pointt, int *n )
{
	CALLS && printf("called %s:makepointtable_nuc()\n", __FILE__);
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  1024;
	point += *n++ *   256;
	point += *n++ *    64;
	point += *n++ *    16;
	point += *n++ *     4;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 1024;
		point *= 4;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}

void makepointtable( int *pointt, int *n )
{
	CALLS && printf("called %s:makepointtable()\n", __FILE__);
	int point;
	register int *p;

	if( *n == -1 )
	{
		*pointt = -1;
		return;
	}

	p = n;
	point  = *n++ *  7776;
	point += *n++ *  1296;
	point += *n++ *   216;
	point += *n++ *    36;
	point += *n++ *     6;
	point += *n++;
	*pointt++ = point;

	while( *n != END_OF_VEC )
	{
		point -= *p++ * 7776;
		point *= 6;
		point += *n++;
		*pointt++ = point;
	}
	*pointt = END_OF_VEC;
}


static void *msadistmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	CALLS && printf("called %s:msadistmtxthread()\n", __FILE__);
	msadistmtxthread_arg_t *targ = (msadistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	float *selfscore = targ->selfscore;
	float **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;


	float ssi, ssj, bunbo, iscoretmp;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthread ) pthread_mutex_lock( targ->mutex );
#endif
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
#ifdef enablemultithread
				if( nthread ) pthread_mutex_unlock( targ->mutex );
#endif
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
#ifdef enablemultithread
		if( nthread ) pthread_mutex_unlock( targ->mutex );
#endif


		if( nthread )
		{
			if( j==i+1 && i % 10 == 0 ) {
				FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
				fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
			}
		}
		else
		{
			if( j==i+1 && i % 10 == 0 ) {
				FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
				fprintf( stderr, "\r% 5d / %d", i, njob );
			}
		}
		ssi = selfscore[i];
		ssj = selfscore[j];
		bunbo = MIN( ssi, ssj );
//fprintf( stderr, "bunbo = %f\n", bunbo );
//fprintf( stderr, "naivepairscorefast() = %f\n", naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) );
		if( bunbo == 0.0 )
			iscoretmp = 2.0; // 2013/Oct/17
		else
		{
			iscoretmp = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast 
			if( iscoretmp > 10 ) iscoretmp = 10.0;  // 2015/Mar/17

		}
		iscore[i][j-i] = iscoretmp;
	}
}

#ifdef enablemultithread
static void *distancematrixthread( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int *jobpospt = targ->jobpospt;
	int **pointt = targ->pointt;
	float **mtx = targ->mtx;

	int *table1;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = *jobpospt;
		if( i == njob )
		{
			pthread_mutex_unlock( targ->mutex );
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		*jobpospt = i+1;
		pthread_mutex_unlock( targ->mutex );

		table1 = (int *)calloc( tsize, sizeof( int ) );
		if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
		if( i % 10 == 0 )
		{
			reporterr(       "\r% 5d / %d (thread %4d)", i+1, njob, thread_no );
		}
		makecompositiontable_p( table1, pointt[i] );

		for( j=i; j<njob; j++ ) 
		{
			mtx[i][j-i] = (float)commonsextet_p( table1, pointt[j] );
		} 
		free( table1 );
	}
}


static void *treebasethread( void *arg )
{
	treebasethread_arg_t *targ = (treebasethread_arg_t *)arg;
	int thread_no = targ->thread_no;
	int *nrunpt = targ->nrunpt;
	int njob = targ->njob;
	int *nlen = targ->nlen;
	int *jobpospt = targ->jobpospt;
	int ***topol = targ->topol;
	Treedep *dep = targ->dep;
	char **aseq = targ->aseq;
	double *effarr = targ->effarr;
	int *alloclen = targ->alloclenpt;
	int *fftlog = targ->fftlog;
	char *mergeoralign = targ->mergeoralign;
	float **newdistmtx = targ->newdistmtx;
	float *selfscore = targ->selfscore;

	char **mseq1, **mseq2;
	char **localcopy;
	int i, m, j, l;
	int immin, immax;
	int len1, len2;
	int clus1, clus2;
	float pscore, tscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	float dumfl = 0.0;
	int ffttry;
	int m1, m2;
	double **dynamicmtx;
	float ssi, ssm, bunbo;
	int tm, ti;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif

	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	localcopy = calloc( njob, sizeof( char * ) );
	for( i=0; i<njob; i++ ) localcopy[i] = NULL;
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );


#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		l = *jobpospt;
		if( l == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
			Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
			free( mseq1 );
			free( mseq2 );
			free( localcopy );
			free( effarr1 );
			free( effarr2 );
			free( indication1 );
			free( indication2 );
			FreeDoubleMtx( dynamicmtx );
			return( NULL );
		}
		*jobpospt = l+1;

		if( dep[l].child0 != -1 )
		{
			while( dep[dep[l].child0].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
		if( dep[l].child1 != -1 ) 
		{
			while( dep[dep[l].child1].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
		while( *nrunpt >= nthread )
			pthread_cond_wait( targ->treecond, targ->mutex );
		(*nrunpt)++;


		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
			dep[l].done = 1;
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
			free( topol[l][0] ); topol[l][0] = NULL;
			free( topol[l][1] ); topol[l][1] = NULL;
			free( topol[l] ); topol[l] = NULL;
			pthread_mutex_unlock( targ->mutex );
			continue;
		}



		m1 = topol[l][0][0];
		m2 = topol[l][1][0];

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip - 0.5 );
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );

		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen <= len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}

		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
		}

		if( alg == 'M' ) // hoka no thread ga M ni shitakamo shirenainode
		{
//			reporterr(       "Freeing commonIP (thread %d)\n", thread_no );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

		pthread_mutex_unlock( targ->mutex );

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( topol[l][0], localcopy, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], localcopy, mseq2, effarr2, effarr, indication2, 0.0 );
#else
		clus1 = fastconjuction_noweight( topol[l][0], localcopy, mseq1, effarr1,  indication1 );
		clus2 = fastconjuction_noweight( topol[l][1], localcopy, mseq2, effarr2,  indication2 );
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		reporterr(       "\rSTEP % 5d / %d (thread %4d)", l+1, njob-1, thread_no );

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/



//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (float)len1/fftlog[m1], clus1, (float)len2/fftlog[m2], clus2 );

		if( force_fft || ( use_fft && ffttry ) )
		{
			reporterr(       "f" );
			if( alg == 'M' )
			{
				reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
			}
		}
		else
		{
			reporterr(       "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					reporterr(       "m" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

		if( disp ) display( localcopy, njob );

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
				ti = topol[l][0][i];
				ssi = selfscore[topol[l][0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[topol[l][1][m]];
					tm = topol[l][1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0.0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}



		pthread_mutex_lock( targ->mutex );
		dep[l].done = 1;
		(*nrunpt)--;
		pthread_cond_broadcast( targ->treecond );

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		pthread_mutex_unlock( targ->mutex );



		for( i=0; (j=topol[l][0][i])!=-1; i++ )
		{
			if(localcopy[j] ) free( localcopy[j] );
			localcopy[j] = NULL;
		}
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
		{
			if( localcopy[j] ) free( localcopy[j] );
			localcopy[j] = NULL;
		}


		if( topol[l][0] ) free( topol[l][0] );
		topol[l][0] = NULL;
		if( topol[l][1] ) free( topol[l][1] );
		topol[l][1] = NULL;
		if( topol[l] ) free( topol[l] );
		topol[l] = NULL;


//		reporterr(       "\n" );
	}
#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
}
#endif

static int treebase( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, double *effarr, float **newdistmtx, float *selfscore, int *alloclen, int (*callback)(int, int, char*) )
{
	CALLS && printf("called %s:treebase()\n", __FILE__);
	int l, len1, len2, i, m, immin, immax;
	int len1nocommongap, len2nocommongap;
	int clus1, clus2;
	float pscore, tscore;
	char *indication1 = NULL, *indication2 = NULL;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	int *fftlog = NULL; // fixed at 2006/07/26
	float dumfl = 0.0;
	int ffttry;
	int m1, m2;
	int *gaplen = NULL;
	int *gapmap = NULL;
	int *alreadyaligned = NULL;
	double **dynamicmtx = NULL;
	float ssi, ssm, bunbo;
	int tm, ti;
#if SKIP
	int **skiptable1 = NULL, **skiptable2 = NULL;
#endif
#if 0
	int i, j;
#endif

//	reporterr( "treebase newdistmtx=%p\n", newdistmtx );

	if( effarr1 == NULL ) 
	{
		effarr1 = AllocateDoubleVec( njob );
		effarr2 = AllocateDoubleVec( njob );
		indication1 = AllocateCharVec( 150 );
		indication2 = AllocateCharVec( 150 );
		fftlog = AllocateIntVec( njob );
		gaplen = AllocateIntVec( *alloclen+10 );
		gapmap = AllocateIntVec( *alloclen+10 );
		alreadyaligned = AllocateIntVec( njob );
		dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
	}
	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	if( callback && callback( 0, 50, "Progressive alignment" ) ) goto chudan_tbfast;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0
	reporterr(       "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		reporterr(       "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

//	for( i=0; i<njob; i++ ) strcpy( aseq[i], seq[i] );


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{

//		reporterr(       "\ndistfromtip = %f\n", dep[l].distfromtip );
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );
		if( mergeoralign[l] == 'n' )
		{
//			reporterr(       "SKIP!\n" );
			free( topol[l][0] ); topol[l][0] = NULL;
			free( topol[l][1] ); topol[l][1] = NULL;
			free( topol[l] ); topol[l] = NULL;
			continue;
		}

		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
		len1 = strlen( aseq[m1] );
		len2 = strlen( aseq[m2] );
		if( *alloclen < len1 + len2 )
		{
			reporterr(       "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  ); 
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				reporterr(       "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				reporterr(       "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			reporterr(       "done. *alloclen = %d\n", *alloclen );
		}

#if 1 // CHUUI@@@@
		clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, 0.0 );
		clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, 0.0 );
#else
		clus1 = fastconjuction_noweight( topol[l][0], aseq, mseq1, effarr1,  indication1 );
		clus2 = fastconjuction_noweight( topol[l][1], aseq, mseq2, effarr2,  indication2 );
#endif
		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' )
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";

		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap );
			commongappick( clus2, mseq2 );
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap );
			commongappick( clus1, mseq1 );
			len1nocommongap = strlen( mseq1[0] );
		}

#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after conjuction) !\n" ); 
            exit( 1 );
        }
    }
#endif


#if 0
    for( i=0; i<clus1; i++ ) 
    {
        if( strlen( mseq1[i] ) != len1 ) 
        {
            reporterr(       "i = %d / %d\n", i, clus1 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
    for( j=0; j<clus2; j++ )
    {
        if( strlen( mseq2[j] ) != len2 ) 
        {
            reporterr(       "j = %d / %d\n", j, clus2 ); 
            reporterr(       "hairetsu ga kowareta (in treebase, after free topol) !\n" ); 
            exit( 1 );
        }
    }
#endif


//		fprintf( trap_g, "\nSTEP-%d\n", l );
//		fprintf( trap_g, "group1 = %s\n", indication1 );
//		fprintf( trap_g, "group2 = %s\n", indication2 );

//		reporterr(       "\rSTEP % 5d / %d %d-%d", l+1, njob-1, clus1, clus2 );
		reporterr(       "\rSTEP % 5d / %d ", l+1, njob-1 );
		if( callback && callback( 0, 50+50*l/(njob-1), "Progressive alignment" ) ) goto chudan_tbfast;

#if 0
		reporterr(       "STEP %d /%d\n", l+1, njob-1 );
		reporterr(       "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
		reporterr(       "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) reporterr(       "..." );
		reporterr(       "\n" );
#endif

/*
		reporterr(       "before align all\n" );
		display( aseq, njob );
		reporterr(       "\n" );
		reporterr(       "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		reporterr(       "\n" );
		reporterr(       "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		reporterr(       "\n" );
*/


		if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000  ) ) )
		{
			reporterr(       "\nlen1=%d, len2=%d, Switching to the memsave mode\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000);
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000); // v6.708
//		reporterr(       "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (float)len1/fftlog[m1], clus1, (float)len2/fftlog[m2], clus2 );

		if( force_fft || ( use_fft && ffttry ) )
		{
			reporterr(       "f" );
			if( alg == 'M' )
			{
				reporterr(       "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
			{
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
//				reporterr(       "######### mseq1[0] = %s\n", mseq1[0] );
			}
		}
		else
		{
			reporterr(       "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					reporterr(       "m" );
//					reporterr(       "%d-%d", clus1, clus2 );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'A' ):
					if( clus1 == 1 && clus2 == 1 )
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = G__align11( dynamicmtx, mseq1, mseq2, *alloclen, outgap, outgap );
					}
					else
					{
//						reporterr(       "%d-%d", clus1, clus2 );
						pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumfl, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					}
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}
#if SCOREOUT
		reporterr(       "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );
//		reporterr(       "\n" );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
			adjustgapmap( strlen( mseq2[0] )-len2nocommongap+len2, gapmap, mseq2[0] );
			restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
			findnewgaps( clus2, 0, mseq2, gaplen );
			insertnewgaps( njob, alreadyaligned, aseq, topol[l][1], topol[l][0], gaplen, gapmap, *alloclen, alg, '-' );
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; (m=topol[l][0][i])>-1; i++ ) alreadyaligned[m] = 1;
		}
		if( mergeoralign[l] == '2' )
		{
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP0 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP0 mseq2[%d] = \n%s\n", i, mseq2[i] );
			adjustgapmap( strlen( mseq1[0] )-len1nocommongap+len1, gapmap, mseq1[0] );
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP1 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP1 mseq2[%d] = \n%s\n", i, mseq2[i] );
			restorecommongaps( njob, aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP2 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP2 mseq2[%d] = \n%s\n", i, mseq2[i] );
			findnewgaps( clus1, 0, mseq1, gaplen );
			insertnewgaps( njob, alreadyaligned, aseq, topol[l][0], topol[l][1], gaplen, gapmap, *alloclen, alg, '-' );
//			for( i=0; i<clus1; i++ ) reporterr(       ">STEP3 mseq1[%d] = \n%s\n", i, mseq1[i] );
//			for( i=0; i<clus2; i++ ) reporterr(       ">STEP3 mseq2[%d] = \n%s\n", i, mseq2[i] );
#if 0
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; i<clus1; i++ ) 
			{
				reporterr( "mseq1[%d] bef change = %s\n", i, mseq1[i] );
				eq2dash( mseq1[i] );
				reporterr( "mseq1[%d] aft change = %s\n", i, mseq1[i] );
			}
			for( i=0; i<clus2; i++ ) 
			{
				reporterr( "mseq2[%d] bef change = %s\n", i, mseq2[i] );
				eq2dash( mseq2[i] );
				reporterr( "mseq2[%d] aft change = %s\n", i, mseq2[i] );
			}
			for( i=0; i<clus1; i++ ) eq2dash( mseq1[i] );
			for( i=0; i<clus2; i++ ) eq2dash( mseq2[i] );
#else
			eq2dashmatometehayaku( mseq1, clus1 );
			eq2dashmatometehayaku( mseq2, clus2 );
#endif
			for( i=0; (m=topol[l][1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

		if( newdistmtx ) // tsukawanai
		{
#if 0
			reporterr( "group1 = " );
			for( i=0; i<clus1; i++ ) reporterr( "%d ", topol[l][0][i] );
			reporterr( "\n" );
			reporterr( "group2 = " );
			for( m=0; m<clus2; m++ ) reporterr( "%d ", topol[l][1][m] );
			reporterr( "\n" );
#endif
#if SKIP
			skiptable1 = AllocateIntMtx( clus1, 0 );
			skiptable2 = AllocateIntMtx( clus2, 0 );
			makeskiptable( clus1, skiptable1, mseq1 ); // allocate suru.
			makeskiptable( clus2, skiptable2, mseq2 ); // allocate suru.
#endif
			for( i=0; i<clus1; i++ ) 
			{
#if SKIP
//				makeskiptable( 1, skiptable1, mseq1+i ); // allocate suru.
#endif
				ti = topol[l][0][i];
				ssi = selfscore[topol[l][0][i]];
				for( m=0; m<clus2; m++ )
				{
					ssm = selfscore[topol[l][1][m]];
					tm = topol[l][1][m];
					if( ti<tm )
					{
						immin = ti;
						immax = tm;
					}
					else
					{
						immin = tm;
						immax = ti;
					}
					bunbo = MIN( ssi, ssm );
					if( bunbo == 0.0 )
						newdistmtx[immin][immax-immin] = 2.0; // 2013/Oct/17
					else
#if SKIP
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscorefast( mseq1[i], mseq2[m], skiptable1[i], skiptable2[m], penalty_dist ) / bunbo ) * 2.0;
#else
						newdistmtx[immin][immax-immin] = ( 1.0 - naivepairscore11( mseq1[i], mseq2[m], penalty_dist ) / bunbo ) * 2.0;
#endif
				}
			}
#if SKIP
			FreeIntMtx( skiptable1 ); skiptable1 = NULL;
			FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif
		}

		free( topol[l][0] ); topol[l][0] = NULL;
		free( topol[l][1] ); topol[l][1] = NULL;
		free( topol[l] ); topol[l] = NULL;


//		reporterr(       ">514\n%s\n", aseq[514] );
	}
#if SCOREOUT
	reporterr(       "totalscore = %10.2f\n\n", tscore );
#endif
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( fftlog );
	free( gaplen );
	free( gapmap );
	FreeDoubleMtx( dynamicmtx );
	free( alreadyaligned );
	effarr1 = NULL;
	return( 0 );

	chudan_tbfast:

	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
	A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
	if( effarr1 ) free( effarr1 ); effarr1 = NULL;
	if( effarr2 ) free( effarr2 ); effarr2 = NULL;
	if( indication1 ) free( indication1 ); indication1 = NULL;
	if( indication2 ) free( indication2 ); indication2 = NULL;
	if( fftlog ) free( fftlog ); fftlog = NULL;
	if( gaplen ) free( gaplen ); gaplen = NULL;
	if( gapmap ) free( gapmap ); gapmap = NULL;
	if( alreadyaligned ) free( alreadyaligned ); alreadyaligned = NULL;
	if( dynamicmtx ) FreeDoubleMtx( dynamicmtx ); dynamicmtx = NULL;
#if SKIP
	if( skiptable1 ) FreeIntMtx( skiptable1 ); skiptable1 = NULL;
	if( skiptable2 ) FreeIntMtx( skiptable2 ); skiptable2 = NULL;
#endif

	return( 1 );
}

static void WriteOptions( FILE *fp )
{
	CALLS && printf("called %s:WriteOptions()\n", __FILE__);
	FILES && printf("file write %d %s:%d\n", fp, __FILE__, __LINE__);

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    reporterr(       "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = UPGMA (average).\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree = Minimum linkage.\n" );
	else
		fprintf( fp, "Unknown tree.\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
        else
        {
            if( fftscore )
                fprintf( fp, "Basis : Polarity and Volume\n" );
            else
                fprintf( fp, "Basis : 20 amino acids\n" );
        }
        fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
        fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
    }
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}


int disttbfast( int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv, int (*callback)(int, int, char*))
{
	int  *nlen = NULL;	
	int  *nogaplen = NULL;	
	char **name = NULL, **seq = NULL;
	char **mseq1 = NULL, **mseq2 = NULL;
	char **bseq = NULL;
	double *eff = NULL;
	int i, j;
	int ***topol = NULL;
	int *addmem = NULL;
	Treedep *dep = NULL;
	float **len = NULL;
	FILE *infp = NULL;
//	FILE *adfp;
	char c;
	int alloclen;
	float longer, shorter;
	float lenfac;
	float bunbo;

	FILE *orderfp = NULL, *hat2p = NULL;
	int *grpseq = NULL;
	char *tmpseq = NULL;
	int  **pointt = NULL;
	float **mtx = NULL; // by D. Mathog
	int *table1 = NULL;
	char b[B];
	int ien, nlim;
	int includememberres0, includememberres1;
	double unweightedspscore;
	int alignmentlength;
	char *mergeoralign = NULL;
	int foundthebranch;
	int nsubalignments, maxmem;
	int **subtable = NULL;
	int *insubtable = NULL;
	int *preservegaps = NULL;
	char ***subalnpt = NULL;
	int val;
	char **tmpargv = NULL;
	int iguidetree;
	float *selfscore = NULL;
	int calcpairdists;
	int	**skiptable = NULL;
	char algbackup;


	if( ngui )
	{
		BRANCHES && printf("branch %d\n", __LINE__); // no
		initglobalvariables();
		njob = ngui;
		nlenmax = 0;
		for( i=0; i<njob; i++ )
		{
			ien = strlen( seqgui[i] );
			if( ien > nlenmax ) nlenmax = ien;
		}
		infp = NULL;
//		stderr = fopen( "/dev/null", "a" ); // Windows????
		tmpargv = AllocateCharMtx( argc, 0 );
		for( i=0; i<argc; i++ ) tmpargv[i] = argv[i];
		gmsg = 1;
	}
	else {
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		gmsg = 0; // iranai
	}

	arguments( argc, argv );
	algbackup = alg; // tbfast wo disttbfast ni ketsugou shitatame.
#ifndef enablemultithread
	BRANCHES && printf("branch %d\n", __LINE__); // no
	nthread = 0;
#endif


	if( ngui )
	{
		BRANCHES && printf("branch %d\n", __LINE__); // no (Don't need to mark ngui branches further; they're never taken.)
		for( i=0; i<argc; i++ ) 
		{
//			free( tmpargv[i] );
			argv[i] = tmpargv[i];
		}
		free( tmpargv );
	}
	else
	{
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		if( inputfile )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			infp = fopen( inputfile, "r" );
			FILES && printf("file open r \"%s\" (%d) %s:%d\n", inputfile, infp, __FILE__, __LINE__);
			if( !infp )
			{
				reporterr(       "Cannot open %s\n", inputfile );
				exit( 1 );
			}
		}
		else {
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			infp = stdin;
		}
	
		getnumlen( infp ); // yes
		rewind( infp );
	}
	
	if( njob > 1000000 )
	{
		reporterr(       "The number of sequences must be < %d\n", 1000000 );
		reporterr(       "Please try the --parttree option for such large data.\n" );
		exit( 1 );
	}

	if( njob < 2 )
	{
		reporterr(       "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	if( specificityconsideration != 0.0 && nlenmax)
	{
		BRANCHES && printf("branch %d\n", __LINE__); // no
		if( nlenmax > 100000 )
		{
			reporterr( "\n" );
			reporterr( "Too long to apply --allowshift or --unalignlevel>0\n" );
			reporterr( "Please use the normal mode.\n" );
			reporterr( "Please also note that MAFFT does not assume genomic rearrangements.\n" );
			reporterr( "\n" );
			exit( 1 );
		}
	}


	if( subalignment )
	{
		BRANCHES && printf("branch %d\n", __LINE__); // no
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem );
		reporterr(       "nsubalignments = %d\n", nsubalignments );
		reporterr(       "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		insubtable = AllocateIntVec( njob );
		preservegaps = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		for( i=0; i<njob; i++ ) preservegaps[i] = 0;
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );
	}


	seq = AllocateCharMtx( njob, nlenmax*1+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	eff = AllocateDoubleVec( njob );
	mergeoralign = AllocateCharVec( njob );

	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );

	if( nadd ) addmem = AllocateIntVec( nadd+1 );

#if 0
	Read( name, nlen, seq );
	readData( infp, name, nlen, seq );
#else
    name = AllocateCharMtx( njob, B+1 );
    nlen = AllocateIntVec( njob ); 
    nogaplen = AllocateIntVec( njob ); 
	if( ngui )
	{
		if( copydatafromgui( namegui, seqgui, name, nlen, seq ) )
			exit( 1 );
	}
	else
	{
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		readData_pointer( infp, name, nlen, seq );
		FILES && printf("file close %d %s:%d\n", infp, __FILE__, __LINE__);
		fclose( infp );
	}
#endif


	constants( njob, seq );


#if 0
	reporterr(       "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	c = seqcheck( seq );
	if( c )
	{
		reporterr(       "Illegal character %c\n", c );
		exit( 1 );
	}

	reporterr(       "\n" );

	reporterr(       "tuplesize = %d, dorp = %c\n", tuplesize, dorp );
	if( dorp == 'p' && tuplesize != 6 )
	{
		reporterr(       "tuplesize must be 6 for aa sequence\n" );
		exit( 1 );
	}
	if( dorp == 'd' && tuplesize != 6 && tuplesize != 10 )
	{
		reporterr(       "tuplesize must be 6 or 10 for dna sequence\n" );
		exit( 1 );
	}

	if( !treein )
	{
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		reporterr(       "\n\nMaking a distance matrix ..\n" );
		if( callback && callback( 0, 0, "Distance matrix" ) ) {
			BRANCHES && printf("branch %d\n", __LINE__); // no
			goto chudan;
		}

	    tmpseq = AllocateCharVec( nlenmax+1 );
		grpseq = AllocateIntVec( nlenmax+1 );
		pointt = AllocateIntMtx( njob, nlenmax+1 );
    	mtx = AllocateFloatHalfMtx( njob ); 
		if( dorp == 'd' ) tsize = (int)pow( 4, tuplesize );
		else              tsize = (int)pow( 6, 6 );

		if( dorp == 'd' && tuplesize == 6 )
		{
			lenfaca = D6LENFACA;
			lenfacb = D6LENFACB;
			lenfacc = D6LENFACC;
			lenfacd = D6LENFACD;
		}
		else if( dorp == 'd' && tuplesize == 10 )
		{
			lenfaca = D10LENFACA;
			lenfacb = D10LENFACB;
			lenfacc = D10LENFACC;
			lenfacd = D10LENFACD;
		}
		else    
		{
			lenfaca = PLENFACA;
			lenfacb = PLENFACB;
			lenfacc = PLENFACC;
			lenfacd = PLENFACD;
		}

		maxl = 0;
		for( i=0; i<njob; i++ ) 
		{
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			gappick0( tmpseq, seq[i] );
			nogaplen[i] = strlen( tmpseq );
			if( nogaplen[i] < 6 )
			{
					BRANCHES && printf("branch %d\n", __LINE__); // no
//				reporterr(       "Seq %d, too short, %d characters\n", i+1, nogaplen[i] );
//				reporterr(       "Please use mafft-ginsi, mafft-linsi or mafft-ginsi\n\n\n" );
//				exit( 1 );
			}
			if( nogaplen[i] > maxl ) maxl = nogaplen[i];
			if( dorp == 'd' ) /* nuc */
			{
				BRANCHES && printf("branch %d\n", __LINE__); // yes
				seq_grp_nuc( grpseq, tmpseq );
//				makepointtable_nuc( pointt[i], grpseq );
//				makepointtable_nuc_octet( pointt[i], grpseq );
				if( tuplesize == 10 )
					makepointtable_nuc_dectet( pointt[i], grpseq );
				else if( tuplesize == 6 ) {
					BRANCHES && printf("branch %d\n", __LINE__);
					makepointtable_nuc( pointt[i], grpseq );
				} else
				{
					reporterr(       "tuplesize=%d: not supported\n", tuplesize );
					exit( 1 );
				}
			}
			else                 /* amino */
			{
				BRANCHES && printf("branch %d\n", __LINE__);  // no
				seq_grp( grpseq, tmpseq );
				makepointtable( pointt[i], grpseq );
			}
		}
		if( nunknown ) reporterr(       "\nThere are %d ambiguous characters.\n", nunknown );
#ifdef enablemultithread
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		if( nthread > 0 )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			distancematrixthread_arg_t *targ; 
			int jobpos;
			pthread_t *handle;
			pthread_mutex_t mutex;

			jobpos = 0; 
			targ = calloc( nthread, sizeof( distancematrixthread_arg_t ) ); 
			handle = calloc( nthread, sizeof( pthread_t ) ); 
			pthread_mutex_init( &mutex, NULL );

			for( i=0; i<nthread; i++ )
			{
				targ[i].thread_no = i;
				targ[i].njob = njob;
				targ[i].jobpospt = &jobpos;
				targ[i].pointt = pointt;
				targ[i].mtx = mtx;
				targ[i].mutex = &mutex;

				pthread_create( handle+i, NULL, distancematrixthread, (void *)(targ+i) );
			}
		
			for( i=0; i<nthread; i++ )
			{
				pthread_join( handle[i], NULL );
			}
			pthread_mutex_destroy( &mutex );
			free( handle );
			free( targ );
		}
		else
#endif
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		{
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			// njob == number of input sequences, in at least one case.
			for( i=0; i<njob; i++ )
			{
				BRANCHES && printf("branch %d\n", __LINE__); // yes
				table1 = (int *)calloc( tsize, sizeof( int ) );
				if( !table1 ) ErrorExit( "Cannot allocate table1\n" );
				if( i % 100 == 0 )
				{
					BRANCHES && printf("branch %d\n", __LINE__); // yes
					reporterr(       "\r% 5d / %d", i+1, njob );
					if( callback && callback( 0, i*25/njob, "Distance matrix" ) ) goto chudan;
				}
				makecompositiontable_p( table1, pointt[i] );
		
				for( j=i; j<njob; j++ ) 
				{
					mtx[i][j-i] = (float)commonsextet_p( table1, pointt[j] );
				} 
				free( table1 ); table1 = NULL;
			}
		}
		reporterr(       "\ndone.\n\n" );
		ien = njob-1;

		for( i=0; i<ien; i++ )
		{
			for( j=i+1; j<njob; j++ ) 
			{
				if( nogaplen[i] > nogaplen[j] )
				{
					longer=(float)nogaplen[i];
					shorter=(float)nogaplen[j];
				}
				else
				{
					longer=(float)nogaplen[j];
					shorter=(float)nogaplen[i];
				}
//				if( tuplesize == 6 )
				lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//				else
//					lenfac = 1.0;
//				reporterr(       "lenfac = %f (%.0f,%.0f)\n", lenfac, longer, shorter );
				bunbo = MIN( mtx[i][0], mtx[j][0] );
				if( bunbo == 0.0 )
					mtx[i][j-i] = 2.0; // 2013/Oct/17 -> 2bai
				else
					mtx[i][j-i] = ( 1.0 - mtx[i][j-i] / bunbo ) * lenfac * 2.0; // 2013/Oct/17 -> 2bai
//				reporterr(       "##### mtx = %f, mtx[i][0]=%f, mtx[j][0]=%f, bunbo=%f\n", mtx[i][j-i], mtx[i][0], mtx[j][0], bunbo );
			}
		}
		if( disopt )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			for( i=0; i<njob; i++ ) 
			{
				sprintf( b, "=lgth = %04d", nogaplen[i] );
				strins( b, name[i] );
			}
		}
		BRANCHES && printf("branch %d\n", __LINE__);
		free( grpseq ); grpseq = NULL;
		free( tmpseq ); tmpseq = NULL;
		FreeIntMtx( pointt ); pointt = NULL;
		commonsextet_p( NULL, NULL );

#if 0 // writehat2 wo kakinaosu -> iguidetree loop nai ni idou
		if( distout )
		{
			hat2p = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
			fclose( hat2p );
		}
#endif

	}
#if 0 
	else 
	{
		reporterr(       "Loading 'hat2' ... " );
		prep = fopen( "hat2", "r" );
		if( prep == NULL ) ErrorExit( "Make hat2." );
		readhat2_float( prep, njob, name, mtx ); // name chuui
		fclose( prep );
		reporterr(       "done.\n" );
	}
#endif


	for( iguidetree=0; iguidetree<nguidetree; iguidetree++ )
//	for( iguidetree=0; ; iguidetree++ )
	{
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		alg = algbackup; // tbfast wo disttbfast ni ketsugou shitatame.



		topol = AllocateIntCub( njob, 2, 0 );
		len = AllocateFloatMtx( njob, 2 );

		if( iguidetree == nguidetree - 1 ) calcpairdists = 0;
		else                               calcpairdists = 1;


		if( calcpairdists ) selfscore = AllocateFloatVec( njob );

		if( callback && callback( 0, 25, "Guide tree" ) ) {
			BRANCHES && printf("branch %d\n", __LINE__); // no
			goto chudan;
		}
	
		if( treein )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			nguidetree = 1; //  iranai
			reporterr(       "Loading a tree ... " );
			loadtree( njob, topol, len, name, nogaplen, dep );
//			loadtop( njob, topol, len, name, NULL, dep ); // 2015/Jan/13, not yet checked
		}
		else if( topin )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			reporterr(       "Loading a topology ... " );
			reporterr(       "--topin has been disabled\n" );
			exit( 1 );
//			loadtop( njob, mtx, topol, len );
//			FreeFloatHalfMtx( mtx, njob );
		}
		else
		{
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			if( distout )
			{
				BRANCHES && printf("branch %d\n", __LINE__); // no
				hat2p = fopen( "hat2", "w" );
				FILES && printf("file open w \"hat2\" (%d) %s:%d\n", hat2p, __FILE__, __LINE__);
				WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, mtx );
				// writehat2 wo kakinaosu
				fclose( hat2p );
			}
			if( subalignment ) // merge error no tame
			{
				BRANCHES && printf("branch %d\n", __LINE__); // no
				reporterr(       "Constructing a UPGMA tree ... " );
				fixed_supg_float_realloc_nobk_halfmtx_treeout_constrained( njob, mtx, topol, len, name, nlen, dep, nsubalignments, subtable, !calcpairdists );
				if( !calcpairdists ) 
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else if( treeout ) // merge error no tame
			{
				BRANCHES && printf("branch %d\n", __LINE__); // no
				reporterr(       "Constructing a UPGMA tree ... " );
				fixed_musclesupg_float_realloc_nobk_halfmtx_treeout( njob, mtx, topol, len, name, nogaplen, dep, !calcpairdists );
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
			else
			{
				BRANCHES && printf("branch %d\n", __LINE__); // yes
				reporterr(       "Constructing a UPGMA tree ... " );
				fixed_musclesupg_float_realloc_nobk_halfmtx( njob, mtx, topol, len, dep, 1, !calcpairdists );
				if( !calcpairdists )
				{
					FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
				}
			}
		}
//		else 
//			ErrorExit( "Incorrect tree\n" );
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 50, "Guide tree" ) ) {
			BRANCHES && printf("branch %d\n", __LINE__); // no
			goto chudan;
		}

		if( sparsepickup && iguidetree == nguidetree-1 )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			reporterr(       "Sparsepickup! \n" );
			pickup( njob, nogaplen, topol, name, seq );
			reporterr(       "done. \n" );
			SHOWVERSION;
			goto chudan;
		}

	
		orderfp = fopen( "order", "w" );
		FILES && printf("file open w \"order\" (%d) %s:%d\n", orderfp, __FILE__, __LINE__);
		if( !orderfp )
		{
			reporterr(       "Cannot open 'order'\n" );
			exit( 1 );
		}
		for( i=0; (j=topol[njob-2][0][i])!=-1; i++ )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			FILES && printf("file write \"order\" %s:%d\n", __FILE__, __LINE__);
			fprintf( orderfp, "%d\n", j );
		}
		for( i=0; (j=topol[njob-2][1][i])!=-1; i++ )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			FILES && printf("file write \"order\" %s:%d\n", __FILE__, __LINE__);
			fprintf( orderfp, "%d\n", j );
		}
		fclose( orderfp );


	
		if( ( treeout || distout )  && noalign ) 
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			writeData_pointer( stdout, njob, name, nlen, seq );
			reporterr(       "\n" );
			SHOWVERSION;
			goto chudan;
//			return( 0 );
		}
		
		if( tbrweight )
		{
			weight = 3; 
#if 0
			utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			counteff_simple_float_nostatic( njob, topol, len, eff );
#endif
		}
		else
		{
			for( i=0; i<njob; i++ ) eff[i] = 1.0;
		}
	
#if 0
		for( i=0; i<njob; i++ )
			reporterr(       "eff[%d] = %20.16f\n", i, eff[i] );
		exit( 1 );
#endif
	
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		FreeFloatMtx( len ); len = NULL;
	
		bseq = AllocateCharMtx( njob, nlenmax*2+1 );
		alloclen = nlenmax*2+1;
	
	
		if( nadd )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			alignmentlength = strlen( seq[0] );
			for( i=0; i<njob-nadd; i++ )
			{
				if( alignmentlength != strlen( seq[i] ) )
				{
					reporterr(       "#################################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# The original %d sequences must be aligned\n", njob-nadd );
					reporterr(       "# alignmentlength = %d, but strlen(seq[%d])=%d\n", alignmentlength, i, (int)strlen( seq[i] ) );
					reporterr(       "#################################################################################\n" );
					exit( 1 );
				}
			}
			if( addprofile )
			{
				BRANCHES && printf("branch %d\n", __LINE__); // no
				alignmentlength = strlen( seq[njob-nadd] );
				for( i=njob-nadd; i<njob; i++ )
				{
					if( alignmentlength != strlen( seq[i] ) )
					{
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# The %d additional sequences must be aligned\n", nadd );
						reporterr(       "# Otherwise, try the '--add' option, instead of '--addprofile' option.\n" );
						reporterr(       "###############################################################################\n" );
						exit( 1 );
					}
				}
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				foundthebranch = 0;
				for( i=0; i<njob-1; i++ )
				{
					BRANCHES && printf("branch %d\n", __LINE__); // no
					if( samemember( topol[i][0], addmem ) ) // jissainiha nai
					{
						mergeoralign[i] = '1';
						foundthebranch = 1;
					}
					else if( samemember( topol[i][1], addmem ) ) // samemembern ni henkou kanou
					{
						mergeoralign[i] = '2';
						foundthebranch = 1;
					}
					else
					{
						mergeoralign[i] = 'n';
					}
				}
				if( !foundthebranch )
				{
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# There is no appropriate position to add the %d sequences in the guide tree.\n", nadd );
					reporterr(       "# Check whether the %d sequences form a monophyletic cluster.\n", nadd );
					reporterr(       "# If not, try the '--add' option, instead of the '--addprofile' option.\n" );
					reporterr(       "############################################################################### \n" );
					exit( 1 );
				}
				BRANCHES && printf("branch %d\n", __LINE__); // no
				commongappick( nadd, seq+njob-nadd );
				for( i=njob-nadd; i<njob; i++ ) strcpy( bseq[i], seq[i] );
			}
			else
			{
				for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'n';
				for( j=njob-nadd; j<njob; j++ )
				{
					addmem[0] = j;
					addmem[1] = -1;
					for( i=0; i<njob-1; i++ )
					{
						BRANCHES && printf("branch %d\n", __LINE__); // no
						if( samemembern( topol[i][0], addmem, 1 ) ) // arieru
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '1';
						}
						else if( samemembern( topol[i][1], addmem, 1 ) )
						{
//							reporterr(       "HIT!\n" );
							if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
							else mergeoralign[i] = '2';
						}
					}
				}
		
				for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
				addmem[nadd] = -1;
				nlim = njob-1;
//				for( i=0; i<njob-1; i++ )
				for( i=0; i<nlim; i++ )
				{
					BRANCHES && printf("branch %d\n", __LINE__); // no
					includememberres0 = includemember( topol[i][0], addmem );
					includememberres1 = includemember( topol[i][1], addmem );
//					if( includemember( topol[i][0], addmem ) && includemember( topol[i][1], addmem ) )
					if( includememberres0 && includememberres1 )
					{
						mergeoralign[i] = 'w';
					}
					else if( includememberres0 )
					{
						mergeoralign[i] = '1';
					}
					else if( includememberres1 )
					{
						mergeoralign[i] = '2';
					}
				}
#if 0
				for( i=0; i<njob-1; i++ )
				{
					reporterr(       "mem0 = " );
					for( j=0; topol[i][0][j]>-1; j++ )	reporterr(       "%d ", topol[i][0][j] );
					reporterr(       "\n" );
					reporterr(       "mem1 = " );
					for( j=0; topol[i][1][j]>-1; j++ )	reporterr(       "%d ", topol[i][1][j] );
					reporterr(       "\n" );
					reporterr(       "i=%d, mergeoralign[] = %c\n", i, mergeoralign[i] );
				}
#endif
				for( i=njob-nadd; i<njob; i++ ) {
					BRANCHES && printf("branch %d\n", __LINE__); // no
					gappick0( bseq[i], seq[i] );
				}
			}

			BRANCHES && printf("branch %d\n", __LINE__); // no
			commongappick( njob-nadd, seq );
			for( i=0; i<njob-nadd; i++ ) strcpy( bseq[i], seq[i] );
		}
//--------------- kokokara ----
		else if( subalignment )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
			for( i=0; i<nsubalignments; i++ )
			{
				reporterr(       "Checking subalignment %d:\n", i+1 );
				alignmentlength = strlen( seq[subtable[i][0]] );
//				for( j=0; subtable[i][j]!=-1; j++ )
//					reporterr(       " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
				for( j=0; subtable[i][j]!=-1; j++ )
				{
					if( subtable[i][j] >= njob ) // check sumi
					{
						reporterr(       "No such sequence, %d.\n", subtable[i][j]+1 );
						exit( 1 );
					}
					if( alignmentlength != strlen( seq[subtable[i][j]] ) )
					{
						reporterr(       "\n" );
						reporterr(       "###############################################################################\n" );
						reporterr(       "# ERROR!\n" );
						reporterr(       "# Subalignment %d must be aligned.\n", i+1 );
						reporterr(       "# Please check the alignment lengths of following sequences.\n" );
						reporterr(       "#\n" );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][0]+1, name[subtable[i][0]]+1, alignmentlength );
						reporterr(       "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j]+1, name[subtable[i][j]]+1, (int)strlen( seq[subtable[i][j]] ) );
						reporterr(       "#\n" );
						reporterr(       "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" );
						if( subalignmentoffset )
						{
							reporterr(       "#\n" );
							reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
							reporterr(       "# In this case, the rule of numbering is:\n" );
							reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
							reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
						}
						reporterr(       "###############################################################################\n" );
						reporterr(       "\n" );
						exit( 1 );
					}
					insubtable[subtable[i][j]] = 1;
				}
				for( j=0; j<njob-1; j++ )
				{
					BRANCHES && printf("branch %d\n", __LINE__); // no
					if( includemember( topol[j][0], subtable[i] ) && includemember( topol[j][1], subtable[i] ) )
					{
						mergeoralign[j] = 'n';
					}
				}
				foundthebranch = 0;
				for( j=0; j<njob-1; j++ )
				{
					BRANCHES && printf("branch %d\n", __LINE__); // no
					if( samemember( topol[j][0], subtable[i] ) || samemember( topol[j][1], subtable[i] ) )
					{
						foundthebranch = 1;
						reporterr(       " -> OK\n" );
						break;
					}
				}
				if( !foundthebranch )
				{
					BRANCHES && printf("branch %d\n", __LINE__); // no
					FILES && printf("file cp \"infile.tree\" \"GuideTree\" %s:%d\n", __FILE__, __LINE__);
					system( "cp infile.tree GuideTree" ); // tekitou
					reporterr(       "\n" );
					reporterr(       "###############################################################################\n" );
					reporterr(       "# ERROR!\n" );
					reporterr(       "# Subalignment %d does not seem to form a monophyletic cluster\n", i+1 );
					reporterr(       "# in the guide tree ('GuideTree' in this directory) internally computed.\n" );
					reporterr(       "# If you really want to use this subalignment, pelase give a tree with --treein \n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/treein.html\n" );
					reporterr(       "# http://mafft.cbrc.jp/alignment/software/merge.html\n" );
					if( subalignmentoffset )
					{
						reporterr(       "#\n" );
						reporterr(       "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
						reporterr(       "# In this case, the rule of numbering is:\n" );
						reporterr(       "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
						reporterr(       "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
					}
					reporterr(       "############################################################################### \n" );
					reporterr(       "\n" );
					exit( 1 );
				}
//				commongappick( seq[subtable[i]], subalignment[i] ); // irukamo
			}
#if 0
			for( i=0; i<njob-1; i++ )
			{
				reporterr(       "STEP %d\n", i+1 );
				reporterr(       "group1 = " );
				for( j=0; topol[i][0][j] != -1; j++ )
					reporterr(       "%d ", topol[i][0][j]+1 );
				reporterr(       "\n" );
				reporterr(       "group2 = " );
				for( j=0; topol[i][1][j] != -1; j++ )
					reporterr(       "%d ", topol[i][1][j]+1 );
				reporterr(       "\n" );
				reporterr(       "%d -> %c\n\n", i, mergeoralign[i] );
			}
#endif
			BRANCHES && printf("branch %d\n", __LINE__); // no
			for( i=0; i<njob; i++ ) 
			{
				if( insubtable[i] ) strcpy( bseq[i], seq[i] );
				else gappick0( bseq[i], seq[i] );
			}
	
			for( i=0; i<nsubalignments; i++ ) 
			{
				for( j=0; subtable[i][j]!=-1; j++ ) subalnpt[i][j] = bseq[subtable[i][j]];
				if( !preservegaps[i] ) {
					BRANCHES && printf("branch %d\n", __LINE__); // no
					commongappick( j, subalnpt[i] );
				}
			}
	
#if 0 // --> iguidetree loop no soto he
			FreeIntMtx( subtable );
			free( insubtable );
			for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
			free( subalnpt );
			free( preservegaps );
#endif
		}
//--------------- kokomade ----
		else
		{
			for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
		}

		if ( calcpairdists ) {
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			for ( i=0; i<njob; i++ ) {
				selfscore[i] = naivepairscore11( seq[i], seq[i], penalty_dist );
			}
		}
	
		reporterr(       "Progressive alignment %d/%d... \n", iguidetree+1, nguidetree );
	
#ifdef enablemultithread
		if( nthread > 0 && nadd == 0 )
		{
			BRANCHES && printf("branch %d\n", __LINE__); // no
			treebasethread_arg_t *targ; 
			int jobpos;
			pthread_t *handle;
			pthread_mutex_t mutex;
			pthread_cond_t treecond;
			int *fftlog;
			int nrun;
			int nthread_yoyu;
	
			nthread_yoyu = nthread * 1;
			nrun = 0;
			jobpos = 0; 
			targ = calloc( nthread_yoyu, sizeof( treebasethread_arg_t ) ); 
			fftlog = AllocateIntVec( njob );
			handle = calloc( nthread_yoyu, sizeof( pthread_t ) ); 
			pthread_mutex_init( &mutex, NULL );
			pthread_cond_init( &treecond, NULL );
	
			for( i=0; i<njob; i++ ) dep[i].done = 0; 
			for( i=0; i<njob; i++ ) fftlog[i] = 1; 
	
			for( i=0; i<nthread_yoyu; i++ )
			{
				targ[i].thread_no = i;
				targ[i].njob = njob;
				targ[i].nrunpt = &nrun;
				targ[i].nlen = nlen;
				targ[i].jobpospt = &jobpos;
				targ[i].topol = topol;
				targ[i].dep = dep;
				targ[i].aseq = bseq;
				targ[i].effarr = eff;
				targ[i].alloclenpt = &alloclen;
				targ[i].fftlog = fftlog;
				targ[i].mergeoralign = mergeoralign;
#if 1 // tsuneni SEPARATELYCALCPAIRDISTS
				targ[i].newdistmtx = NULL;
				targ[i].selfscore = NULL;
#else
				if( calcpairdists ) // except for last cycle
				{
					targ[i].newdistmtx = mtx;
					targ[i].selfscore = selfscore;
				}
				else
				{
					targ[i].newdistmtx = NULL;
					targ[i].selfscore = NULL;
				}
#endif
				targ[i].mutex = &mutex;
				targ[i].treecond = &treecond;
	
				pthread_create( handle+i, NULL, treebasethread, (void *)(targ+i) );
			}
			
			for( i=0; i<nthread_yoyu; i++ )
			{
				pthread_join( handle[i], NULL );
			}
			pthread_mutex_destroy( &mutex );
			pthread_cond_destroy( &treecond );
			free( handle );
			free( targ );
			free( fftlog );
		}
		else
#endif
		{
#if 0
			if( calcpairdists ) // except for last
			{
				if( treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, mtx, selfscore, &alloclen, callback ) ) goto chudan;
			}
			else
#endif
			{
				BRANCHES && printf("branch %d\n", __LINE__); // yes
				if( treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, NULL, NULL, &alloclen, callback ) ) goto chudan;
			}
		}
		reporterr(       "\ndone.\n\n" );
		if( callback && callback( 0, 100, "Progressive alignment" ) ) goto chudan;
		free( topol[njob-1][0] ); topol[njob-1][0]=NULL;
		free( topol[njob-1][1] ); topol[njob-1][1]=NULL;
		free( topol[njob-1] ); topol[njob-1]=NULL;
		free( topol ); topol=NULL;



#if 1
// Distance matrix from MSA SEPARATELYCALCPAIRDISTS
//		if( iguidetree < nguidetree-1 )
#ifdef enablemultithread
//		if( nthread>0 && nadd==0 ) if( calcpairdists )
		if( calcpairdists )
#else
//		if( 0 && nadd==0 ) if( calcpairdists ) // zettai nai
		if( calcpairdists )
#endif
		{
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			reporterr( "Making a distance matrix from msa.. \n" );
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, bseq ); // allocate suru.
#ifdef enablemultithread
			if( nthread > 0 )
			{
				BRANCHES && printf("branch %d\n", __LINE__); // no
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
	
				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( nthread, sizeof( msadistmtxthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
	
				for( i=0; i<nthread; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].iscore = mtx;
					targ[i].seq = bseq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, msadistmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthread; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
//				reporterr( "Check source!\n" );
//				exit( 1 );

#if 1
				msadistmtxthread_arg_t *targ;
				Jobtable jobpos;

				jobpos.i = 0;
				jobpos.j = 0;
	
				targ = calloc( 1, sizeof( msadistmtxthread_arg_t ) );
	
				{
					BRANCHES && printf("branch %d\n", __LINE__); // yes
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].iscore = mtx;
					targ[0].seq = bseq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
	
					msadistmtxthread( targ );
				}
	
				free( targ );
#endif
			}
			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;
			reporterr(       "\ndone.\n\n" );
		}
// Distance matrix from MSA end
#endif

		if( calcpairdists ) 
		{
			BRANCHES && printf("branch %d\n", __LINE__); // yes
			free( selfscore );
			selfscore = NULL;
			FreeCharMtx( bseq );
			bseq = NULL;
		}
	}
#if DEBUG
	reporterr(       "closing trap_g\n" );
#endif
//	fclose( trap_g );

	if( scoreout )
	{
		BRANCHES && printf("branch %d\n", __LINE__); // no
		unweightedspscore = plainscore( njob, bseq );
		reporterr(       "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore );
		reporterr(       "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
		reporterr(       "\n\n" );
	}

#if DEBUG
	reporterr(       "writing alignment to stdout\n" );
#endif


	val = 0;
	if( ngui ) 
	{
		ien = strlen( bseq[0] );
		if( ien > lgui )
		{
			reporterr( "alignmentlength = %d, gui allocated %d", ien, lgui );
			val = GUI_LENGTHOVER;
		}
		else
		{
			for( i=0; i<njob; i++ ) 
			{
#if 1
				strcpy( seqgui[i], bseq[i] );
#else
				free( seqgui[i] );
				seqgui[i] =  bseq[i];
#endif
			}
		}
	}
	else
	{
		BRANCHES && printf("branch %d\n", __LINE__); // yes
		writeData_pointer( stdout, njob, name, nlen, bseq );
	} 

	if( spscoreout ) reporterr( "Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore( njob, bseq ) );
	SHOWVERSION;

	if( subalignment )
	{
		BRANCHES && printf("branch %d\n", __LINE__); // no
		FreeIntMtx( subtable );
		free( insubtable );
		for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
		free( subalnpt );
		free( preservegaps );
	}

#if 1 // seqgui[i] =  bseq[i] no toki bseq ha free shinai
	BRANCHES && printf("branch %d\n", __LINE__); // yes
	FreeCharMtx( bseq );
#endif
	FreeCharMtx( name );
    free( nlen );

	free( mergeoralign );
	FreeCharMtx( seq );
    free( nogaplen );

	free( mseq1 );
	free( mseq2 );
//	FreeIntCub( topol ); // 
//	FreeFloatMtx( len ); //
//	free( mergeoralign ); //
	free( dep );

	if( nadd ) free( addmem );
	free( eff );
	BRANCHES && printf("branch %d\n", __LINE__); // yes
	freeconstants();
	closeFiles();
	FreeCommonIP();
	return( val );

chudan:
	// This seems like some kind of cleanup routine.
	BRANCHES && printf("branch chudan\n");
	BRANCHES && printf("branch %d\n", __LINE__); // no

	if( nlen ) free( nlen ); nlen = NULL;
	if( seq ) FreeCharMtx( seq ); seq = NULL;
	if( mseq1 ) free( mseq1 ); mseq1 = NULL;
	if( mseq2 ) free( mseq2 ); mseq2 = NULL;
	if( topol ) 
	{
		BRANCHES && printf("branch %d\n", __LINE__); // no
		for( i=0; i<njob; i++ )
		{
			if( topol[i] && topol[i][0] ) 
			{
				free( topol[i][0] ); topol[i][0] = NULL;
			}
			if( topol[i] && topol[i][1] ) 
			{
				free( topol[i][1] ); topol[i][1] = NULL;
			}
			if( topol[i] ) free( topol[i] ); topol[i] = NULL;
		}
		free( topol ); topol = NULL;
	}
	if( len ) FreeFloatMtx( len ); len = NULL;
	if( eff ) free( eff ); eff = NULL;
	if( mergeoralign ) free( mergeoralign ); mergeoralign = NULL;
	if( dep ) free( dep ); dep = NULL;
	if( addmem ) free( addmem ); addmem = NULL;
	if( name ) FreeCharMtx( name ); name = NULL;
	if( nogaplen ) free( nogaplen ); nogaplen = NULL;

	if( tmpseq ) free( tmpseq ); tmpseq = NULL;
	if( grpseq ) free( grpseq ); grpseq = NULL;
	if( pointt ) FreeIntMtx( pointt ); pointt = NULL;
	if( mtx ) FreeFloatHalfMtx( mtx, njob ); mtx = NULL;
	if( table1 ) free( table1 ); table1 = NULL;

	if( bseq ) FreeCharMtx( bseq ); bseq = NULL;
	if( selfscore ) free( selfscore ); selfscore = NULL;
	if( skiptable ) FreeIntMtx( skiptable ); skiptable = NULL;

	freeconstants();
	closeFiles();
	FreeCommonIP();

	return( GUI_CANCEL );
}

int main( int argc, char **argv )
{
	int res = disttbfast( 0, 0, NULL, NULL, argc, argv, NULL );
	if( res == GUI_CANCEL ) res = 0; // treeout de goto chudan wo riyousuru
	return res;
}
