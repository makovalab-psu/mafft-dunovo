#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "mtxutl.h"

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

void MtxuntDouble( double **mtx, int n ) // yes
{
		CALLS && printf("called %s:MtxuntDouble()\n", __FILE__);
    int i, j;
    for( i=0; i<n; i++ ) for( j=0; j<n; j++ ) mtx[i][j] = 0.0;
    for( i=0; i<n; i++ ) mtx[i][i] = 1.0;
}

void MtxmltDouble( double **mtx1, double **mtx2, int n ) // yes
{
		CALLS && printf("called %s:MtxmltDouble()\n", __FILE__);
    int i, j, k;
    double s, *tmp;

	tmp = (double *)calloc( n, sizeof( double ) );
    for( i=0; i<n; i++ ) 
    {
        for( k=0; k<n; k++ ) tmp[k] = mtx1[i][k];
        for( j=0; j<n; j++ ) 
        {
            s = 0.0;
            for( k=0; k<n; k++ ) s += tmp[k] * mtx2[k][j];
            mtx1[i][j] = s;
        }
    }
	free( tmp );
}

char *AllocateCharVec( int l1 ) // yes
{
	CALLS && printf("called %s:AllocateCharVec()\n", __FILE__);
	char *cvec;
	
	cvec = (char *)calloc( l1, sizeof( char ) );
	if( cvec == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Cannot allocate %d character vector.\n", l1 );
		exit( 1 );
	}
	return( cvec );
}

void ReallocateCharMtx( char **mtx, int l1, int l2 )
{
	CALLS && printf("called %s:ReallocateCharMtx()\n", __FILE__);
	int i;
	for( i=0; i<l1; i++ ) 
	{
		mtx[i] = (char *)realloc( mtx[i], (l2+1) * sizeof( char ) );
		if( mtx[i] == NULL )
		{
			fprintf( stderr, "Cannot reallocate %d x %d character matrix.\n", l1, l2 );
		}
	}
}

char **AllocateCharMtx( int l1, int l2 ) // yes
{
	CALLS && printf("called %s:AllocateCharMtx()\n", __FILE__);
	int i;
	char **cmtx;
	
	cmtx = (char **)calloc( l1+1, sizeof( char * ) );
	if( cmtx == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Cannot allocate %d x %d character matrix.\n", l1, l2 );
		exit( 1 );
	}   
	if( l2 )
	{
		for( i=0; i<l1; i++ ) 
		{
			cmtx[i] = AllocateCharVec( l2 );
		}
	}
	cmtx[l1] = NULL;
	return( cmtx );
} 

void FreeCharMtx( char **mtx ) // yes
{
	CALLS && printf("called %s:FreeCharMtx()\n", __FILE__);
/*
	char **x;
	x = mtx;
	while( *x != NULL ) free( *x++ );
	free( mtx );
*/
	int i;
	for( i=0; mtx[i]; i++ ) 
	{
		free( mtx[i] );
	}
	free( mtx );
}

float *AllocateFloatVec( int l1 ) // yes
{
	CALLS && printf("called %s:AllocateFloatVec()\n", __FILE__);
	float *vec;

	vec = (float *)calloc( (unsigned int)l1, sizeof( float ) );
	if( vec == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Allocation error ( %d fload vec )\n", l1 );
		exit( 1 );
	}
	return( vec );
}

void FreeFloatVec( float *vec ) // yes
{
	CALLS && printf("called %s:FreeFloatVec()\n", __FILE__);
	free( (char *)vec );
}

float **AllocateFloatHalfMtx( int ll1 ) // yes
{
	CALLS && printf("called %s:AllocateFloatHalfMtx()\n", __FILE__);
	float **mtx;
	int i;

	mtx = (float **)calloc( (unsigned int)ll1+1, sizeof( float * ) );
	if( mtx == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Allocation error ( %d fload halfmtx )\n", ll1 );
		exit( 1 );
	}
	for( i=0; i<ll1; i++ )
	{
		mtx[i] = (float *)calloc( ll1-i, sizeof( float ) );
		if( !mtx[i] )
		{
			FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
			fprintf( stderr, "Allocation error( %d floathalfmtx )\n", ll1 );
			exit( 1 );
		}
	}
	mtx[ll1] = NULL;
	return( mtx );
}

float **AllocateFloatMtx( int ll1, int ll2 ) // yes
{
	CALLS && printf("called %s:AllocateFloatMtx()\n", __FILE__);
	float **mtx;
	int i;

	mtx = (float **)calloc( (unsigned int)ll1+1, sizeof( float * ) );
	if( mtx == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Allocation error ( %d x %d fload mtx )\n", ll1, ll2 );
		exit( 1 );
	}
	if( ll2 )
	{
		for( i=0; i<ll1; i++ )
		{
			mtx[i] = (float *)calloc( ll2, sizeof( float ) );
			if( !mtx[i] )
			{
				FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
				fprintf( stderr, "Allocation error( %d x %d floatmtx )\n", ll1, ll2 );
				exit( 1 );
			}
		}
	}
	mtx[ll1] = NULL;
	return( mtx );
}

void FreeFloatHalfMtx( float **mtx, int n ) // yes
{
	CALLS && printf("called %s:FreeFloatHalfMtx()\n", __FILE__);
	int i;

	for( i=0; i<n; i++ ) 
	{
		if( mtx[i] ) FreeFloatVec( mtx[i] ); mtx[i] = NULL;
	}
	free( mtx );
}
void FreeFloatMtx( float **mtx ) // yes
{
	CALLS && printf("called %s:FreeFloatMtx()\n", __FILE__);
	int i;

	for( i=0; mtx[i]; i++ ) 
	{
		if( mtx[i] ) FreeFloatVec( mtx[i] ); mtx[i] = NULL;
	}
	free( mtx );
}

int *AllocateIntVec( int ll1 ) // yes
{
	CALLS && printf("called %s:AllocateIntVec()\n", __FILE__);
	int *vec;

	vec = (int *)calloc( ll1, sizeof( int ) );
	if( vec == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Allocation error( %d int vec )\n", ll1 );
		exit( 1 );
	}
	return( vec );
}	

void FreeIntVec( int *vec ) // yes
{
	CALLS && printf("called %s:FreeIntVec()\n", __FILE__);
	free( (char *)vec );
}

float **AllocateFloatTri( int ll1 )
{
	CALLS && printf("called %s:AllocateFloatTri()\n", __FILE__);
	float **tri;
	int i;

	tri = (float **)calloc( (unsigned int)ll1+1, sizeof( float * ) );
	if( !tri )
	{
		fprintf( stderr, "Allocation error ( float tri )\n" );
		exit( 1 );
	}
	for( i=0; i<ll1; i++ ) 
	{
		tri[i] = AllocateFloatVec( i+3 );
	}
	tri[ll1] = NULL;
		
	return( tri );
}

void FreeFloatTri( float **tri )
{
	CALLS && printf("called %s:FreeFloatTri()\n", __FILE__);
/*
	float **x;
	x = tri;
	while( *tri != NULL ) free( *tri++ );
	free( x );
*/
	int i;
	for( i=0; tri[i]; i++ ) 
		free( tri[i] );
	free( tri );
}
		
int **AllocateIntMtx( int ll1, int ll2 ) // yes
{
	CALLS && printf("called %s:AllocateIntMtx()\n", __FILE__);
	int i;
	int **mtx;

	mtx = (int **)calloc( ll1+1, sizeof( int * ) );
	if( !mtx )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "Allocation error( %d x %d int mtx )\n", ll1, ll2 );
		exit( 1 );
	}
	if( ll2 )
	{
		for( i=0; i<ll1; i++ ) mtx[i] = AllocateIntVec( ll2 );
	}
	else
	{
		for( i=0; i<ll1; i++ ) mtx[i] = NULL;
	}
	mtx[ll1] = NULL;
	return( mtx );
}

/*
void FreeIntMtx( int **mtx )
{
*
	int **x;
	x = mtx;
	while( !*mtx ) free( *mtx++ );
	free( x );
*
	int i;
	for( i=0; mtx[i] != NULL; i++ ) 
		free( (char *)mtx[i] );
	free( (char *)mtx );
}
*/

char ***AllocateCharCub( int ll1, int ll2, int  ll3 )
{
	CALLS && printf("called %s:AllocateCharCub()\n", __FILE__);
	int i;
	char ***cub;

	cub = (char ***)calloc( ll1+1, sizeof( char ** ) );
	if( !cub ) 
	{
		fprintf( stderr, "Allocation error( %d x %d x %d char cube\n", ll1, ll2, ll3 );
		exit( 1 );
	}
	if( ll2 )
	{
		for( i=0; i<ll1; i++ ) 
		{
			cub[i] = AllocateCharMtx( ll2, ll3 );
		}
	}
	cub[ll1] = NULL;
	return( cub );
}

void FreeCharCub( char ***cub )
{
	CALLS && printf("called %s:FreeCharCub()\n", __FILE__);
	int i;

	for( i=0; cub[i]; i++ ) 
	{
		FreeCharMtx( cub[i] );
	}
	free( cub );
}

void freeintmtx( int **mtx, int ll1, int ll2 )
{
		CALLS && printf("called %s:freeintmtx()\n", __FILE__);
    int i;

    for( i=0; i<ll1; i++ ) 
        free( (char *)mtx[i] );
    free( (char *)mtx );
}
      
void FreeIntMtx( int **mtx ) // yes
{
	CALLS && printf("called %s:FreeIntMtx()\n", __FILE__);
	int i;

	for( i=0; mtx[i]; i++ ) 
	{
		if( mtx[i] ) free( (char *)mtx[i] ); mtx[i] = NULL;
	}
	free( (char *)mtx );
}

char ****AllocateCharHcu( int ll1, int ll2, int ll3, int ll4 )
{
	CALLS && printf("called %s:AllocateCharHcu()\n", __FILE__);
	int i;
	char ****hcu;

	hcu = (char ****)calloc( ll1+1, sizeof( char *** ) );
	if( hcu == NULL ) exit( 1 );
	for( i=0; i<ll1; i++ ) 
		hcu[i] = AllocateCharCub( ll2, ll3, ll4 );
	hcu[ll1] = NULL;
	return( hcu );
}

void FreeCharHcu( char ****hcu )
{
	CALLS && printf("called %s:FreeCharHcu()\n", __FILE__);
	int i;
	for( i=0; hcu[i]; i++ )
	{
		FreeCharCub( hcu[i] );
	}
	free ( (char *)hcu );
}

double *AllocateDoubleVec( int ll1 ) // yes
{
	CALLS && printf("called %s:AllocateDoubleVec()\n", __FILE__);
	double *vec;

	vec = (double *)calloc( ll1, sizeof( double ) ); // filled with 0.0
	return( vec );
}

void FreeDoubleVec( double *vec ) // yes
{
	CALLS && printf("called %s:FreeDoubleVec()\n", __FILE__);
	free( vec );
}

int ***AllocateIntCub( int ll1, int ll2, int ll3 ) // yes
{
	CALLS && printf("called %s:AllocateIntCub()\n", __FILE__);
	int i;
	int ***cub;

	cub = (int ***)calloc( ll1+1, sizeof( int ** ) );
	if( cub == NULL )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "cannot allocate IntCub\n" );
		exit( 1 );
	}
	for( i=0; i<ll1; i++ ) 
		cub[i] = AllocateIntMtx( ll2, ll3 );
	cub[ll1] = NULL;

	return cub;
}

void FreeIntCub( int ***cub )
{
	CALLS && printf("called %s:FreeIntCub()\n", __FILE__);
	int i;
	for( i=0; cub[i]; i++ ) 
	{
		if( cub[i] ) FreeIntMtx( cub[i] ); cub[i] = NULL;
	}
	free( cub );
}

double **AllocateDoubleMtx( int ll1, int ll2 ) // yes
{
	CALLS && printf("called %s:AllocateDoubleMtx()\n", __FILE__);
	int i;
	double **mtx;
	mtx = (double **)calloc( ll1+1, sizeof( double * ) );
	if( !mtx )
	{
		FILES && printf("file write stderr %s:%d\n", __FILE__, __LINE__);
		fprintf( stderr, "cannot allocate DoubleMtx\n" );
		exit( 1 );
	}
	if( ll2 )
	{
		for( i=0; i<ll1; i++ ) 
			mtx[i] = AllocateDoubleVec( ll2 );
	}
	mtx[ll1] = NULL;

	return mtx;
}

void FreeDoubleMtx( double **mtx ) // yes
{
	CALLS && printf("called %s:FreeDoubleMtx()\n", __FILE__);
	int i;
	for( i=0; mtx[i]; i++ )
		FreeDoubleVec( mtx[i] );
	free( mtx );
}

float ***AllocateFloatCub( int ll1, int ll2, int  ll3 )
{
	CALLS && printf("called %s:AllocateFloatCub()\n", __FILE__);
	int i;
	float ***cub;

	cub = (float ***)calloc( ll1+1, sizeof( float ** ) );
	if( !cub ) 
	{
		fprintf( stderr, "cannot allocate float cube.\n" );
		exit( 1 );
	}
	for( i=0; i<ll1; i++ ) 
	{
		cub[i] = AllocateFloatMtx( ll2, ll3 );
	}
	cub[ll1] = NULL;
	return( cub );
}

void FreeFloatCub( float ***cub )
{
	CALLS && printf("called %s:FreeFloatCub()\n", __FILE__);
	int i;

	for( i=0; cub[i]; i++ ) 
	{
		FreeFloatMtx( cub[i] );
	}
	free( cub );
}

double ***AllocateDoubleCub( int ll1, int ll2, int  ll3 )
{
	CALLS && printf("called %s:AllocateDoubleCub()\n", __FILE__);
	int i;
	double ***cub;

	cub = (double ***)calloc( ll1+1, sizeof( double ** ) );
	if( !cub ) 
	{
		fprintf( stderr, "cannot allocate double cube.\n" );
		exit( 1 );
	}
	for( i=0; i<ll1; i++ ) 
	{
		cub[i] = AllocateDoubleMtx( ll2, ll3 );
	}
	cub[ll1] = NULL;
	return( cub );
}

void FreeDoubleCub( double ***cub )
{
	CALLS && printf("called %s:FreeDoubleCub()\n", __FILE__);
	int i;

	for( i=0; cub[i]; i++ ) 
	{
		FreeDoubleMtx( cub[i] );
	}
	free( cub );
}


short *AllocateShortVec( int ll1 )
{
	CALLS && printf("called %s:AllocateShortVec()\n", __FILE__);
	short *vec;

	vec = (short *)calloc( ll1, sizeof( short ) );
	if( vec == NULL )
	{	
		fprintf( stderr, "Allocation error( %d short vec )\n", ll1 );
		exit( 1 );
	}
	return( vec );
}	

void FreeShortVec( short *vec )
{
	CALLS && printf("called %s:FreeShortVec()\n", __FILE__);
	free( (char *)vec );
}

short **AllocateShortMtx( int ll1, int ll2 )
{
	CALLS && printf("called %s:AllocateShortMtx()\n", __FILE__);
	int i;
	short **mtx;


	mtx = (short **)calloc( ll1+1, sizeof( short * ) );
	if( !mtx )
	{
		fprintf( stderr, "Allocation error( %d x %d short mtx ) \n", ll1, ll2 );
		exit( 1 );
	}
	for( i=0; i<ll1; i++ ) 
	{
		mtx[i] = AllocateShortVec( ll2 );
	}
	mtx[ll1] = NULL;
	return( mtx );
}

void FreeShortMtx( short **mtx )
{
	CALLS && printf("called %s:FreeShortMtx()\n", __FILE__);
	int i;

	for( i=0; mtx[i]; i++ ) 
		free( (char *)mtx[i] );
	free( (char *)mtx );
}

