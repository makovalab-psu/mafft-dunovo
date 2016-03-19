#include <stdio.h>
#include "dp.h"
#include "mltaln.h"

int TLS commonAlloc1 = 0;
int TLS commonAlloc2 = 0;
int TLS **commonIP = NULL;
int TLS **commonJP = NULL;
int nthread = 1;
int randomseed = 0;
int parallelizationstrategy = BAATARI1;


char modelname[500];
int njob, nlenmax;
int amino_n[0x80];
char amino_grp[0x80];
int amino_dis[0x80][0x80];
double **n_disLN;
double amino_dis_consweight_multi[0x80][0x80];
int **n_dis;
int **n_disFFT;
double **n_dis_consweight_multi;
char amino[0x80];
double polarity[0x80];
double volume[0x80];
int ribosumdis[37][37];

int ppid;
double thrinter;
double fastathreshold;
int pslocal, ppslocal;
int constraint;
int divpairscore;
int fmodel; // 1-> fmodel 0->default -1->raw
int nblosum; // 45, 50, 62, 80
int kobetsubunkatsu;
int bunkatsu;
int dorp;
int niter;
int contin;
int calledByXced;
int devide;
int scmtd;
int weight;
int utree;
int tbutree;
int refine;
int check;
double cut;
int cooling;
int trywarp = 0;
int penalty, ppenalty, penaltyLN;
int penalty_dist, ppenalty_dist;
int RNApenalty, RNAppenalty;
int RNApenalty_ex, RNAppenalty_ex;
int penalty_ex, ppenalty_ex, penalty_exLN;
int penalty_EX, ppenalty_EX;
int penalty_OP, ppenalty_OP;
int penalty_shift, ppenalty_shift;
double penalty_shift_factor = 100.0;
int RNAthr, RNApthr;
int offset, poffset, offsetLN, offsetFFT;
int scoremtx;
int TMorJTT;
char use_fft;
char force_fft;
int nevermemsave;
int fftscore;
int fftWinSize;
int fftThreshold;
int fftRepeatStop;
int fftNoAnchStop;
int divWinSize;
int divThreshold;
int disp;
int outgap = 1;
char alg;
int cnst;
int mix;
int tbitr;
int tbweight;
int tbrweight;
int disopt;
int pamN;
int checkC;
float geta2;
int treemethod;
int kimuraR;
char *swopt;
int fftkeika;
int score_check;
int makedistmtx;
char *inputfile;
char *addfile;
int addprofile = 1;
int rnakozo;
char rnaprediction;
int scoreout = 0;
int spscoreout = 0;
int outnumber = 0;
int legacygapcost = 0;

char *signalSM;
FILE *prep_g;
FILE *trap_g;
char **seq_g;
char **res_g;

float consweight_multi = 1.0;
float consweight_rna = 0.0;
char RNAscoremtx = 'n';

char TLS *newgapstr = "-";

int nalphabets = 26;
int nscoredalphabets = 20;

double specificityconsideration = 0.0;
int ndistclass = 10;
int maxdistclass = -1;

int gmsg = 0;

double sueff_global = SUEFF;

void initglobalvariables()
{
	commonAlloc1 = 0;
	commonAlloc2 = 0;
	commonIP = NULL;
	commonJP = NULL;
	nthread = 1;
	randomseed = 0;
	parallelizationstrategy = BAATARI1;

	trywarp = 0;
	penalty_shift_factor = 100.0;
	outgap = 1;
	addprofile = 1;
	scoreout = 0;
	outnumber = 0;
	legacygapcost = 0;
	consweight_multi = 1.0;
	consweight_rna = 0.0;
	RNAscoremtx = 'n';

	newgapstr = "-";

	nalphabets = 26;
	nscoredalphabets = 20;

	specificityconsideration = 0.0;
	ndistclass = 10;
	maxdistclass = -1;
	
	gmsg = 0;
}
