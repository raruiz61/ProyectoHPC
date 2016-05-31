/* M.Nielsen April 2002 mniel@cbs.dtu.dk */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"

static double squash( a ) double a; { return (1.0 / (1.0 + exp( -a ))); }

int	p_ptarget;
int	p_sdev;
int	p_synfile;
int	p_verbose;
FILENAME	p_use;
int	p_printall;
FILENAME	p_prefix;
int	p_seqinp;
WORD	p_alphabet;
int	p_x;

PARAM   param[] = {
	"-pt", 	VINT	p_ptarget,	"Print target value", "1",
	"-sdev", VINT	p_sdev,		"Print sdev",	"0",
	"-sfil", VINT	p_synfile,	"[0] Use synlist, [1] Use single synfile",	"0",
	"-v", VSWITCH	p_verbose,	"Verbose mode", "0",
	"-s", VSWITCH	p_synfile, 	"Use single synfile [default is synlist]", "0",
	"-sd", VSWITCH  p_sdev,		"Include standart deviation in output", "0",
	"-u", VFNAME	p_use,		"List of networks to use in prediction (default all)", "",
	"-a", VSWITCH	p_printall, 	"Print prediction of all networks", "0",
	"-px", VFNAME	p_prefix, 	"Prefix for synaps files", "",
	"-seqinp", VSWITCH p_seqinp, 	"Use sequence input (only for rescaled ANNs)", "0",
	"-alpha", VWORD	p_alphabet, 	"Alphabet for sequence encoding (only for rescaled ANNs)", "ARNDCQEGHILKMFPSTWYV",
	"-x", VSWITCH	p_x, 		"Include extra neuron to handle X's", "0",
	0
};

/* GLOBAL variables */

float	*inp;

int	scan_input_fp( FILE *fp, int n )

{
	int 	i;
	char	ch;
	LINE	line;

	while ( ( ch = fgetc( fp ) ) && ungetc( ch, fp ) && ch == '#' ) { /* comment line */
                fgets( line, sizeof(line), fp );
                if ( p_verbose )
                        printf( "# %s", line );
        }

	for ( i = 0; i <= n; i++)
        	if ( fscanf(fp, "%f", &inp[i]) != 1 ) 
			return( 0 );

	return( 1 );
}

void    setnet(float **x, SYNAPS * syn)

{
        int     l;
        for (l = 0; l < syn->nlay; l++)
                x[l][syn->nnlay[l]] = 1.0;
}

main( int argc, char *argv[] )

{

	SYNAPS	*syn, *s;
	FILE    *fp;
	int	i, l, k, j;
	float	**x;
	float	a, *p1, *p2, *sdev;
	int	nlines, ns, nout, nin, n, nnu, nu;
	int	fc, ff;
	float	*use;
	int	synnlay, maxnper;
	LONGLINE	tmpfilename;
	int	*ix=NULL;
	LINE            line, restofline;
        LINE            seq;
	int	len, len0, alen;
	float	**wgt_l,*wgt_lk;
	float	*x_l, *x_l1;
	int	c;

	pparse( &argc, &argv, param, 2, "synapslist/synapsfile inputfile" );

	if ( strlen( p_use ) > 0 ) {
		use = fvector_split( p_use, &nu );
	}
	else
		use = NULL;

	if ( p_synfile ) {
		syn = synaps_read_nfile(argv[1], &n);
		printf( "# Number of synaps file is %s %i\n", argv[1], n );
	}
	else
		syn = synaps_read_pfix(argv[1], p_prefix);

	for ( i=0; i<syn->nlay; i++ )
		printf( "# Nlayer %i NNlay %i\n", i, syn->nnlay[i] );

	for ( synnlay=-99, s=syn; s; s=s->next )
			synnlay = ( s->nlay > synnlay ? s->nlay : synnlay );

	for ( maxnper=-99, s=syn; s; s=s->next )
		for ( i=0; i<s->nlay; i++ )
			if ( s->nnlay[i] > maxnper )
				maxnper = s->nnlay[i];

        if ( ( fp = stream_input( argv[2], &fc, &ff ) ) == NULL ) {
                printf( "Error Cannot open file %s for input. Exit\n", argv[2] );
                exit( 1 );
        }

	nout = syn->nnlay[syn->nlay-1];
	nin = syn->nnlay[0];

	x = fmatrix( 0, synnlay-1, 0, maxnper );

	inp = fvector( 0, nin+nout-1 );
	p1 = fvector( 0, nout-1 );
	p2 = fvector( 0, nout-1 );
	sdev = fvector( 0, nout-1 );

	printf( "# Ninput %i. Noutput %i\n", nin, nout );

	nlines = 0;

	if ( p_seqinp ) {

		printf( "# Using sequence input\n" );
		alen = strlen( p_alphabet );

		if ( p_x )
			alen++;

		while ( fgets( line, sizeof line, fp ) != NULL ) {

                	if ( strlen( line ) <= 1 ) continue;

                	if ( strncmp( line, "#", 1 ) == 0 ) continue;

			restofline[0] = 0;
                	if ( sscanf( line, "%s %[^\n]", seq, restofline ) < 1 ) {
                        	printf( "Error reading line %s. Exit\n", line );
                        	exit( 1 );
                	}

			len = strlen(seq);

			if ( ix == NULL ) {
				ix = ivector( 0, len-1 );
				len0 = len;
			}

			if ( len != len0 ) {
				printf( "Inconsistent input length %s %i %i. Exit\n", seq, len0, len );
				exit( 1 );
			}

			for ( i=0; i<len; i++ ) {

				if ( p_x && seq[i] == 'X' )
					ix[i] = 20;
				else {
					c = strpos( p_alphabet, seq[i] );
					ix[i] = c;
				}
			}

			for ( i=0; i<nout; i++ ) {
				p1[i] = 0.0;
				p2[i] = 0.0;
				sdev[i] = 0.0;
			}

			for ( ns=0, nnu=0,s=syn; s; s=s->next, nnu++ ) {

				if ( use && nnu >= nu ) {
					printf( "Error. Usage list %s is longer %i than number of Networks\n", p_use, nu );
					exit( 1 );
				}

				if ( use && ! use[nnu] ) continue;

				setnet( x, s );

				x_l = x[1];
				wgt_l = s->wgt[0];

				for ( k = 0; k < s->nnlay[1]; k++ ) {
					a = 0.0;
					wgt_lk = wgt_l[k];
					for (i = 0; i<len; i++ )
						if ( ix[i] >= 0 )
							a += wgt_lk[i*alen+ix[i]];
					a += wgt_lk[s->nnlay[0]];
					x_l[k] = squash( a );
#ifdef DEBUG
					printf( "TEST %i %i %f %f\n", 1, k, a, x[1][k] );
#endif
				}

				for (  l = 2; l < s->nlay; l++ ) {
					x_l = x[l];
					wgt_l = s->wgt[l-1];
					x_l1 = x[l-1];
					for ( k = 0; k < s->nnlay[l]; k++ ) {
						a = 0.0;
						wgt_lk = wgt_l[k];
						for (i = 0; i < s->nnlay[l-1]+1; i++ ) {
							a += x_l1[i] * wgt_lk[i];
						}
						x_l[k] = squash( a );
					}
				}

				if ( p_printall )
					printf( "%3i ", ns+1 );

				for ( i=0; i<nout; i++ ) {

					if ( p_printall ) 
						printf( "%6.4f ", x[s->nlay-1][i] );

					p1[i] += x[s->nlay-1][i];
					p2[i] += x[s->nlay-1][i] * x[s->nlay-1][i];
				}

				ns++;

			}

			for ( i=0; i<nout; i++ ) {
				p1[i] = ( ns> 0 ? p1[i]/ns : 0.0 );
				p2[i] = ( ns> 0 ? p2[i]/ns : 0.0 );
				sdev[i] = p2[i] - p1[i]*p1[i];

				sdev[i] = ( sdev[i] > 0.0 ? sqrt( sdev[i] ) : 0.0 );

				printf( "%8.6f ", p1[i] );

				if ( p_ptarget ) 
					printf( "Target= %8s ", restofline  );

				if ( p_sdev )
					printf( "Sdev= %8.6f ", sdev[i] );

			}

			printf( "\n" );

			nlines++;
		}

	}
	else {
		while ( ! p_seqinp && scan_input_fp( fp, nin+nout-1 ) ) {

			for ( i=0; i<nout; i++ ) {
				p1[i] = 0.0;
				p2[i] = 0.0;
				sdev[i] = 0.0;
			}

			for ( ns=0, nnu=0,s=syn; s; s=s->next, nnu++ ) {

				if ( use && nnu >= nu ) {
					printf( "Error. Usage list %s is longer %i than number of Networks\n", p_use, nu );
					exit( 1 );
				}

				if ( use && ! use[nnu] ) continue;

				setnet( x, s );

				for ( i=0; i<s->nnlay[0]; i++ ) 
					x[0][i] = inp[i];

				for (  l = 1; l < s->nlay; l++ ) {
					for ( k = 0; k < s->nnlay[l]; k++ ) {
						a = 0.0;
						for (i = 0; i < s->nnlay[l-1]+1; i++ ) {
							a+= x[l-1][i] * s->wgt[l-1][k][i];
						}
						x[l][k] = squash( a );
					}
				}

				if ( p_printall )
					printf( "%3i ", ns+1 );

				for ( i=0; i<nout; i++ ) {

					if ( p_printall ) 
						printf( "%6.4f ", x[s->nlay-1][i] );

					p1[i] += x[s->nlay-1][i];
					p2[i] += x[s->nlay-1][i] * x[s->nlay-1][i];
				}

				ns++;

			}

			for ( i=0; i<nout; i++ ) {
				p1[i] = ( ns> 0 ? p1[i]/ns : 0.0 );
				p2[i] = ( ns> 0 ? p2[i]/ns : 0.0 );
				sdev[i] = p2[i] - p1[i]*p1[i];

				sdev[i] = ( sdev[i] > 0.0 ? sqrt( sdev[i] ) : 0.0 );

				printf( "%8.6f ", p1[i] );

				if ( p_ptarget ) 
					printf( "Target= %8.6f ", inp[syn->nnlay[0]+i] );

				if ( p_sdev )
					printf( "Sdev= %8.6f ", sdev[i] );

			}

			printf( "\n" );

			nlines++;
		}
	}

        stream_close( fp, fc, argv[2] );

	printf( "# Nlines %i processed thorugh %i Neural Networks\n", nlines, ns );

	exit( 1 );
}
