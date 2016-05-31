/* M.Nielsen June 2008 mniel@cbs.dtu.dk */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

float   p_gapf;
float   p_gapn;
int     p_minall;
FILENAME	p_refpos;
FILENAME	p_refseq;
FILENAME        p_matrix;
int	p_newnumber;
int	p_verbose;

PARAM   param[] =
{
	"-gf", VFLOAT p_gapf, "penalty for first gap", "11",
	"-gn", VFLOAT p_gapn, "penalty for next gap", "1",
	"-mal", VINT p_minall, "Minimum alignment lenght", "3",
	"-p", VFNAME p_refpos, "Pseudo positions file name", "$NETMHCpan/data/all_varcontacts.nlist",
	"-r", VFNAME p_refseq, "Align reference file name", "$NETMHCpan/data/B0702.fsa",
	"-m", VFNAME    p_matrix, "Alignment matrix", "$NETMHCpan/data/matrices/BLOSUM50",
	"-n", VSWITCH p_newnumber, "Use new numbers", "0",
	"-v", VSWITCH p_verbose, "Verbose mode", "0",
	0
};

ALN    *agap_align(float **m, PRF * q, PRF * b, float self_score, float gapf, float gapn)

{
	float **sco;
	float   score;
	int     firsti, firstj;
	int     alen;
	ALN    *new = NULL;
	int     i;
	char	qpal[3000], dpal[3000];
	float	s1, s2;

	sco = fmatrix(0, q->chlen - 1, 0, b->chlen - 1);

#ifdef NOTNOW
        nw_sco2(m, q->chlen, b->chlen, gapf, gapn, 0, &score, sco );
        nw_ali2(q->chlen, b->chlen, gapf, gapn, 0, sco, q->seq, b->seq,
		qpal, dpal, &alen );

        firsti = 0;
        firstj = 0;

#else
	sw_sco(m, q->chlen, b->chlen, gapf, gapn, &score, sco, &firsti, &firstj);
	sw_ali(q->chlen, b->chlen, gapf, gapn, sco, q->seq, b->seq,
                        qpal, dpal, &alen, firsti, firstj);
#endif		

	fmatrix_free(sco, 0, q->chlen - 1, 0, b->chlen - 1);

	if (alen < p_minall)
		return( NULL );

	if ((new = aln_alloc()) == NULL) {
		printf("Cannot alloc new ALN\n");
		exit(1);
	}

	new->alen = alen;

	new->score = score;

	new->mlen = 0;
	new->ngap = 0;
	new->nid = 0;

	new->sscore = self_score;

	for (i = 0; i < new->alen; i++) {
		if (qpal[i] != '-' && dpal[i] != '-') {
			new->mlen++;
		} else {
			new->ngap++;
		}
		if (qpal[i] == dpal[i])
			new->nid++;
	}

	new->qof = firsti;
	new->qal = cvector(0, strlen(qpal));
	strcpy(new->qal, qpal);

	new->dof = firstj;
	new->dal = cvector(0, strlen(dpal));
	strcpy(new->dal, dpal);

	strcpy(new->qname, q->name);
	new->qlen = q->len;
	new->qchlen = q->chlen;

	new->qseq = (q->seq ? cvector(0, strlen(q->seq)) : NULL);
	strcpy(new->qseq, q->seq);

	strcpy(new->dname, b->name);
	new->dlen = b->len;
	new->dchlen = b->chlen;

	new->dseq = (b->seq ? cvector(0, strlen(b->seq)) : NULL);
	strcpy(new->dseq, b->seq);

	strcpy(new->type, "prol_list");

	new->rscore = -new->score;

	return (new);

}

ALN    *align(PRF *d, PRF *q, float gapf, float gapn)

{

	float 	**m;
	float   self_score, ss_sec, ss_seq;
	ALN    	*new;
	int     i,np;

	m = sm_b62(q, d);

#ifdef NOTNOW
	self_score = ss_b62( q );
#endif

	new = agap_align(m, q, d, self_score, gapf, gapn);

	fmatrix_free(m, 0, q->chlen-1, 0, d->chlen-1 );

	return( new );

}

main(int argc, char *argv[])

{
	PRF     	*q, *d;
	FSALIST		*fsa, *fsalist;
	float   	gapf, gapn;
	ALN		*aln;
	int		*pos, np, i;
	LINELIST	*linelist;
	int		ix;

	pparse(&argc, &argv, param, 1, "fsa");

	bl_init_file( p_matrix );

	if ( ( fsa = fsalist_read( p_refseq ) ) == NULL ) {
                printf("Cannot read fasta file %s\n", p_refseq );
                exit(1);
        }
         
        if ( ( d = prf_alloc() ) == NULL ) {
                printf( "Cannot allocate Q PRF\n" );
                exit( 1 );
        }
         
        d->seq = fsa->seq;
        d->chlen = fsa->len;
        strncpy( d->name, fsa->name, 6 );
        d->name[6] = 0;
        
        prf_iassign_profile_order( d );

	gapf = p_gapf;
	gapn = p_gapn;

	if ( ( fsalist = fsalist_read( argv[1] ) ) == NULL ) {
		printf("Cannot read fasta file %s\n", argv[1] );
		exit(1);
	}

	fsalist_checkalpha( fsalist );

	if ( ( q = prf_alloc() ) == NULL ) {
		printf( "Cannot allocate Q PRF\n" );
		exit( 1 );
	}

	linelist = linelist_read( p_refpos );

	if ( ! linelist ) {
		printf( "Error. Cannot read reference position list from file %s\n", p_refpos );
		exit( 1 );
	}

	pos = ivector_split( linelist->line, &np );

	for ( fsa=fsalist; fsa; fsa = fsa->next ) {

		q->seq = fsa->seq;
        	q->chlen = fsa->len;
        	strncpy( q->name, fsa->name, 6 );
        	q->name[6] = 0;

        	prf_iassign_profile_order( q );

		aln = align( d, q, gapf, gapn);

		if ( ! aln ) {
			printf( "Error. No alignment was made. Exit\n" );
			exit( 1 );
		}

		if ( p_verbose )
			aln_write_single( aln );

		ix = 0;
		for ( i=0;i<aln->alen; i++ ) {
			if ( aln->dal[i] != '-' ) {
				aln->dal[ix] = aln->dal[i];
				aln->qal[ix] = aln->qal[i];
				ix++;
			}
		}

		aln->alen = ix;

#ifdef NOTNOW
		printf( "%s ", fsa->name );
#endif

		for ( i=0;i<np;i++ ) {
			if ( p_newnumber ) 
				ix = pos[i] - 1 - aln->dof;
			else
				ix = pos[i] - aln->dof + 23;

			if ( ix < 0 || ix > aln->alen-1 )
				printf( "%c", 'X' );
			else
				printf( "%c", aln->qal[ix] );
		}

		printf( "\n" );

		aln_free( aln );
	}

	exit(0);
}
