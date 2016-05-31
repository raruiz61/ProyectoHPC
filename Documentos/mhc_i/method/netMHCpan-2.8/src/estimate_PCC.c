/* M.Nielsen April 2002 mniel@cbs.dtu.dk */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

FILENAME        p_matrix;
float		p_alpha;
float		p_beta;
FILENAME	p_pseudolist;

PARAM   param[] = {
	"-a", VFLOAT	p_alpha, "Alpha for fit", "-1.1612",
	"-b", VFLOAT	p_beta, "Beat for fit", "0.8530",
	"-m", VFNAME    p_matrix, "Alignment matrix", "$NETMHCpan/data/matrices/BLOSUM62",
	"-p", VFNAME	p_pseudolist, "List of reference pseudo sequences", "$NETMHCpan/data/training.pseudo",
	0
};

main( int argc, char *argv[] )

{

	PAIRLIST	*ref_list, *neighbor, *ref;
	char	*pseudoseq;
	float	mindist;
	float	pcc;
	float	sco, ssco1, ssco2, dist;
	int	len;
	int	ia, ib, i;

	pparse( &argc, &argv, param, 1, "pseudoseq" );

	ref_list = pairlist_read( p_pseudolist );

	if ( ! ref_list ) {
		printf( "Error. Cannot read REF_LIST from file %s. Exit\n", p_pseudolist );
		exit( 1 );
	}

	bl_init_file( p_matrix );

	pseudoseq = argv[1];

	if ( strlen( ref_list->name2) !=  strlen( pseudoseq ) ) {
		printf( "Error. Inconsistent pseudo sequence length %s %s\n", ref_list->name2, pseudoseq );
		exit( 1 );
	}

	mindist = 2.0;
	len = strlen( pseudoseq );

	for ( ref= ref_list; ref; ref=ref->next ) {

		sco = 0.0;
		ssco1 = 0.0;
		ssco2 = 0.0;

		for ( i=0;i<len; i++ ) {

			ia = strpos(PROFILE_ORDER, ref->name2[i] );
			ib = strpos(PROFILE_ORDER, pseudoseq[i] );

			if ( ia < 0 || ia > 19 ) {
				printf( "Error. Wrong pseudo seq character %c %s. Exit\n", 
					ref->name2[i], PROFILE_ORDER );
				exit( 1 );
			}

			if ( ib < 0 || ib > 19 ) {
				printf( "Error. Wrong pseudo seq character %c %s. Exit\n", 
					pseudoseq[i], PROFILE_ORDER );
				exit( 1 );
			}

			sco += bl_s( ia, ib );
			ssco1 += bl_s( ia, ia );
			ssco2 += bl_s( ib, ib );

		}

		dist = 1.0 - sco/( sqrt( ssco1 * ssco2 ) );

		if ( dist < mindist ) {
			
			mindist = dist;
			neighbor = ref;
		}
	}

	pcc = p_alpha * mindist + p_beta;

	printf( "%6.4f %s %6.4f\n", mindist, neighbor->name1, pcc );

	exit( 0 );
}
