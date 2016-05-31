/* M.Nielsen August 2002 mniel@cbs.dtu.dk */

#include <stdio.h>
#include <math.h>
#include "utils.h"
#include <sys/types.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

FILENAME	p_spsyn;
FILENAME        p_blsyn;
FILENAME	p_seq2inp;
int		p_verbose;
int		p_writeseq;
int		p_dirty;
int		p_block;
int		p_sdev;
FILENAME	p_tmpdir;
int		p_printall;
FILENAME        p_bdir;

PARAM   param[] = {
	"-spsyn", VFNAME p_spsyn, "Sparse synaps list", 
		"/home/people/mniel/epipred/data/useall_hmm+seq.4.synlist",
	"-blsyn", VFNAME p_blsyn, "Blosum synaps list",
		"/home/people/mniel/epipred/data/useall_hmm+bl50.3.synlist",
	"-seq2inp", VFNAME p_seq2inp, "Sequence to NN input program", "seq2inp",
	"-writeseq", VSWITCH p_writeseq, "Write sequence in output", "1",
	"-v", VSWITCH p_verbose, "Verbose mode", "0",
	"-dirty", VSWITCH p_dirty, "Dirty mode, leave tmp dir+files", "0",
	"-block", VINT p_block, "Block size of predictions", "50000",
	"-sdev", VSWITCH p_sdev, "Print sdev", "0",
	"-tdir", VFNAME p_tmpdir, "Temporary directory (Default $$)", "",
	"-bdir",  VFNAME p_bdir, "Binary directory", "",
	"-pa", VSWITCH	p_printall, "Print all network predictions", "0",
	0
};

void	doit( FILENAME cmd )

{
	if ( p_verbose )
		printf( "# %s", cmd );

	system( cmd );
}

void	epipred( FILENAME filename, FILENAME dir, int noff )

{
	int		n;
	FILENAME	cmd;
	float		bl, seq, o;
	float		sseq, sbl;
	LINELIST	*linelist_pep, *ll, *ll_pep;
	LINELIST	*linelist_bl, *ll_bl;
	LINELIST        *linelist_sp, *ll_sp;
	LINE		seq_all, bl_all;
	WORD		pep;
	LINE		restofline;
	float		sdev;
	int		nseq, nbl, npep;
	int		i,len;
	int		offset;

#ifdef NOTNOW
	sprintf( cmd, "%s%s -min 0 -max 1 %s | grep -v \"#\" > %s/seq.inp\n", 
		p_bdir, p_seq2inp, filename, dir );
	
	doit( cmd );
#endif


#ifdef NOTNOW
	sprintf( cmd, "%snnlinplayer -x %s %s %s %s/seq.inp | grep -v \"#\" > %s/bl.res\n", 
#endif
	sprintf( cmd, "%snnlinplayer -x -seqinp %s %s %s %s | grep -v \"#\" > %s/bl.res\n",
		p_bdir,
		( p_sdev ? "-sdev 1 " : " " ), 
		( p_printall ? "-a " : " " ), 
		p_blsyn, filename, 
		dir );

	doit( cmd );

#ifdef NOTNOW
	sprintf( cmd, "%snnlinplayer -x %s %s %s %s/seq.inp | grep -v \"#\" > %s/sp.res\n", 
#endif
	sprintf( cmd, "%snnlinplayer -x -seqinp %s %s %s %s | grep -v \"#\" > %s/sp.res\n",
		p_bdir,
		( p_sdev ? "-sdev 1 " : " " ),
		( p_printall ? "-a " : " " ),
		p_spsyn, filename,  
		dir );

	doit( cmd );

	sprintf( cmd, "cat %s/sp.res", dir );

	linelist_sp = linelist_pread( cmd );

	sprintf( cmd, "cat %s/bl.res", dir );

	linelist_bl = linelist_pread( cmd );

	if ( p_writeseq ) 
		linelist_pep = linelist_read( filename );
	else
		linelist_pep = NULL;

	n=0;

	for ( ll=linelist_sp, nseq = 0; ll; ll=ll->next, nseq++ );
	for ( ll=linelist_bl, nbl = 0; ll; ll=ll->next, nbl++ );
	for ( ll=linelist_pep, npep = 0; ll; ll=ll->next, npep++ );

	if ( nseq != nbl ) {
		printf( "Error. nseq != nbl %i %i\n", nseq, nbl );
		exit( 1 );
	}

	if ( p_writeseq && nseq != npep ) {
		printf( "Error. nseq != npep %i %i\n", nseq, npep );
		exit( 1 );
	}

	if ( p_writeseq && nbl != npep ) {
		printf( "Error. nbl != npep %i %i\n", nbl, npep );
		exit( 1 );
	}

	for ( ll_sp = linelist_sp, ll_bl = linelist_bl, ll_pep=linelist_pep; 
		ll_sp && ll_bl;
		ll_sp = ll_sp->next, ll_bl = ll_bl->next ) {

		if ( p_printall ) {
			i = 0;
			len = strlen( ll_sp->line );
			while ( i<len && strncmp( ll_sp->line+i, "Target=", 7 ) != 0 ) 
				i++;

			i-=10;

			strncpy( seq_all, ll_sp->line, i );
			seq_all[i] = 0;

			offset = i;

		}
		else {
			offset = 0;
		}

		if ( p_sdev )
			sscanf( ll_sp->line+offset, "%f %*s %*f %*s %f", &seq, &sseq );
		else
			sscanf( ll_sp->line+offset, "%f", &seq );

		if ( p_printall ) {
			i = 0;
			len = strlen( ll_bl->line );
			while ( i<len && strncmp( ll_bl->line+i, "Target=", 7 ) != 0 )
				i++;   

			i-=10;
			
			strncpy( bl_all, ll_bl->line, i );
			bl_all[i] = 0;
			
			offset = i;
		}
		else {
			offset = 0;
		}

		if ( p_sdev )
			sscanf( ll_bl->line+offset, "%f %*s %*f %*s %f", &bl, &sbl );
		else
			sscanf( ll_bl->line+offset, "%f", &bl );

		o = 0.5*(seq + bl);

		printf( "%5i ", n+noff );	

		if ( p_writeseq ) {

			if ( !ll_pep ) {
				printf( "Error. Linelist empty\n" );
				exit( 1 );
			}

			if ( sscanf( ll_pep->line, "%s %[^\n]", pep, restofline ) == 1 ) {
				strcpy( restofline, "0.0" );
			}

			printf( "%14s %8s ", pep, restofline );

			ll_pep = ll_pep->next;
		}

		printf( "%8.6f", o );

		if ( p_sdev ) {
			sdev = sqrt( 0.5 * ( sseq * sseq + sbl * sbl ) );
			printf( " %8.6f", sdev );
		}

		if ( p_printall ) {
			printf( " Sparse %s", seq_all );
			printf( " BL %s", bl_all );
		}

		printf( "\n" );

		n++;

	}

	if ( linelist_pep )
		linelist_free( linelist_pep );
	
	if ( linelist_sp )
		linelist_free( linelist_sp );

	if ( linelist_bl )
		linelist_free( linelist_bl );

	printf( "# %i prediction done from file %s\n", n, filename );
		
}

main( int argc, char *argv[] )

{
	int		pid;
	FILENAME	filename, cmd;
	LINE		line;
	FILE		*fp, *fp1;
	int		fc, ff;
	int		nlines, n, np;
	FILENAME		tmpdir;

	pparse( &argc, &argv, param, 1, "filename" );

	if ( strlen( p_tmpdir ) == 0 ) {
		pid = getpid();
		sprintf( tmpdir, "%i", pid );
	}
	else {
		sprintf( tmpdir, "%s", p_tmpdir );
	}

	sprintf( cmd, "mkdir -p %s\n", tmpdir);

	doit( cmd );

        if ( ( fp = stream_input( argv[1], &fc, &ff ) ) == NULL ) {
                printf( "Error. Cannot open file %s\n", argv[1] );
                exit( 1 );
        }

	n = 0;
	np = 0;
	nlines = 0;

	sprintf( filename, "%s/%i.dat", tmpdir, np );

	if ( ( fp1 = fopen( filename, "w" ) ) == NULL ) {
		printf( "Cannot open file %s\n", filename );
		exit( 1 );
	}

	printf( "%5s %14s %8s %8s\n", "#   n", "pep", "Aff", "NN_pred" );

        while ( fgets( line, sizeof line, fp ) != NULL ) {

                if ( strlen( line ) <= 1 ) continue;

                if ( strncmp( line, "#", 1 ) == 0 ) continue;

		fprintf( fp1, "%s", line );

		n++;

		if ( n >= p_block ) {

			n=0;
			np++;

			fclose( fp1 );

			epipred( filename, tmpdir, (np-1)*p_block );

			sprintf( filename, "%s/%i.dat", tmpdir, np );

			if ( ( fp1 = fopen( filename, "w" ) ) == NULL ) {
				printf( "Cannot open file %s\n", filename );
				exit( 1 );
			}
		}

		nlines++;
	}

	if ( n > 0 ) {
		n=0;
		np++;

		fclose( fp1 );

		epipred( filename, tmpdir, (np-1)*p_block );
	}

	if ( ! p_dirty ) {
		sprintf( cmd, "rm -rf %s\n", tmpdir );
		doit( cmd );
	}

	exit( 0 );
}
