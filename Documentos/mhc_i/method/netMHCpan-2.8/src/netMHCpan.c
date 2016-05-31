/* M.Nielsen August 2002 mniel@cbs.dtu.dk */

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"

FILENAME	p_rdir;
FILENAME	p_syn;
int		p_verbose;
int		p_dirty;
int		p_sdev;
FILENAME	p_tmpdir;
FILENAME	p_hlapseudo;
FILENAME	p_hlaseq;
LINE		p_hla;
int		p_w;
FILENAME	p_f;
int		p_sort;
int		p_pep;
float		p_hthr;
float		p_lthr;
float           p_hthr_r;
float           p_lthr_r;
WORD		p_len;
int             p_xlsdump;
FILENAME        p_xlsfile;
float		p_thr;
int		p_list;
FILENAME	p_thrfmt;
int		p_ic50;
FILENAME        p_whitelist;
int		p_exprefix;
FILENAME	p_version;
int		p_inptype;

PARAM   param[] = {
	"-rdir",  VFNAME p_rdir, "Home directory for NetMHpan", "$NETMHCpan",
	"-syn", VFNAME p_syn, "Synaps file", 
		"$NETMHCpan/data/syn/synaps",
	"-v", VSWITCH p_verbose, "Verbose mode", "0",
	"-dirty", VSWITCH p_dirty, "Dirty mode, leave tmp dir+files", "0",
	"-sdev", VSWITCH p_sdev, "Print sdev", "0",
	"-tdir", VFNAME p_tmpdir, "Temporary directory (Default $$)", "$TMPDIR",
	"-hlapseudo", VFNAME p_hlapseudo, "File with HLA pseudo sequences", "$NETMHCpan/data/MHC_pseudo.dat",
	"-hlaseq", VFNAME p_hlaseq, "File with full length HLA sequences", "",
	"-a", VLINE	p_hla, "HLA allele", "HLA-A02:01",
	"-f", VFNAME    p_f, "File name with input", "",
        "-w", VSWITCH   p_w, "w option for webface", "0",
	"-s", VSWITCH	p_sort, "Sort output on descending affinity", "0",
	"-p", VSWITCH	p_pep, "Use peptide input", "0",
	"-th", VFLOAT	p_hthr, "Threshold for high binding peptides", "50.0",
	"-lt", VFLOAT   p_lthr, "Threshold for low binding peptides", "500.0",
	"-rth", VFLOAT	p_hthr_r, "Rank Threshold for high binding peptides", "0.5",
	"-rlt", VFLOAT   p_lthr_r, "Rank Threshold for low binding peptides", "2.0",
	"-l",  VWORD	p_len,	"Peptide length [8-11] (multiple length with ,)", "9",
        "-xls", VSWITCH p_xlsdump, "Save output to xls file", "0",
        "-xlsfile", VFNAME p_xlsfile, "Filename for xls dump", "NetMHCpan_out.xls",
	"-t", VFLOAT	p_thr, "Threshold for output", "-99.9",
	"-thrfmt", VFNAME p_thrfmt, "Format for threshold filenames", "$NETMHCpan/data/threshold/%s.thr",
	"-ic50", VSWITCH p_ic50, "Print IC50 values for all alleles", "0",
	"-whitelist", VFNAME p_whitelist, "List of alleles to include IC50 values", "$NETMHCpan/data/threshold/whitelist",
	"-expfix", VSWITCH p_exprefix, "Exclude prefix from synlist", "0",
	"-version", VFNAME p_version, "File with version information", "$NETMHCpan/data/version",
	"-inptype", VINT p_inptype, "Input type [0] FASTA [1] Peptide", "0",
	"-listMHC", VSWITCH p_list, "Print list of alleles included in netMHCpan", "0",
	0
};

void	doit( char *cmd )

{
	if ( p_verbose )
		printf( "# %s", cmd );

	system( cmd );
}

typedef struct thrlist {
        struct thrlist *next;
        WORD    all;
	float	thr;
	float	score;
}THRLIST;

THRLIST        *thrlist_alloc()

{
        THRLIST        *n;

        if ( ( n = ( THRLIST * ) malloc ( sizeof(THRLIST))) != NULL ) {
                strcpy( n->all, "UNDEF" );
		n->thr = -99.9;
		n->score = -99.9;
                n->next = NULL;
        }

        return( n );
}

THRLIST        *thrlist_read( char *filename )

{
        THRLIST        *first, *last, *new;
        FILE            *fp;
        int             ff, fc;
        WORD            all;
        LINE            line;
        int             n;
	float		thr, score;
	WORD		pep;

        first = NULL;
        n = 0;

        if ( ( fp = stream_input( filename, &fc, &ff )) == NULL ) {
                printf( "Error. Cannot read THRLIST from file %s. Exit\n",
                        filename );
                exit( 1 );
        }

        while ( fgets(line, sizeof line, fp) != NULL ) {

                if ( strncmp( line, "#", 1 ) == 0 )
                        continue;

                if ( strlen( line ) < 1 )
                        continue;

                if ( sscanf( line, "%s %f %s %f", all, &thr, pep, &score ) != 4 )
			continue;

                if ( ( new = thrlist_alloc()) == NULL ) {
                        printf( "Error. Cannot allocate THRLIST. Exit\n" );
                        exit( 1 );
                }

                strcpy( new->all, all );
		new->thr = thr;
		new->score = score;

                if ( first == NULL )
                        first = new;
                else
                        last->next = new;

                last = new;
                n++;
        }

        stream_close( fp, fc, filename );

	return( first );

}

void	thrlist_free( THRLIST *thrlist )

{
	THRLIST	*l, *n;

	for ( l=thrlist; l; l=n ) {
                n=l->next;
                free( l );
        }
} 

THRLIST	*find_thr( THRLIST *thrlist, float score )

{
	THRLIST *thrl;

	for ( thrl=thrlist; thrl; thrl=thrl->next )
		if ( score > thrl->score )
			return thrl;

	return( NULL );
}

void peplist_print_file_my( PEPLIST *peplist, FILENAME filename, char *allele, char *name )

{
        FILE    *fp;
        PEPLIST *pl;
	WORD	shortname;
	int	i,n;

	strncpy( shortname, name, 15 );
        shortname[15] = 0;
        for ( i=0;i<15;i++ )
                if ( shortname[i] == ' ' )
                        shortname[i] = '_';
               
        if ( ( fp = fopen( filename, "w" ) ) == NULL ) {
                printf( "Cannot open file %s\n", filename );
                exit( 1 );
        }      
               
        for ( n=1,pl=peplist; pl; pl=pl->next, n++ ) 
                fprintf( fp, "%i %s %f %f %f %s %s\n", pl->nn, pl->pep, pl->score, pl->nM, pl->rank, shortname, allele );
               
        fclose( fp );
}         

PEPLIST	*nnmultipred( PEPLIST *peplist, FILENAME dir, char *syn )

{
	LINE		cmd;
	FILENAME	filename;
	PEPLIST		*peplist_pred;
	LINE		pfixcmd;

	if ( p_exprefix )
		strcpy( pfixcmd, "" );
	else
		sprintf( pfixcmd, "-px %s", p_rdir );

	sprintf( filename, "%s/0.dat",  dir );

	peplist_print_file( peplist, filename, 0 );

	sprintf( cmd, "cat %s | gawk '{print substr($1,1,9)}' > %s/pep\n", filename, dir );

        doit( cmd );

	sprintf( cmd, "%s/bin/nnlinplayer -x -seqinp %s %s -pt 0 -s %s %s | grep -v \"#\" > %s/pred\n", 
		p_rdir, pfixcmd, ( p_sdev ? "-sdev 1" : ""), syn, filename, dir );
	
	doit( cmd );

	if ( p_sdev )
		sprintf( cmd, "paste %s/pep %s/pred | gawk '{print $1, $2, $4}'",
			dir, dir );
	else
		sprintf( cmd, "paste %s/pep %s/pred | gawk '{print $1, $2}'", 
			dir, dir );

	peplist_pred = peplist_pread( cmd );

	return( peplist_pred );
}

void	peplist_assign_ic50( PEPLIST  *peplist )

{

	PEPLIST *pl;
	float	aff;


	for ( pl = peplist; pl; pl=pl->next ) {

                aff = pl->score;

                pl->nM = exp( (1-aff)*log(50000));
	
	}

}

void	 peplist_assign_rank( PEPLIST *peplist, THRLIST *thrlist )

{

	PEPLIST *pl;
	THRLIST *thrl;

	for ( pl = peplist; pl; pl=pl->next ) {

		if ( thrlist ) {
                        thrl = find_thr( thrlist, pl->score );

                        if ( thrl ) {

                                pl->rank = thrl->thr;
                        }
                        else {

                                pl->rank = 50.0;
                        }

		}
	}	
}

void	print_results( PEPLIST	*peplist, char *name, char *hla_allele, THRLIST *thrlist, int ic50 )

{
	int             n, nSB, nWB;
        WORD            shortname;
	float		aff;
	PEPLIST		*pl;
	float		nM;
	float		sdev;
	int		i;
        WORD            hlaname;
        WORD            hlashort;
	WORD		allele_freq;
	WORD		name1;
	THRLIST		*thrl;
	char		**words;
	int		nw;
	char		loci;

	nSB = 0;
	nWB = 0;

	strncpy( shortname, name, 15 );
	shortname[15] = 0;
	for ( i=0;i<15;i++ )
		if ( shortname[i] == ' ' )
			shortname[i] = '_';

	if ( strncmp( hla_allele, "HLA-", 4 ) == 0 ) {

		words = cmatrix( 0, 2, 0, WORDSIZE );

		strncpy( words[0], hla_allele, 5 );
		words[0][5] = 0;

		strncpy( words[1], hla_allele+5, strlen( hla_allele ) - 5 );
		words[1][strlen( hla_allele ) - 5] = 0;

		sprintf( hlaname, "%s*%s", words[0], words[1] );

		cmatrix_free( words, 0, 2, 0, WORDSIZE );
	}
	else if ( strncmp( hla_allele, "Mamu", 4 ) == 0 ) {

		def_sepchars( ":" );
                words = split( hla_allele, &nw );

		if ( nw != 2 ) {
                        printf( "Error. Inconsistent allele name %s\n", hla_allele );
                        exit( 1 );
                }

		sprintf( hlaname, "%s*%s", words[0], words[1] );

		cmatrix_free( words, 0, nw-1, 0, WORDSIZE );
	}
	else if ( strncmp( hla_allele, "Patr", 4 ) ==  0 || strncmp( hla_allele, "Gogo", 4 ) ==  0 ) {

		words = cmatrix( 0, 1, 0, WORDSIZE );

		strncpy( words[0], hla_allele, 6 );
		words[0][6] = 0;
		strncpy( words[1], hla_allele+6, 4 );
		words[1][4] = 0;

		sprintf( hlaname, "%s*%s", words[0], words[1] );

		cmatrix_free( words, 0, 1, 0, WORDSIZE );
        }
	else if ( strncmp( hla_allele, "SLA-", 4 ) == 0 || strncmp( hla_allele, "BoLA-", 5 ) == 0 ) {

		for ( i=0; i<strlen(hla_allele); i++)
			if ( hla_allele[i] == ':' )
				hlaname[i] = '*';
			else
				hlaname[i] = hla_allele[i];

		hlaname[strlen(hla_allele)] = 0;

        }      
	else if ( strncmp( hla_allele, "H-2", 3 ) == 0 ) {
                sprintf( hlaname, "%s", hla_allele );
        }      
	else if ( strcmp( hla_allele, "USER_DEF" ) == 0 ) {
		sprintf( hlaname, "%s", hla_allele );
	}
	else {
		printf( "Unrecognized HLA allele %s. Exit\n", hla_allele );
		exit( 1 );
	}

	if ( ic50 ) {
		printf( "# Affinity Threshold for Strong binding peptides %7.3f\n", p_hthr );
		printf( "# Affinity Threshold for Weak binding peptides %7.3f\n", p_lthr );
	}

	if ( thrlist ) {
		printf( "# Rank Threshold for Strong binding peptides %7.3f\n", p_hthr_r );
		printf( "# Rank Threshold for Weak binding peptides %7.3f\n", p_lthr_r );
	}

	printf( "-----------------------------------------------------------------------------------\n" );

	printf( "%5s %12s %12s %15s %13s",
                        "pos", "HLA", "peptide", "Identity", "1-log50k(aff)" );

	if ( ic50 )
		printf( " %12s", "Affinity(nM)" );

	if ( thrlist ) 
		printf( " %8s", "%Rank" );

	if ( peplist->aff > -9 ) 
		printf( " %12s", "Exp_Binding" );

	printf( " %10s\n", "BindLevel" );
       
        printf( "-----------------------------------------------------------------------------------\n" );

	for ( n=1,pl = peplist; pl; pl=pl->next,n++ ) {

		aff = pl->score;

		if ( aff < p_thr)
			continue;

		printf( "%5i ", pl->nn );

		printf( "%12s %12s %15s %13.3f", hlaname, pl->pep, shortname, aff );

		if ( ic50 ) 
			printf( " %12.2f", pl->nM );

		if ( p_sdev )
			printf( " %8.3f", sdev );

		if ( thrlist ) 
			printf( " %8.2f", pl->rank );

		if ( pl->aff > -9 )
			printf( " %12.3f", pl->aff );

		if ( ( ic50 && pl->nM <= p_hthr ) || ( thrlist && pl->rank <= p_hthr_r ) ) {

			nSB++;
			printf( " <= SB" );

		}
		else if ( ( ic50 && pl->nM <= p_lthr ) || ( thrlist && pl->rank <= p_lthr_r ) ) {

			nWB++;
			printf( " <= WB" );

		}

		printf( "\n" );

	}

	printf( "-----------------------------------------------------------------------------------\n" );

	printf( "\nProtein %s. Allele %s. Number of high binders %i. Number of weak binders %i. Number of peptides %i\n\n",
                        shortname, hlaname, nSB, nWB, n-1 );

	if ( p_w ) {

		if ( strncmp( hlaname, "HLA-", 4 ) == 0 ) {

			strncpy( allele_freq, hlaname+4, strlen(hlaname) - 4 );

			allele_freq[strlen(hlaname) - 4] = 0;

			printf( "Link to Allele Frequencies in Worldwide Populations <a href=\"%s%s\" target=_blank >HLA-%s</a>\n",
				"http://www.allelefrequencies.net/hla6006a.asp?hla_selection=", allele_freq, allele_freq );
		}

	}

	printf( "-----------------------------------------------------------------------------------\n" );
}

void	peplist_addhlapseudo( PEPLIST *peplist, char *hlapseudo )

{
	PEPLIST *pl;

	for ( pl=peplist; pl; pl=pl->next ) {

		sprintf( pl->pep, "%s%s", pl->pep, hlapseudo );

	}
}

main( int argc, char *argv[] )

{
	int		pid;
	FILENAME	filename, filename2;
	LINE		cmd;
	int		nlines, np;
	FILENAME	tmpdir;
	PEPLIST		*hla_pseudolist, *entry;
	FILENAME	hla_pseudofile;
	FILENAME	hlapseudo;
	FSALIST		*fsalist, *fsa;
	char		**alleles;
	int		nall, na;
	FILENAME	syn;
	PEPLIST		*pl, *pl3, *peplist, *peplist2, *peplist3, *peplist_all;
	float		aff;
	int		n, nn;
	FILENAME        user;
	int		fullMHC;
	LINELIST	*linelist;
	WORD		neighbor;
	float		pcc, dist;
	FILENAME	thrfilename;
	THRLIST		*thrlist;
	NAMELIST	*whitelist;
	int		ic50;
	LINELIST	*versionlist, *vl;
	int		*lvec, nlen, nl, plen, lmin, lmax;

	if ( strcmp( envir("USER", user), "www" ) == 0 ) 
		spparse( &argc, &argv, param, 0, "[fastafile/peptidefile]" );
	else
		pparse( &argc, &argv, param, 0, "[fastafile/peptidefile]" );

	if ( strlen(p_hla) == 0 ) {

		printf( "An MHC allele must be specified using the -a option. Exit\n" );
		exit( 1 );
	}
		
       	if ( strlen(p_f) == 0 && ! p_list && argc < 2 ) {
                printf( "You most specify input file either by using -f or as %s [args] fastafile\n", argv[0] );
                printf( "Usage: netMHCpan [-h] [args] [fastafile]\n" );
                exit( 1 );
        }

	set_pep_verbose( 0 );
	set_list_verbose( 0 );
	set_fsa_verbose( 0 );

	pid = getpid();

        if ( p_w && p_xlsdump )
                sprintf( p_xlsfile, "%s/%i_%s", "/usr/opt/www/pub/CBS/services/NetMHCpan/tmp", pid, "NetMHCpan.xls" );

	versionlist = linelist_read( p_version );

	printf( "\n" );
	for ( vl = versionlist; vl; vl=vl->next )
		printf( "# %s\n", vl->line );
	printf( "\n" );

	linelist_free( versionlist );

	if ( strlen( p_tmpdir ) == 0 ) {
		sprintf( tmpdir, "%i", pid );
	}
	else {
		sprintf( tmpdir, "%s/%i", p_tmpdir, pid );
	}

	sprintf( cmd, "mkdir -p %s\n", tmpdir);

	doit( cmd );

	np = 0;
	nlines = 0;

	if ( strlen( p_hlaseq ) > 0 ) {

		sprintf( cmd, "%s/bin/mhcfsa2psseq %s", p_rdir, p_hlaseq );

		if ( p_verbose ) {
			printf( "# %s\n", cmd );
		}

		linelist = linelist_pread( cmd );

		if ( ! linelist ) {
			printf( "Error. Cannot make pseudo sequence from file %s\n", p_hlaseq );
			exit( 1 );
		}

		sscanf( linelist->line, "%s", hlapseudo );

		if ( strpos( hlapseudo, '-' ) >= 0 ) {
			printf( "Error. MHC input sequence incomplete. Pseudo sequence %s contains gaps\n", 
					hlapseudo );
			exit( 1 );
		}

		fullMHC = 1;
		alleles = split( "USER_DEF", &nall );

		printf( "# Use user MHC pseudo sequence %s %s\n", hlapseudo, alleles[0] );
	}
	else {
		fullMHC = 0;
		alleles = split( p_hla, &nall );

        	if ( p_w && nall > 20 ) {
                	printf( "# Worning. Only 20 alleles allowed per submission. Remaining alleles skipped\n" );
                	nall = 20;
        	}

		sprintf( hla_pseudofile, "%s", p_hlapseudo );
      
		hla_pseudolist = peplist_read( hla_pseudofile );

	}

	if ( p_list ) {

		for ( pl = hla_pseudolist; pl; pl=pl->next ) 
			printf( "%s\n", pl->pep );

		exit( 0 );

	}

	whitelist = namelist_read( p_whitelist );

	if ( ! whitelist ) {
		printf( "Error. Cannot read whitelist from file %s\n", p_whitelist );
		exit( 1 );
	}

	fsalist = NULL;
	peplist = NULL;

	if ( p_inptype == 1 )
		p_pep = 1;

	if ( argc == 2 ) {
		if ( p_pep ) {
			if ( ( peplist = peplist_read_XX( argv[1] ) ) == NULL ) {
				printf( "Error. No Peptide entries  read from %s\n", argv[1] );
				exit( 1 );
			}
	}
		else {
			if ( ( fsalist = fsalist_read( argv[1] ) ) == NULL ) {
				printf( "Error. No FASTA entries read from %s\n", argv[1] );
                       		exit( 1 );
			}
		}
	}
	else {
		if ( p_pep ) {
			if ( ( peplist = peplist_read_XX( p_f ) ) == NULL ) {
				printf( "Error. No Peptide entries  read from %s\n", p_f );
				exit( 1 );
			}
		}
		else {
			if ( ( fsalist = fsalist_read( p_f ) ) == NULL ) {
				printf( "Error. No FASTA entries read from %s\n", p_f );
                       		exit( 1 );
			}
		}
	}

	if ( peplist ) {
        	printf( "# Input is in PEPTIDE format\n" );
		peplist_check_replace( peplist, "ARNDCQEGHILKMFPSTWYVX", 'X' );

		lvec = ivector( 0, 0 );
		nlen = 1;

		lvec[0] = peplist->len;
	}
	else {
		printf( "# Input is in FSA format\n" );
		fsalist = fsalist_check_names( fsalist );
		fsalist_checkalphabet_replace( fsalist, "ARNDCQEGHILKMFPSTWYVX", 'X' );

		lvec = ivector_split( p_len, &nlen ); 

		lmin = ivector_min( lvec, nlen );
		lmax = ivector_max( lvec, nlen );

#ifdef NOTNOW
		if ( lmin < 8 || lmax > 11 ) {
			printf( "Error peptide length %s -l option muste be between 8-11\n", p_len );
			exit( 1 );
		}
#endif
        
		printf( "\n# Peptide length %s\n", p_len );
	}

	for ( na=0; na<nall; na++ ) {

		if ( ! fullMHC ) {

			if ( (entry = peplist_find( hla_pseudolist, alleles[na] ) ) == NULL ) {
				printf( "Error. %s cannot be found in hla_pseudo list %s\n", 
					alleles[na], hla_pseudofile );
				exit( 1 );
			}
       
			sscanf( entry->line, "%*s %s", hlapseudo );
		}
		
		if ( ! fullMHC ) {

			sprintf( thrfilename, p_thrfmt, alleles[na] );

			if ( filereadable( thrfilename ) ) {

				thrlist = thrlist_read( thrfilename );

				if ( ! thrlist ) {
					printf( "Error. Cannot read Threshold file %s. Exit\n", thrfilename );
					exit( 1 );
				}
			}
			else {
				printf( "Error. Threshold file %s does not exist. Exit\n", thrfilename );
				exit( 1 );
			}
		}
		else
			thrlist = NULL;

		sprintf( cmd, "%s/bin/estimate_PCC %s", p_rdir, hlapseudo );

		linelist = linelist_pread( cmd );

		if ( ! linelist ) {
			printf( "Error. Cannot read line from %s. Exit\n", cmd );
			exit( 1 );
		}

		sscanf( linelist->line, "%f %s %f", &dist, neighbor, &pcc );	

		linelist_free( linelist );

		printf( "\n%s : Estimated prediction accuracy %6.3f (using nearest neighbor %s)\n\n", 
			alleles[na], pcc, neighbor );

		ic50 = ( p_ic50 || namelist_find( whitelist, alleles[na] ));

		sprintf( syn, p_syn );

		if ( ! peplist ) {

			if ( p_xlsdump ) {
                                sprintf( cmd, "touch %s/%s.all; rm -f %s/%s.all", tmpdir, alleles[na], tmpdir, alleles[na] );

                                doit( cmd );
                        }
	
        		for ( np=0,fsa=fsalist; fsa; fsa=fsa->next, np++ ) {

				peplist_all = NULL;

				for ( nl = 0; nl<nlen; nl++ ) {

					plen = lvec[nl];

					peplist = fsalist2pep_single( fsa, plen, 1 );

					if ( ! peplist ) 
						continue;
				
					sprintf( filename, "%s/%i.pep", tmpdir, np );

					peplist_print_file( peplist, filename, 0 );
			
					nn = 0;

					if ( plen != 9 ) {

						if ( plen == 8 )
							sprintf( cmd, "cat %s | %s/bin/pep28mer", filename, p_rdir );
						else
							sprintf( cmd, "cat %s | %s/bin/pep2lmer", filename, p_rdir );

						peplist2 = peplist_pread( cmd );

						peplist_addhlapseudo( peplist2, hlapseudo );

						peplist3 = nnmultipred( peplist2, tmpdir, syn );

						pl = peplist;
						aff = 0;
						n=0;

						if ( plen > 9 ) {

							for ( pl3 = peplist3; pl3 && pl; pl3=pl3->next ) {

								aff += pl3->score;
								n++;

								if ( n == 6 ) {
									pl->score = aff/6.0;
									pl->nn = nn++;
									pl = pl->next;
									n = 0;
									aff = 0.0;
								}
							}
						}
						else {

							for ( pl3 = peplist3; pl3 && pl; pl3=pl3->next ) {

								aff += pl3->score;
								n++;

								if ( n == 5 ) {
									pl->score = aff/5.0;
									pl->nn = nn++;
									pl = pl->next;
									n = 0;
									aff = 0.0;
								}
							}
						}

						peplist_free( peplist3 );
						peplist_free( peplist2 );
					}
					else {

						peplist2 = peplist;

						peplist_addhlapseudo( peplist2, hlapseudo );

						peplist3 = nnmultipred( peplist2, tmpdir, syn );

						peplist_free( peplist );

						peplist = peplist3;

						for ( pl = peplist; pl; pl=pl->next ) 
							pl->nn = nn++;

					}

					if ( ! peplist_all )
						peplist_all = peplist;
					else {
						for ( pl = peplist_all; pl->next; pl=pl->next );

						pl->next = peplist;
					}

				}

				peplist = peplist_all;

#ifdef NOTNOW
				for ( pl = peplist_all; pl->next; pl=pl->next )
					printf( "TEST %s %f %f\n", pl->pep, pl->score, pl->aff );
#endif

				if ( ! peplist )
					continue;

				if ( ic50 )
					peplist_assign_ic50( peplist );

				peplist_assign_rank( peplist, thrlist );

				if ( p_xlsdump ) {
					sprintf( filename2, "%s/%s_%i.pred", tmpdir, alleles[na], np );
					peplist_print_file_my( peplist, filename2, alleles[na], fsa->name );

					sprintf( cmd, "cat %s/%s_%i.pred >> %s/%s.all", tmpdir, alleles[na], np, tmpdir, alleles[na] );

					doit( cmd );
				}

				if ( p_sort )
					peplist = peplist_sort( peplist, 1 );

				print_results( peplist, fsa->name, alleles[na], thrlist, ic50 );

				nlines++;

				peplist_free( peplist );
				peplist = NULL;
			}

#ifdef NOTNOW
			if ( p_xlsdump ) {
				sprintf( cmd, "cat %s/%s_*.pred > %s/%s.all", tmpdir, alleles[na], tmpdir, alleles[na] );

				doit( cmd );
			}
#endif
		}
		else {

			plen = lvec[0];

			for ( pl=peplist; pl; pl=pl->next ) {

				if ( pl->len != plen ) {
					printf( "Error. Peptide lenght must be equal for all peptides %i %s. Exit\n", pl->len, p_len );
					exit( 1 );
				}

				if ( na == 0 )
					pl->aff = pl->score;

			}

			sprintf( filename, "%s/%i.pep", tmpdir, np );

			peplist_print_file( peplist, filename, 0 );

			if ( plen != 9 ) {

				if ( plen > 9 )
					sprintf( cmd, "cat %s | %s/bin/pep2lmer", filename, p_rdir );
				else
					sprintf( cmd, "cat %s | %s/bin/pep28mer", filename, p_rdir );

				peplist2 = peplist_pread( cmd );

				peplist_addhlapseudo( peplist2, hlapseudo );

				peplist3 = nnmultipred( peplist2, tmpdir, syn );

				pl = peplist;
				aff = 0;
				n=0;

				if ( plen > 9 ) {

					for ( pl3 = peplist3; pl3 && pl; pl3=pl3->next ) {

						aff += pl3->score;
						n++;

						if ( n == 6 ) {
							pl->score = aff/6.0;
							pl = pl->next;
							n = 0;
							aff = 0.0;
						}
					}
				}
				else {

					for ( pl3 = peplist3; pl3 && pl; pl3=pl3->next ) {

						aff += pl3->score;
						n++;

						if ( n == 5 ) {
							pl->score = aff/5.0;
							pl = pl->next;
							n = 0;
							aff = 0.0;
						}
					}
				}

				peplist_free( peplist3 );
				peplist_free( peplist2 );
			}
			else {

				peplist2 = peplist;

				peplist_addhlapseudo( peplist2, hlapseudo );

				peplist3 = nnmultipred( peplist2, tmpdir, syn );

				if ( peplist->aff > -9 ) {
					for ( pl=peplist, pl3=peplist3; pl && pl3; pl=pl->next, pl3=pl3->next )
						pl3->aff = pl->aff; 
				}

				peplist_free( peplist );

				peplist = peplist3;

			}


			if ( ic50 )
				peplist_assign_ic50( peplist );

			peplist_assign_rank( peplist, thrlist );

			if ( p_xlsdump ) {
				sprintf( filename2, "%s/%s.all", tmpdir, alleles[na] );
				peplist_print_file_my( peplist, filename2, alleles[na], "PEPLIST" );
			}

			if ( p_sort )
				peplist = peplist_sort( peplist, 1 );

			print_results( peplist, "PEPLIST", alleles[na], thrlist, ic50 );

		}

		if ( thrlist )
			thrlist_free( thrlist );

	}

#define GAWKCMD "gawk -v thr2=2.0 -v thr1=500 '{ printf( \"%i\\t%s\\t%s\", $1, $2, $6); ave=0.0; nb=0;n=0; for ( i=0; i<NF; i+=7 ) {printf( \"\\t%6.4f\t%6.4f\t%6.4f\", $(i+3), $(i+4), $(i+5)); ave +=$(i+3); nb += ( ( $(i+4) > 0 && $(i+4) < thr1 ) || ( $(i+5) > 0 && $(i+5) <= thr2)); n++} if ( n==0 ) { n=1} printf( \"\\t%6.4f\\t%6i\\n\", ave/n, nb) }'"

#define GAWKCMD_HEADER "gawk '{ printf( \"\\t\\t\" ); for ( i=0; i<NF; i+=7 ) {printf( \"\\t\t%s\t\", $(i+7) )} printf( \"\\t\\n\") }'"

#define GAWKCMD_HEADER2 "gawk '{ printf( \"%s\\t%s\\t%s\", \"Pos\", \"Peptide\", \"ID\" ); for ( i=0; i<NF; i+=7 ) {printf( \"\\t1-log50k\tnM\tRank\")} printf( \"\\tAve\\tNB\\n\") }'"

	if ( p_xlsdump ) {

		sprintf( cmd, "paste %s/*.all | head -1 |%s > %s", tmpdir, GAWKCMD_HEADER, p_xlsfile );

		doit( cmd );

		sprintf( cmd, "paste %s/*.all | head -1 |%s >> %s", tmpdir, GAWKCMD_HEADER2, p_xlsfile );

		doit( cmd );

		sprintf( cmd, "paste %s/*.all | %s >> %s", tmpdir, GAWKCMD, p_xlsfile );

		doit( cmd );

		if ( p_w )
			printf( "Link to output xls file <a href=\"%s/%i_%s\">NetMHCpan_out.xls</a>\n",
				"http://www.cbs.dtu.dk/services/NetMHCpan/tmp/", pid, "NetMHCpan.xls" );
        }

	if ( ! p_dirty ) {
		sprintf( cmd, "rm -rf %s\n", tmpdir );
		doit( cmd );
	}

	exit( 0 );
}
