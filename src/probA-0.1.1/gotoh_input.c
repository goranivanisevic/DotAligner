#include "pfgoto.h"
#include "config.h"


char prog_name[10] = "probA";



float BETA=1;            /* default value of BETA */
float ENDGAP=0;          /* default value of endgaps */
char MAT_SER[20]="gon";  /* variable stores the selected matrix series */
			 /* default: gonnet_series */
float DISTANCE=-1;       /* user can specify a pam distance (gon, pam) or an
			    observed identity (blo) to determine which matrix
			    out of a scoring matrix series to use:
			    default is -1 => probA specifies the score matrix
			    if DISTANCE is positive the scoring matrix
			    closed to the value of DISTANCE is selected */ 
long Nr=1;                /* number of alignments generated by stochastic
			    backtracking */
/* char TRACK[30];           stores the filename containing the alignment 
			     encoded as one string of symbols;
			     option`-calc_s' */
/* char *track;              stores the alignment encoded as one string of
			     symbols; option`-calc_s' */



sequ *read_seq(FILE *ptr);
/*  char *get_line(FILE *fp); now in pfgoto.h*/ 
sequ *lies_seq(FILE *spt);


/*--------------------------------------------------------------------------*/
void start(int argc, char **argv)
{
  int c;
  float Temp;
  float end;
  float di;
  char buf[20];
  /*  FILE *Track; */

  
  while ((c = getopt_long(argc, argv, shortopts, long_options, &option_index)) != -1)

    switch (c)
      {
      case 0:
	/* If this option set a flag and accepts no argument,
	   do nothing else now. */
	if ((long_options[option_index].flag != 0) &&
	    (long_options[option_index].has_arg == 0))
	  break;
	
	if(strcmp(long_options[option_index].name,"score_matrix")==0)
	  {
	    if(sscanf(optarg, "%s", buf) != EOF){
	      strcpy(MAT_SER,buf);
	      /*  printf("start: matrix series = %s\n",MAT_SER); */ 
	    }
	    /*  printf ("option %s", long_options[option_index].name); */
/*  	    if (optarg) */
/*  	      printf (" with arg %s", optarg); */
/*  	    printf ("\n"); */
/*  	    break;  */ 
	  } 
	
	
	if(strcmp(long_options[option_index].name,"pam")==0)
	  {
	    if(sscanf(optarg, "%f", &di) != EOF)
	      {
		DISTANCE = di;
		/*  printf("start: pam or observed id = %e\n",DISTANCE); */
	      }
	    /*  printf ("option %s", long_options[option_index].name); */
/*  	    if (optarg) */
/*  	      printf (" with arg %s", optarg); */
/*  	    printf ("\n"); */
/*  	    break; */ 
	  }
       

	if(strcmp(long_options[option_index].name,"endgaps")==0)
	  {
	    if(sscanf(optarg, "%f", &end) != EOF)
	      {
		ENDGAP = end;
		/*  printf("start: endgap penalty = %e\n",ENDGAP); */ 
	      }
	    /*  printf ("option %s", long_options[option_index].name); */
/*  	    if (optarg) */
/*  	      printf (" with arg %s", optarg); */
/*  	    printf ("\n"); */
/*  	    break; */ 
	  }
	
	/*  if(strcmp(long_options[option_index].name,"calc_s")==0) */ 
/*  	  { */ 
/*  	    if(sscanf(optarg, "%s", buf) != EOF) */ 
/*  	      { */ 
/*  		strcpy(TRACK,buf); */ 
/*  		printf("start: TRACK = %s\n",TRACK);  */ 
/*  		Track = fopen(TRACK,"r"); */ 
/*  		track = get_line(Track); */ 
/*  		fclose(Track); */ 
		
/*  	      } */ 
/*  	  }  */ 
	
	break;

      case 'T':
	if(sscanf(optarg, "%f", &Temp) != EOF) {
	  BETA = 1.0/(double) Temp;
	  
	  /*  printf("start: T= %e\n",Temp); */
	  Temp=0;
	}
	break;

      case 'N':
	if(sscanf(optarg, "%ld", &Nr) != EOF) {
	  /*  printf("start: N = %ld\n",Nr); */
	}
	break;
	/*new*/
      case 'h':
          usage(0);
	  break;

      case 'V':
          printf ("%s %s\n", PACKAGE, VERSION);
          exit(0);

	/*new*/
	  
      case '?':
	/* getopt_long already printed an error message. */
	break;
	
      default:
	printf("default");
	usage(0);
      }


  /* Instead of reporting `--verbose'
     and `--brief' as they are encountered,
     we report the final status resulting from them. */
  if (verbose_flag)
    puts ("verbose flag is set");
  
  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }
  
 
}
/*--------------------------------------------------------------------------*/
sequ *input(int argc, char **argv)
{
  sequ *sq;  
  int istty;

  istty=isatty(fileno(stdout))&&isatty(fileno(stdin));
  if(istty)
    {
      sq=lies_seq(stdin);
    }
  else
    { 
      sq=read_seq(stdin);
    }
  return(sq);  
}

/*--------------------------------------------------------------------------*/
void *myspace(unsigned size)
{/*
  void *pointer;

   if ( (pointer = (void*) calloc(1,size)) == NULL)
     {
       fprintf(stderr,"SPACE: requested size: %d\n",size);
     }
     return pointer;*/
  void *pointer;
    
  if ( (pointer = (void *) calloc(1, size) ) == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"SPACE: requested size: %d\n", size);
      nrerror("SPACE allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      nrerror("SPACE allocation failure -> no memory");
  }
  return  pointer;
  
}
/*--------------------------------------------------------------------------*/

char *get_line(FILE *fp)
{
  char s[512], *line, *cp;

  line = NULL;
  do {
    if (fgets(s, 512, fp)==NULL) break;
    cp = strchr(s, '\n');
    if (cp != NULL) *cp = '\0';
    if (line==NULL)
      line = (char*) calloc((strlen(s)+1), sizeof(char));
    else
      line = (char*) realloc(line, (strlen(s)+strlen(line)+1));
    strcat(line,s);
  } while(cp==NULL);

  return(line);
}

/*--------------------------------------------------------------------------*/
sequ *lies_seq(FILE *spt)
{
  sequ *twos;
  char *sname;
  char *seq;


  twos = (sequ*) calloc(2, sizeof(sequ));

    
  printf("name of the first sequence:");
  sname=get_line(stdin);
  twos[0].name = strdup(sname);
  free(sname);

  printf("enter the first sequence:\n");
  seq=get_line(stdin);
  twos[0].seq = strdup(seq);
  free(seq);
  printf("\n");
   
  printf("name of the secound sequence:");
  sname=get_line(stdin);
  twos[1].name = strdup(sname);
  free(sname);
  
  printf("enter the secound sequence:\n");
  seq=get_line(stdin);
  twos[1].seq = strdup(seq);
  free(seq);
  printf("\n");
  return(twos);
  
}

/*--------------------------------------------------------------------------*/
sequ *read_seq(FILE *ptr)
{
  int i;
  char *line, workstr[100], *strg;
  sequ *seqs;

  seqs=(sequ*) calloc( 2, sizeof(sequ));
  
  
  i=0;
 
  line=get_line(ptr);
  while (i<2) {  /* until we have 2 sequences */

    while(line && ((line[0]=='\0') || isspace(line[0]))){
      free(line);
      line=get_line(ptr);
    } 

    if (!line) break;
    
    if (line[0] == '>'){
      sscanf(line,">%99s", workstr);
      seqs[i].name = strdup(workstr);
      free(line);
      line=get_line(ptr);
    }

    if (line && isalpha(line[0])) {
      strg=(char*) calloc((strlen(line)+1), sizeof(char));
      sscanf(line,"%s", strg);
      seqs[i].seq = strdup(strg);
      free(line);
      free(strg);
      line=get_line(ptr);
      ++i; 
    }
  }
  if(line)
  free(line);
   /*
  printf("read_seq:%s\n%s\n\n%s\n%s\n",seqs[0].name,seqs[0].seq,seqs[1].name,seqs[1].seq);
  */ 
  if (i<2) nrerror("too few sequences in input\n");
  
  return(seqs);
}

/*--------------------------------------------------------------------------*/

void usage(int status)
{
   /*  nrerror("usage: Gotoh\n\t\t[-T input temperature for the alignement\n\n\t\t[-score_matrix select a scoring matrix series:\n\t\tblo blosum series;\n\t\tblo_p blosum series all positv matrices;\n\t\tpam pam series;\n\t\tpam_p pam series all positv matrices;\n\n\t\t-endgaps treat endgap differenly, the value for the endgap penalty can be modified(optional);");  */

printf("%s - Compute the partition function over all possible alignments between the two input sequences, generate optimal and suboptimal alignments by stochastic backtracking\n", prog_name);

  printf("Usage: %s [OPTION]... [FILE]\n", prog_name);
  printf(
         "Options:\n"
         /*  "-q, --quiet, --silent      be quiet, inhibit PS output\n" */
/*           "--verbose                  print more information\n" */
         "-h, --help       display this help and exit\n"
         "-V, --version    output version information and exit\n"
         "-T               governs relative weight of alignment paths with different scores\n"
	 "-N               number of alignments generated by stochastic backtracking\n"
	 "--endgaps        set the score for terminal gaps, default is 0\n"
	 "--noEg           terminal gaps are treated as gaps inside the alignment\n"
	 "--noPS           suppress the generation of a dot plot\n"
	 "--DNA            input polymers will be treated as nucleic acids\n"
	 "--prot           input polymers will be treated as amino acids\n"
	 "--score_matrix   determines the substitution matrix used for protein alignments\n"
	 "--pam            specify one specific matrix out of a matrix series: set --score_matrix to the matrix of choice and --pam to the pam distance (or identity value) of the respective matrix\n"
         );
  printf("\nFormat of the input file:\n"); 
  printf("> name_of_sequence_1\n");
  printf("sequence_1\n\n");
  printf("> name_of_sequence_2\n");
  printf("sequence_2\n\n");
  printf("matrix series:\n");
  printf("gonnet-series:gonnet_40, gonnet_80, gonnet_120, gonnet_160, gonnet_250, gonnet_300, gonnet_350;\n");
  printf("blosum-series:blosum_30, blosum_62, blosum_50, blosum_80;\n");
  printf("pam-series:pam_20, pam_60, pam_120, pam_350;\n");
  
  exit (status);


}

/*--------------------------------------------------------------------------*/











