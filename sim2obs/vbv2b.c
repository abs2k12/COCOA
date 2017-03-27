/*========================================================*/
/*                                                        */
/*  vbv2b.c          2014.03.06      version 0.1.1        */
/*                                                        */
/*  Copyright (C) 2014 Wojtek Pych, CAMK PAN              */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*  Convert B-V into B in the cluster model database.     */
/*                                                        */
/*========================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LLEN 1024

const int nv=24;	/* number of V column		*/
const int bpos=319;	/* position of B-V in the lines	*/

int main(int argc, char *argv[])
{
	char  *iname,
	      *p,
	      line[LLEN],
	      nline[LLEN],
	      b[LLEN];
	int   i, j;
	float v, bv;
	FILE  *inf;
	
  if (argc != 2) 
  { printf("\n\n\tUsage: vbv2b  catalog\n\n"); return(EXIT_FAILURE); }
  
  iname=argv[1];
  
  if ((inf=fopen(iname, "r")) == NULL)
  { perror(iname); return(EXIT_FAILURE); }
  
  for (i=0; feof(inf)==0; i++)
  {
    if (fgets(line, LLEN, inf) == NULL) break;
    if (feof(inf) != 0) break;
    
    strcpy(nline, line);
    
    p=strtok(line, " ");
    for (j=0; j<nv-1; j++) p=strtok(NULL," ");
    sscanf(p, "%f", &v);
    p=strtok(NULL," ");
    sscanf(p, "%f", &bv);
/*    printf("%10.8E %10.8E %10.8E\n", v, bv, v+bv); */
    sprintf(b, "%+10.8E", v+bv);
    memcpy(nline+bpos, b, strlen(b));
    printf("%s", nline);
  }
  
  fclose(inf);
  
  return(EXIT_SUCCESS);
}