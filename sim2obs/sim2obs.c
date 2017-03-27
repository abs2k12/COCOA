/*========================================================*/
/*                                                        */
/*  sim2obs.c           version 0.5.1   2016.01.19        */
/*                                                        */
/*  Copyright (C) 2014-2016 by Wojtek Pych, CAMK PAN      */
/*  pych@camk.edu.pl                                      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*                                                        */
/*========================================================*/

/***************************************************************************/
/*                                                                         */
/*   This program is free software; you can redistribute it and/or modify  */
/*   it under the terms of the GNU General Public License as published by  */
/*   the Free Software Foundation; either version 2 of the License, or     */
/*   (at your option) any later version.                                   */
/*                                                                         */
/*   This program is distributed in the hope that it will be useful,       */
/*   but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*   GNU General Public License for more details.                          */
/*                                                                         */
/*   You should have received a copy of the GNU General Public License     */
/*   along with this program; if not, write to the                         */
/*   Free Software Foundation, Inc.,                                       */
/*   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             */
/*                                                                         */
/***************************************************************************/

/* It is assumed that the X Y MAG are written in this order	*/
/* nx < ny < nm							*/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "pfitsio.h"

#define LLEN 1024
#define DG 1.0e-6	/* counts on the border for Gaussian	*/
#define DM 1.0e-6	/* counts on the border for Moffat	*/

typedef struct
{
  char	filter,
	psf,		/* PSF model				*/
	object[VALUE_SIZE]; /* object name			*/
  short verbose,	/* verbose level 0,1,..			*/
	naxis1,		/* image size [pixels]			*/
	naxis2,		/* image size [pixels]			*/
	nx,		/* X-coordinate column			*/
	ny,		/* y-coordinate column			*/
	nm,		/* magnitude column			*/
	noise;		/* type of noise			*/
  float	dist,		/* distance				*/
	gain,		/* photons/ADU				*/
	sat,		/* detector saturation level		*/
	pixscale,	/* pixel scale of the camera [arcsec]	*/
	seeing,		/* [arcsec]				*/
	beta,		/* beta parameter for Moffat function	*/
	exp,		/* exposure time			*/
	bkg,		/* background level			*/
	raoff,		/* offset in R.A. [arcseconds]		*/
	decoff;		/* offset in Dec. [arcseconds]		*/  
}
PARAMS;

/*--------------------------------------------------------*/
int read_params(char *pname, PARAMS *par)
{
        char    line[LLEN],
                key[LLEN],
                val[LLEN];
        FILE    *inpf;
        PARAMS  lpar;

/* set default parameter values */
  lpar.filter=' ';
  lpar.naxis1=2048;
  lpar.naxis2=2048;
  lpar.nx=2;
  lpar.ny=3;
  lpar.nm=24;
  lpar.dist=3.3e3;
  lpar.pixscale=0.4;
  lpar.gain=1.0;
  lpar.sat=32000.0;
  lpar.exp=1.0;
  lpar.seeing=1.0;
  lpar.bkg=10.0;
  lpar.psf='M';
  lpar.beta=1.8;
  lpar.noise=1;
  lpar.raoff=0.0;
  lpar.decoff=0.0;
  lpar.verbose=1;

/* read the parameters from input file */
  if ((inpf=fopen(pname, "r")) == NULL)
  { perror(pname); return(EXIT_FAILURE); }

  while (feof(inpf) == 0)
  {
    if (fgets(line, LLEN, inpf) == NULL) break;
    if (feof(inpf) != 0) break;

    sscanf(line, "%s = %s", key, val);

    if (key[0] == '#') continue;
    if (strcasecmp(key, "END") == 0) break;
    else if (strcasecmp(key, "NAXIS1") == 0) sscanf(val, "%hd", &lpar.naxis1);
    else if (strcasecmp(key, "NAXIS2") == 0) sscanf(val, "%hd", &lpar.naxis2);
    else if (strcasecmp(key, "DB-NX") == 0) sscanf(val, "%hd", &lpar.nx);
    else if (strcasecmp(key, "DB-NY") == 0) sscanf(val, "%hd", &lpar.ny);
    else if (strcasecmp(key, "DB-NMAG") == 0) sscanf(val, "%hd", &lpar.nm);
    else if (strcasecmp(key, "OBJECT") == 0) sscanf(val, "%s", lpar.object);
    else if (strcasecmp(key, "DISTANCE") == 0) sscanf(val, "%f", &lpar.dist);
    else if (strcasecmp(key, "FILTER") == 0) sscanf(val, "%c", &lpar.filter);
    else if (strcasecmp(key, "PIXSCALE") == 0) sscanf(val, "%f", &lpar.pixscale);
    else if (strcasecmp(key, "GAIN") == 0) sscanf(val, "%f", &lpar.gain);
    else if (strcasecmp(key, "SATLEVEL") == 0) sscanf(val, "%f", &lpar.sat);
    else if (strcasecmp(key, "EXPOSURE") == 0) sscanf(val, "%f", &lpar.exp);
    else if (strcasecmp(key, "SEEING") == 0) sscanf(val, "%f", &lpar.seeing);
    else if (strcasecmp(key, "BACKGROUND") == 0) sscanf(val, "%f", &lpar.bkg);
    else if (strcasecmp(key, "PSF") == 0) sscanf(val, "%c", &lpar.psf);
    else if (strcasecmp(key, "M_BETA") == 0) sscanf(val, "%f", &lpar.beta);
    else if (strcasecmp(key, "NOISE") == 0) sscanf(val, "%hd", &lpar.noise);
    else if (strcasecmp(key, "RA_OFFSET") == 0) sscanf(val, "%f", &lpar.raoff);
    else if (strcasecmp(key, "DEC_OFFSET") == 0) sscanf(val, "%f", &lpar.decoff);
    else if (strcasecmp(key, "VERBOSE") == 0) sscanf(val, "%hd", &lpar.verbose);
    else printf("\n\tWARNING! read_params(): unknown parameter %s\n", key);
  }

  fclose(inpf);

  if (lpar.verbose > 1)
  {
    printf("Reading %s\n", pname);
    printf("-----------\n");
    printf("NAXIS1     = %d\n", lpar.naxis1);
    printf("NAXIS2     = %d\n", lpar.naxis2);
    printf("DB-NX      = %d\n", lpar.nx);
    printf("DB-NY      = %d\n", lpar.ny);
    printf("DB-NMAG    = %d\n", lpar.nm);
    printf("OBJECT     = %s\n", lpar.object);
    printf("DISTANCE   = %g\n", lpar.dist);
    printf("FILTER     = %c\n", lpar.filter);
    printf("PIXSCALE   = %f\n", lpar.pixscale);
    printf("GAIN       = %f\n", lpar.gain);
    printf("SATLEVEL   = %f\n", lpar.sat);
    printf("EXPOSURE   = %f\n", lpar.exp);
    printf("SEEING     = %f\n", lpar.seeing);
    printf("BACKGROUND = %g\n", lpar.bkg);
    printf("PSF        = %c\n", lpar.psf);
    printf("M_BETA     = %g\n", lpar.beta);
    printf("NOISE      = %d\n", lpar.noise);
    printf("RA_OFFSET  = %g\n", lpar.raoff);
    printf("DEC_OFFSET = %g\n", lpar.raoff);
    printf("VERBOSE    = %d\n", lpar.verbose);
    printf("---------\n");
  }

  *par=lpar;

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int read_db(PARAMS par, char *dbname, float **x, float **y, float **m)
{
      char  line[LLEN],		/* buffor for database read		*/
            *p;			/* pointer to database line element	*/
      int   i, j;		/* loop numerator			*/
      float *lx, *ly, *lm;
      FILE  *inf;

  if (par.verbose > 0) printf("\n  %s\n", dbname);

  if ((inf=fopen(dbname, "r")) == NULL)
  { perror(dbname); return(0); }

  if ((lx=(float *)calloc(1, sizeof(float))) == NULL)
  { perror("\n\tERROR! read_db(): calloc(lx)"); return(0); }

  if ((ly=(float *)calloc(1, sizeof(float))) == NULL)
  { perror("\n\tERROR! read_db(): calloc(ly)"); return(0); }

  if ((lm=(float *)calloc(1, sizeof(float))) == NULL)
  { perror("\n\tERROR! read_db(): calloc(lm)"); return(0); }

  for (i=0; feof(inf) == 0; i++)
  {
    if (fgets(line, LLEN, inf) == NULL) break;
    if (feof(inf) != 0) break;

    if ((lx=(float *)realloc(lx, (i+1)*sizeof(float))) == NULL)
    { perror("\n\tERROR! read_db(): realloc(lx)"); return(0); }

    if ((ly=(float *)realloc(ly, (i+1)*sizeof(float))) == NULL)
    { perror("\n\tERROR! read_db(): reallocly)"); return(0); }

    if ((lm=(float *)realloc(lm, (i+1)*sizeof(float))) == NULL)
    { perror("\n\tERROR! read_db(): realloc(lm)"); return(0); }

    p=strtok(line, " ");
    for (j=0; j<par.nx-1; j++) p=strtok(NULL," ");
    sscanf(p, "%f", &lx[i]);
    for (; j<par.ny-1; j++) p=strtok(NULL," ");
    sscanf(p, "%f", &ly[i]);
    for (; j<par.nm-1; j++) p=strtok(NULL," ");
    sscanf(p, "%f", &lm[i]);
    if (par.verbose > 2) printf("%g %g %g\n", lx[i], ly[i], lm[i]);
  }

  fclose(inf);

  *x=lx;
  *y=ly;
  *m=lm;

  if (par.verbose > 0) printf("  %d database records\n", i);

  return(i);
}
/*--------------------------------------------------------*/
int proc_coomag(PARAMS par, int ndb, float *x, float *y, float *m)
{
        int   i;
	float dm,	/* distance modulus	*/
	      scale;	/* coordinates scale	*/

const float pcau=206264.81;	/* astronomical units in 1 parsec */

  dm=5.0*log10(par.dist)-5.0;
  scale=1.0/par.pixscale/par.dist*pcau;
  if (par.verbose > 0)
    printf("  dist.modulus= %g  coord.scale= %g [pixels/pc]\n", dm, scale);

  for (i=0; i<ndb; i++)
  {
    x[i]*=scale;
    x[i]+=par.raoff/par.pixscale+(float)par.naxis1/2.0;
    y[i]*=scale;
    y[i]+=par.decoff/par.pixscale+(float)par.naxis2/2.0;
    m[i]+=dm;
    m[i]=pow(10.0, -0.4*m[i]+10.0);
  }

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int pixels(PARAMS par, int ndb, float *x, float *y, float *flux, float ***im)
{
	int	i, j, m, n,	/* loop numerator			*/
		dmax,	/* maksimum distance for gaussian		*/
		nout,	/* number of star outside field of view	*/
		tx, ty,
		ix, iy;	/* integer coordinate of the star center	*/
	float	a,	/* peak counts					*/
		fwhm,	/* FWHM						*/
		sigma, sigma2,
		fr,	/* A to flux ratio				*/
		count,
		cx, cy,
		**lim;	/* local image pointer				*/

  if ((lim=(float **)calloc(par.naxis2, sizeof(float *))) == NULL)
  { printf("\n\tERROR! pixels(): calloc(lim)\n"); return(EXIT_FAILURE); }

  for (i=0; i<par.naxis2; i++)
  {
    if ((lim[i]=(float *)calloc(par.naxis1, sizeof(float))) == NULL)
    { printf("\n\tERROR! pixels(): calloc(lim[i[)\n"); return(EXIT_FAILURE); }

    for (j=0; j<par.naxis1; j++) lim[i][j]=0.0;
  }

  fwhm=par.seeing/par.pixscale;
  sigma2=fwhm*fwhm/8.0/log(2.0);
  sigma=sqrt(sigma2);
  fr=2.0*M_PI*sigma2;
  if (par.verbose > 0) printf("  FWHM= %g  sigma= %g\n", fwhm, sigma);

  if ((par.psf != 'G') && (par.psf != 'M'))
  {
    printf("\n\tERROR! pixels(): unsupported PSF type: %c\n", par.psf);
    return(EXIT_FAILURE);
  }

  nout=0;
  for (i=0; i<ndb; i++)
  {
    ix=(int)floor(x[i]+0.5);
    if ((ix < 0) || (ix >=par.naxis1)) { nout++; continue; }
    iy=(int)floor(y[i]+0.5);
    if ((iy < 0) || (iy >=par.naxis2)) { nout++; continue; }

    a=flux[i]/fr;

    if (par.psf == 'G')
    {
      if (a < DG) continue;
      dmax=(int)floor(sigma*sqrt(2.0*log(a/DG))+0.5);
    }
    else if (par.psf == 'M')
    {
      if (a < DM) continue;
      dmax=(int)floor(fwhm/2.0*sqrt(pow(a/DM, 1/par.beta)-1.0));
    }
    else return(EXIT_FAILURE);

    for (m=-dmax; m<=dmax; m++)
    {
      ty=iy+m;
      if ((ty < 0) || (ty >=par.naxis2)) continue;
      for (n=-dmax; n<=dmax; n++)
      {
	tx=ix+n;
	if ((tx < 0) || (tx >=par.naxis1)) continue;

	cx=x[i]-tx;
	cy=y[i]-ty;
	if (par.psf == 'G')
	  count=a*exp(-(cx*cx+cy*cy)/2.0/sigma2);
	else if (par.psf == 'M')
	{
	  count=a/pow(1.0+4.0*(cx*cx+cy*cy)/fwhm/fwhm, par.beta);
	}
	else return(EXIT_FAILURE);

	lim[ty][tx]+=count;
      }
    }
  }

  if (par.verbose > 0) printf("  %d stars outside field of view\n", nout);

/* scale the image according to EXPOSURE */
  for (i=0; i<par.naxis2; i++)
  {
    for (j=0; j<par.naxis1; j++)
    {
      lim[i][j]*=par.exp;
      lim[i][j]+=par.bkg*par.exp;
    }
  }

  *im=lim;

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int add_noise(PARAMS par, time_t rawtime1, float **im)
{
	int i, j;

	const gsl_rng_type	* T;
	gsl_rng			* r;

  gsl_rng_env_setup();

  T=gsl_rng_default;
  r=gsl_rng_alloc(T);
  
  gsl_rng_set(r, rawtime1);

  for (i=0; i<par.naxis2; i++)
  {
    for (j=0; j<par.naxis1; j++)
      im[i][j]=gsl_ran_poisson(r, (int)im[i][j]);
  }

  gsl_rng_free(r);

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int gain_sat(PARAMS par, float **im)
{
      int i, j,		/* loop numerators		*/
	  nsat;		/* number of saturated pixels	*/

  nsat=0;
  for (i=0; i<par.naxis2; i++)
  {
    for (j=0; j<par.naxis1; j++)
    {
      im[i][j]/=par.gain;
      if (im[i][j] > par.sat) { im[i][j]=par.sat; nsat++; }
    }
  }
  
  if (par.verbose > 0) printf("  %d saturated pixels\n", nsat);

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int new_header(PARAMS par, time_t rawtime1, time_t rawtime2, char *name,
  int *ncards, char ***header)
{
	char  **lheader,     /* header               */
	      line[CARD_SIZE];  /* header card          */
	int   i;
	struct tm *t_info;

  rawtime2-=rawtime1;
  t_info=gmtime(&rawtime1);

  if ((lheader=(char **)calloc(RECORD_CARDS, sizeof(char *))) == NULL)
  {
    perror("\n\tERROR! new_header(): calloc(lheader)");
    return(EXIT_FAILURE);
  }

  for (i=0; i<RECORD_CARDS; i++)
  {
    if ((lheader[i]=(char *)calloc(CARD_SIZE, sizeof(char))) == NULL)
    {
      perror("\n\tERROR! new_header(): calloc(lheader[i])");
      return(EXIT_FAILURE);
    }
    memset(lheader[i], ' ', CARD_SIZE);
  }

  sprintf(line, "SIMPLE  =                    T");
  memcpy(lheader[0], line, strlen(line));

  sprintf(line, "BITPIX  =                  -32 /");
  memcpy(lheader[1], line, strlen(line));

  sprintf(line, "NAXIS   =                    2 /");
  memcpy(lheader[2], line, strlen(line));

  sprintf(line, "NAXIS1  =               %6d /", par.naxis1);
  memcpy(lheader[3], line, strlen(line));

  sprintf(line, "NAXIS2  =               %6d /", par.naxis2);
  memcpy(lheader[4], line, strlen(line));

  sprintf(line, "OBJECT  = '%s' /", par.object);
  memcpy(lheader[5], line, strlen(line));

  sprintf(line, "FILTER  = '%c' /", par.filter);
  memcpy(lheader[6], line, strlen(line));

  sprintf(line, "DATE-OBS= '%d-%02d-%02d' /",
    t_info->tm_year+1900, t_info->tm_mon+1, t_info->tm_mday+1);
  memcpy(lheader[7], line, strlen(line));

  sprintf(line, "UT-START= '%02d:%02d:%02d' /",
    t_info->tm_hour, t_info->tm_min, t_info->tm_sec);
  memcpy(lheader[8], line, strlen(line));

  sprintf(line, "EXPTIME = %5ld / SECONDS", rawtime2);
  memcpy(lheader[9], line, strlen(line));

  sprintf(line, "OBSERVAT= 'CAMK' /");
  memcpy(lheader[10], line, strlen(line));

  sprintf(line, "PIXSCALE= %g / [arcsec/pixel]", par.pixscale);
  memcpy(lheader[11], line, strlen(line));

  sprintf(line, "GAIN    = %g / [photons/ADU]", par.gain);
  memcpy(lheader[12], line, strlen(line));

  sprintf(line, "COMMENT = sim2obs Copyright (C) 2014, Wojtek Pych, CAMK PAN");
  memcpy(lheader[13], line, strlen(line));

  sprintf(line, "DB-NMAG = %d / magnitude column in database", par.nm);
  memcpy(lheader[14], line, strlen(line));
  
  sprintf(line, "DISTANCE= %g / [parsec]", par.dist);
  memcpy(lheader[15], line, strlen(line));

  sprintf(line, "PIXSCALE= %g / [arcseconds]", par.pixscale);
  memcpy(lheader[16], line, strlen(line));

  sprintf(line, "SATLEVEL= %g / saturation level [ADU]", par.sat);
  memcpy(lheader[17], line, strlen(line));

  sprintf(line, "EXPOSURE= %g / exposure", par.exp);
  memcpy(lheader[18], line, strlen(line));

  sprintf(line, "SEEING  = %g / seeing [arcseconds]", par.seeing);
  memcpy(lheader[19], line, strlen(line));

  sprintf(line, "BACKGROU= %g / background level", par.bkg);
  memcpy(lheader[20], line, strlen(line));

  sprintf(line, "PSF     = %c / G-Gaussian, M-Moffat", par.psf);
  memcpy(lheader[21], line, strlen(line));

  sprintf(line, "M_BETA  = %g / Moffat beta", par.beta);
  memcpy(lheader[22], line, strlen(line));

  sprintf(line, "NOISE    = %d / 1=Poisson", par.noise);
  memcpy(lheader[23], line, strlen(line));

  sprintf(line, "RA_OFFSET= %g / R.A. offset [arcseconds]", par.raoff);
  memcpy(lheader[24], line, strlen(line));

  sprintf(line, "DECOFFSET= %g / Dec. offset [arcseconds]", par.decoff);
  memcpy(lheader[25], line, strlen(line));

  sprintf(line, "DBNAME   = %s / Input database", name);
  memcpy(lheader[26], line, strlen(line));

  sprintf(line, "END       ");
  memcpy(lheader[RECORD_CARDS-1], line, strlen(line));

  *ncards=RECORD_CARDS;
  *header=lheader;

  return(EXIT_SUCCESS);
}
/*--------------------------------------------------------*/
int main(int argc, char *argv[])
{
        char    *dbname,
                *pname,
                *oname,
		**header;
        int     ndb,            /* number of database records   */
		i,
		ncards;
        float   *x, *y, *m,
		**im;
        PARAMS  par;            /* paremeters           */
        time_t  rawtime1, rawtime2;

  time(&rawtime1);

  if (argc != 4)
  {
    printf("\n\tUSAGE: sim2obs  parameters database output_FITS\n");
    return(EXIT_FAILURE);
  }

  pname=argv[1];
  dbname=argv[2];
  oname=argv[3];

/* read parameters */
  if (read_params(pname, &par) != EXIT_SUCCESS) return(EXIT_FAILURE);

/* read the database */
  if ((ndb=read_db(par, dbname, &x, &y, &m)) < 1) return(EXIT_FAILURE);

/* calculate X,Y coordinates and magnitudes for the stars */
  if (proc_coomag(par, ndb, x, y, m) != EXIT_SUCCESS) return(EXIT_FAILURE);

/* calculate photon flux in the image pixels using PSF model */
  if (pixels(par, ndb, x, y, m, &im) != EXIT_SUCCESS) return(EXIT_FAILURE);

/* add noise */
  if (par.noise == 1)
    if (add_noise(par, rawtime1, im) != EXIT_SUCCESS) return(EXIT_FAILURE);

/* apply camera gain and saturation */
  if (gain_sat(par, im) != EXIT_SUCCESS) return(EXIT_FAILURE);

/* create the FITS header */
  time(&rawtime2);

  if (new_header(par, rawtime1, rawtime2, dbname, &ncards, &header) != EXIT_SUCCESS)
    return(EXIT_FAILURE);

/* write output into FITS file */
  write_FITS_2Dfile(oname, ncards, header, par.naxis1, par.naxis2, sizeof(float), (void **)im);

  for (i=0; i<par.naxis2; i++) free(im[i]);
  free(im);
  for (i=0; i<ncards; i++) free(header[i]);
  free(header);

  return(EXIT_SUCCESS);
}
/*** END ***/
