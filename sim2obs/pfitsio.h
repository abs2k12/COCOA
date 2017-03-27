/*========================================================*/
/*                                                        */
/*  pfitsio.h       version 3.5.0         2015.10.07      */
/*                                                        */
/*  Copyright (C) 2007-2015 by Wojtek Pych, CAMK PAN      */
/*                                                        */
/*  Written for GNU project C and C++ Compiler            */
/*  Private FITS I/O library.                             */
/*                                                        */
/*========================================================*/

#include "pfitsin.h"

extern double scale(const int, float *, int16_t *, double *);
extern int write_FITS_1Dimage(FILE *, const int, const int, void *);
extern int write_FITS_2Dimage(FILE *, const int, const int, const int,
  void **);
extern int write_FITS_1Dfile(char *, const int, char **, const int, const int,
  void *);
extern int write_FITS_2Dfile(char *, const int, char **, const int, const int,
  const int, void **);

/*** END ***/
