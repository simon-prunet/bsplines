/*****************************************************************************
 *	Date: January 5, 2009
 *----------------------------------------------------------------------------
 *	This C program is based on the following paper:
 *		P. Thevenaz, T. Blu, M. Unser, "Interpolation Revisited,"
 *		IEEE Transactions on Medical Imaging,
 *		vol. 19, no. 7, pp. 739-758, July 2000.
 *----------------------------------------------------------------------------
 *	Philippe Thevenaz
 *	EPFL/STI/IMT/LIB/BM.4.137
 *	Station 17
 *	CH-1015 Lausanne VD
 *----------------------------------------------------------------------------
 *	phone (CET):	+41(21)693.51.61
 *	fax:			+41(21)693.37.01
 *	RFC-822:		philippe.thevenaz@epfl.ch
 *	X-400:			/C=ch/A=400net/P=switch/O=epfl/S=thevenaz/G=philippe/
 *	URL:			http://bigwww.epfl.ch/
 *----------------------------------------------------------------------------
 *	This file is best viewed with 4-space tabs (the bars below should be aligned)
 *	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|	|
 *  |...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|...|
 ****************************************************************************/

/*--------------------------------------------------------------------------*/
extern double InterpolatedValue
  (
    double  *Bcoeff,  /* input B-spline array of coefficients */
    long    Width,    /* width of the image */
    long    Height,   /* height of the image */
    double  x,    /* x coordinate where to interpolate */
    double  y,    /* y coordinate where to interpolate */
    long    SplineDegree,   /* degree of the spline model */
    long    BoundariesType /* mirror or periodic */
  );

extern double InterpolatedValue1D
  (
    double  *Bcoeff,  /* input B-spline array of coefficients */
    long    Datasize, /* size of the array */
    double  x,    /* x coordinate where to interpolate */
    long    SplineDegree,   /* degree of the spline model */
    long    BoundariesType  /* mirror or periodic */
  );

extern void ArrayInterpolatedValue
  (
    double  *Bcoeff,  /* input B-spline array of coefficients */
    long    Width,    /* width of the image */
    long    Height,   /* height of the image */
    double  *x,    /* x coordinates where to interpolate */
    double  *y,    /* y coordinates where to interpolate */
    double  *Result, /* Interpolated values, output */
    long    Npoints, /* Number of interpolated values */
    long    SplineDegree,   /* degree of the spline model */
    long    BoundariesType
  );

extern void ArrayInterpolatedValue1D
  (
    double  *Bcoeff,  /* input B-spline array of coefficients */
    long    Datasize, /* size of the array */
    double  *x,    /* x coordinates where to interpolate */
    double  *Result, /* Interpolated values, output */
    long    Npoints, /* number of interpolated values */
    long    SplineDegree,   /* degree of the spline model */
    long    BoundariesType  /* mirror or periodic */
  );

extern void ScatterValue1D
  (
    double  value, /* value at position x, to be scattered on spline coefficients */
    double  x,    /* x coordinate of value */
    double  *Bcoeff,  /* input B-spline array of coefficients, to be updated */
    long    Datasize, /* size of the array */
    long    SplineDegree,   /* degree of the spline model */
    long    BoundariesType  /* mirror or periodic */
  );

extern void ArrayScatterValue1D
  (
    double  *value, /* values at position x, to be scattered on spline coefficients */
    double  *x,    /* x coordinates of values */
    long    Npoints, /* number of values to be scattered */
    double  *Bcoeff,  /* input B-spline array of coefficients, to be updated */
    long    Datasize, /* size of the array */
    long    SplineDegree,   /* degree of the spline model */
    long    BoundariesType  /* mirror or periodic */
  );

extern int SplineWeights
  (
    double x,         /* x position */
    double y,         /* y position */
    double *xWeight, /* Spline weights, x dimension */
    double *yWeight, /* Spline weights, y dimension */
    long *xIndex,  /* Spline index, x dimension */
    long *yIndex,  /* Spline index, y dimension */
    long SplineDegree /* degree of spline model */
  );


extern int SplineWeights1D
  (
    double x, /* x position */
    double *xWeight, /* Spline weights, x dimension */
    long *xIndex, /* Spline index, x dimension */
    long SplineDegree /* degree of spline model */
  );

