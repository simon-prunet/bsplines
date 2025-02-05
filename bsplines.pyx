import numpy as nm
cimport numpy as nm
nm.import_array()
from libc.stdlib cimport *
from libc.stdio cimport *
from libc.string cimport *
from libc.math cimport *

################################################
# Spline interpolation routines ################
################################################

cdef extern from "coeff.h":
  int C_SamplesToCoefficients1D "SamplesToCoefficients1D" (double*,long,long,long)
  int C_SamplesToCoefficients2D "SamplesToCoefficients" (double*,long,long,long,long)
cdef extern from "interpol.h":
  double C_InterpolatedValue1D "InterpolatedValue1D" (double*,long,double,long,long)
  double C_InterpolatedValue2D "InterpolatedValue" (double*,long,long,double,double,long,long)
  void C_ArrayInterpolatedValue1D "ArrayInterpolatedValue1D" (double*,long,double*,double*,long,long,long)
  void C_ArrayInterpolatedValue2D "ArrayInterpolatedValue" (double*,long,long,double*,double*,double*,long,long,long)
  void C_ScatterValue1D "ScatterValue1D" (double,double,double*,long,long,long)
  void C_ArrayScatterValue1D "ArrayScatterValue1D" (double*,double*,long,double*,long,long,long)
  void C_SplineWeights1D "SplineWeights1D" (double*, double*, long)
  void C_SplineWeights2D "SplineWeights2D" (double*, double*, double*, double*, long)

def SamplesToCoefficients1D (nm.ndarray[nm.double_t,ndim=1] Data,
                             long SplineDegree, long BoundariesType):
  '''
  Takes data on a *regular* grid and computes B-spline coefficients
  '''
  cdef int err
  cdef long Datasize
  Datasize = Data.size
  err = C_SamplesToCoefficients1D(<double*>Data.data, <long> Datasize, <long> SplineDegree,
                                  <long> BoundariesType)
  return err

def InterpolatedValue1D (nm.ndarray[nm.double_t,ndim=1] Bcoeff,
                         x, long SplineDegree, long BoundariesType):

  '''
  Takes B-spline coefficient array and scalar (or vector) coordinate(s) x, 
  return value(s) at the coordinate(s) x
  '''

  cdef long Datasize
  cdef nm.ndarray[nm.double_t,ndim=1] xx, Result
  Datasize = Bcoeff.size
  if isinstance(x,(int,float)): 
    y = C_InterpolatedValue1D(<double*>Bcoeff.data, <long> Datasize, <double> x,
                            <long> SplineDegree, <long> BoundariesType)
    return y
  else: #Array
    Result=nm.zeros(x.size,dtype=nm.double)
    xx=x
    C_ArrayInterpolatedValue1D(<double*>Bcoeff.data, <long> Datasize, <double*> xx.data, <double*> Result.data,
                               <long> x.size, <long> SplineDegree, <long> BoundariesType)
    return Result

def ScatterValue1D (value,x,nm.ndarray[nm.double_t,ndim=1] Bcoeff, 
                    long SplineDegree, long BoundariesType):

  '''
  Adjoint (not inverse !) operation of InterpolatedValue1D. Scatters value(s) at coordinate(s) x
  onto B-spline coefficients.
  '''

  cdef long Datasize
  cdef nm.ndarray[nm.double_t,ndim=1] xx, vv
  Datasize = Bcoeff.size
  if isinstance(x,(int,float)):
    C_ScatterValue1D(<double> value, <double> x, <double*> Bcoeff.data, 
                     <long> Datasize, <long> SplineDegree, <long> BoundariesType)
  else: # Array
    xx = x
    vv = value
    C_ArrayScatterValue1D(<double*>vv.data, <double*> xx.data, <long> xx.size, 
                          <double*> Bcoeff.data, <long> Datasize, <long> SplineDegree, <long> BoundariesType)
  return

def CoefficientsToCoefficients1D (nm.ndarray[nm.double_t,ndim=1] Bcoeff, nm.ndarray[nm.double_t,ndim=1] Samples,
                                  long SplineDegree, long BoundariesType):

  '''
  Applies in a sequence InterpolatedValue1D and ScatterValue1D. This will be useful for the inverse
  problem of estimating B-spline coefficients from irregular samples and values.
  Takes B-spline coefficients and Samples as input, returns B-spline coefficients.
  Typically used with more samples than coefficients for the inverse problem to be well posed.
  '''

  cdef nm.ndarray[nm.double_t,ndim=1] SampleValues, CoeffResult
  SampleValues=nm.zeros(Samples.size)
  CoeffResult=nm.zeros(Bcoeff.size)

  cdef long Bcoeff_size, Samples_size
  Bcoeff_size = Bcoeff.size
  Samples_size = Samples.size
  C_ArrayInterpolatedValue1D(<double*>Bcoeff.data, <long>Bcoeff_size, <double*> Samples.data, <double*>SampleValues.data,
                             <long>Samples_size, <long>SplineDegree, <long>BoundariesType)
  C_ArrayScatterValue1D(<double*>SampleValues.data,<double*>Samples.data,<long>Samples_size,
                        <double*>CoeffResult.data, <long>Bcoeff_size, <long> SplineDegree, <long> BoundariesType)

  return CoeffResult
  #####
  
def EstimateCoefficients1D (Samples, Values, n_coeffs, SplineDegree, BoundariesType, epsilon=1e-6):

  '''
  Solves the inverse problem of estimating n_coeffs B-spline coefficients from irregularly spaced 
  coordinates (Samples) and associated values (Values).
  NB: Samples are supposed to be scaled so that they pertain to the [0, n_coeffs] interval
  Returns array of B-spline coefficients.
  '''
  # Compute RHS (A^T.d)
  Coeffs_rhs = nm.zeros(n_coeffs)
  ScatterValue1D(Values,Samples,Coeffs_rhs,SplineDegree,BoundariesType)
  # Define A^T.A
  def Amul(coeffs):
    return CoefficientsToCoefficients1D(coeffs,Samples,SplineDegree,BoundariesType)
  # Solve (A^T.A)^{-1}.A^T.d
  res = nm.zeros(n_coeffs)
  res = cg_solve(Coeffs_rhs, Amul, res, epsilon=epsilon)
  return(res)

def SamplesToCoefficients2D (nm.ndarray[nm.double_t,ndim=2] Data, long Width, long Height,
                             long SplineDegree, long BoundariesType):

  ## Beware that the 2D size is (Height, Width in that order, i.e. (y,x)
  cdef nm.ndarray[nm.double_t,ndim=1] FlatData
  FlatData = Data.reshape(-1,) # Flattened view of the array
  err = C_SamplesToCoefficients2D(<double*>FlatData.data,<long> Width, <long> Height,
                                  <long> SplineDegree, <long> BoundariesType)
  return err

def InterpolatedValue2D(nm.ndarray[nm.double_t,ndim=2] Bcoeff, long Width, long Height,
                        x, y, long SplineDegree, long BoundariesType):
  
  cdef nm.ndarray[nm.double_t,ndim=1] FlatData, xx, yy, Result
  cdef double res

  if isinstance(x,(int,float)):
    FlatData = Bcoeff.reshape(-1,) ## Flattened view
    res = C_InterpolatedValue2D(<double*>FlatData.data,<long> Width, <long> Height,
                                <double> x, <double> y, <long> SplineDegree, <long> BoundariesType)
    return res
  else: # Array
    Result=nm.zeros(x.size,dtype=nm.double)
    FlatData = Bcoeff.reshape(-1)
    xx = x
    yy = y
    C_ArrayInterpolatedValue2D(<double*>FlatData.data, <long> Width, <long> Height, 
                               <double*> xx.data, <double*> yy.data,
                               <double*> Result.data, <long> x.size, <long> SplineDegree, 
                               <long> BoundariesType)
    return Result

## Simple resampling example
def ShiftVector (nm.ndarray[nm.double_t,ndim=1] Vector, shift=0.5, SplDeg=2, BndType=0):

  tmp = Vector.copy()
  si = tmp.shape[0]
  err = SamplesToCoefficients1D(tmp,si,SplDeg,BndType)
  x = nm.arange(si)
  result = InterpolatedValue1D(tmp,si,x+shift,SplDeg,BndType)
  return(result)

def ShiftImage (nm.ndarray[nm.double_t,ndim=2] Image, shift=(0.5,0.5), SplDeg=2, BndType=0):

  tmp = Image.copy()
  Width,Height = tmp.shape
  err = SamplesToCoefficients2D(tmp,Width,Height,SplDeg,BndType)
  x=nm.outer(nm.ones(Height,dtype=nm.double),nm.arange(Width))
  y=nm.outer(nm.arange(Height),nm.ones(Width,dtype=nm.double))
  result = InterpolatedValue2D(tmp,Width,Height,x.reshape(-1)+shift[0],y.reshape(-1)+shift[1],SplDeg,BndType)

  return nm.reshape(result,(Width,Height))

## Exemple image
def get_lena(path='.'):
  filename=path+'/lena.img'
  res=nm.fromfile(filename,dtype=nm.uint8)
  return nm.reshape(res,(256,256))*1.0

def get_tiger(path='.'):
  import astropy.io.fits as fits
  filename=path+'/graytiger.fits'
  res=fits.getdata(filename).astype(nm.double)
  return res

## PCG solver part

def cg_solve (b,Amul,x0,epsilon=1e-6,nstepmax=500,Mdiv=None):

  x=x0
  r = b - Amul(x)
  if (Mdiv is not None):
    z=Mdiv(r)
  else:
    z=r
  p=z*1.
  rz = nm.dot(r,z)

  xold = x*1.
  rold = r*1.
  rzold = rz*1.

  nn=0
  while nn<nstepmax:
    Ap = Amul(p)
    if (rz != 0):
      alpha = rz / nm.dot(p,Ap)
    else:
      alpha = 0.0
    xold = x
    rold = r
    rzold = rz
    x = x + alpha*p
    r = r - alpha*Ap

    print(nn,nm.sqrt(nm.dot(r,r)))
    
    if cg_crit(r,epsilon):
      return x
    if (Mdiv is not None):
      z = Mdiv(r)
    else:
      z=r
    rz = nm.dot(r,z)
    if (rz != 0):
      beta = rz / rzold
    else:
      beta = 0.
    p = z + beta*p
    nn += 1

  print ('CG did not converge in %d iterations '%nstepmax)
  return

def cg_crit(r,epsilon):

  normr = nm.sqrt(nm.dot(r,r))
  if (normr < epsilon):
    return True
  else:
    return False


