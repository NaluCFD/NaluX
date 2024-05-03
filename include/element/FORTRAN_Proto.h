/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/



#ifndef FORTRAN_Proto_h
#define FORTRAN_Proto_h

#include <stk_util/util/Fortran.hpp>

extern "C" void
SIERRA_FORTRAN( hex_scs_det ) ( const int*  nelem,
                                const int*  npe,
                                const int*  nscs,
                                const double* coords,
                                double* area_vec);

extern "C" void
SIERRA_FORTRAN( hex_scv_det ) ( const int*  nelem,
                                const int*  npe,
                                const int*  nscv,
                                const double* coords,
                                double* volume,
                                double* error,
                                int *nerr);

extern "C" void
SIERRA_FORTRAN( quad3d_scs_det ) ( const int*  nelem,
                                   const double* coords,
                                   double* areav );

extern "C" void
SIERRA_FORTRAN( hex_shape_fcn ) ( const int*  npts,
                                  const double* par_coords,
                                  double* shape_fcn );

extern "C" void
SIERRA_FORTRAN( hex_derivative ) ( const int*  npts,
                                   const double* par_coords,
                                   double* deriv );
extern "C" void
SIERRA_FORTRAN( quad_derivative ) ( const int*  npts,
				    const double* par_coords,
				    double* deriv );
				 
extern "C" void
SIERRA_FORTRAN( hex_gradient_operator ) ( const int*  nelem,
                                          const int* npe,
                                          const int* numint,
                                          double *deriv,
                                          const double* coords,
                                          double* gradop,
                                          double *det_j,
                                          double *error,
                                          int* lerr);

extern "C" void
SIERRA_FORTRAN( quad3d_shape_fcn ) ( const int*  npts,
                                     const double* par_coords,
                                     double* shape_fcn );

#endif // FORTRAN_Proto_h
