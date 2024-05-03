/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleGasDynamicsElemFluxAlgorithm_h
#define AssembleGasDynamicsElemFluxAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace sierra{
namespace nalu{

class Realm;
class AssembleGasDynamicsElemFluxAlgorithm : public Algorithm
{
public:

  AssembleGasDynamicsElemFluxAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *density,
    VectorFieldType *momentum,
    VectorFieldType *velocity,
    ScalarFieldType *totalH,
    ScalarFieldType *pressure,
    ScalarFieldType *temperature_,
    ScalarFieldType *speedOfSound,
    ScalarFieldType *viscosity,
    ScalarFieldType *thermalCond,
    GenericFieldType *rhsGasDyn);
  virtual ~AssembleGasDynamicsElemFluxAlgorithm() {}

  virtual void execute();

  ScalarFieldType *density_;
  VectorFieldType *momentum_;
  VectorFieldType *velocity_;
  ScalarFieldType *totalH_;
  ScalarFieldType *pressure_;
  ScalarFieldType *temperature_;
  ScalarFieldType *speedOfSound_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *thermalCond_;
  GenericFieldType *rhsGasDyn_;
  VectorFieldType *velocityRTM_;
  VectorFieldType *edgeAreaVec_;
  VectorFieldType *coordinates_;
};

} // namespace nalu
} // namespace Sierra

#endif
