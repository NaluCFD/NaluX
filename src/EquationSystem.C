/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <EquationSystem.h>
#include <AuxFunctionAlgorithm.h>
#include <AlgorithmDriver.h>
#include <InitialConditions.h>
#include <Realm.h>
#include <Simulation.h>
#include <SolutionOptions.h>
#include <FieldTypeDef.h>
#include <NaluParsing.h>
#include <NaluEnv.h>
#include <ConstantAuxFunction.h>
#include <Enums.h>


// stk
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

EquationSystem::EquationSystem(
  EquationSystems& eqSystems,
  const std::string name,
  const std::string eqnTypeName)
  : equationSystems_(eqSystems),
    realm_(eqSystems.realm_),
    name_(name),
    userSuppliedName_(name),
    eqnTypeName_(eqnTypeName),
    maxIterations_(1),
    convergenceTolerance_(1.0),
    timerAssemble_(0.0),
    timerLoadComplete_(0.0),
    timerSolve_(0.0),
    timerMisc_(0.0),
    timerInit_(0.0),
    timerPrecond_(0.0),
    avgLinearIterations_(0.0),
    maxLinearIterations_(0.0),
    minLinearIterations_(1.0e10),
    nonLinearIterationCount_(0),
    reportLinearIterations_(false),
    firstTimeStepSolve_(true)
{
  // nothing to do
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
EquationSystem::~EquationSystem()
{
  // initial conditions and bc prop
  std::vector<AuxFunctionAlgorithm *>::iterator ii;
  for( ii=bcDataAlg_.begin(); ii!=bcDataAlg_.end(); ++ii )
    delete *ii;

  std::vector<Algorithm *>::iterator iim;
  for( iim=bcDataMapAlg_.begin(); iim!=bcDataMapAlg_.end(); ++iim )
    delete *iim;

  for( iim=copyStateAlg_.begin(); iim!=copyStateAlg_.end(); ++iim )
    delete *iim;

  // prop algs
  std::vector<Algorithm *>::iterator iip;
  for( iip=propertyAlg_.begin(); iip!=propertyAlg_.end(); ++iip ) {
    delete *iip;
  }

  for (auto it = preIterAlgDriver_.begin(); it != preIterAlgDriver_.end(); ++it) {
    delete *it;
  }

  for (auto it = postIterAlgDriver_.begin(); it != postIterAlgDriver_.end(); ++it) {
    delete *it;
  }
}

//--------------------------------------------------------------------------
//-------- system_is_converged ---------------------------------------------
//--------------------------------------------------------------------------
bool
EquationSystem::system_is_converged()
{
  bool isConverged = true;
  // check a linear system
  return isConverged;
}

//--------------------------------------------------------------------------
//-------- pre_timestep_work -----------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::pre_timestep_work()
{
  firstTimeStepSolve_ = true;
}

//--------------------------------------------------------------------------
//-------- provide_scaled_norm ---------------------------------------------
//--------------------------------------------------------------------------
double
EquationSystem::provide_scaled_norm()
{
  // query linsys_
  return ( 0.0 );
}

//--------------------------------------------------------------------------
//-------- provide_norm ---------------------------------------------
//--------------------------------------------------------------------------
double
EquationSystem::provide_norm()
{
  // query linsys_
  return ( 0.0 );
}

//--------------------------------------------------------------------------
//-------- provide_norm_increment ------------------------------------------
//--------------------------------------------------------------------------
double
EquationSystem::provide_norm_increment()
{
  // query linsys_
  return ( 0.0 );
}

//--------------------------------------------------------------------------
//-------- dump_eq_time ----------------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::dump_eq_time()
{
  // subtract out preconditioning time from solve time
  timerSolve_ -= timerPrecond_;

  double l_timer[6] = {timerAssemble_, timerLoadComplete_, timerSolve_, timerMisc_, timerInit_, timerPrecond_};
  double g_min[6] = {};
  double g_max[6] = {};
  double g_sum[6] = {};

  int nprocs = NaluEnv::self().parallel_size();

  NaluEnv::self().naluOutputP0() << "Timing for Eq: " << userSuppliedName_ << std::endl;

  // get max, min, and sum over processes
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &l_timer[0], &g_sum[0], 6);
  stk::all_reduce_min(NaluEnv::self().parallel_comm(), &l_timer[0], &g_min[0], 6);
  stk::all_reduce_max(NaluEnv::self().parallel_comm(), &l_timer[0], &g_max[0], 6);

  // output
  NaluEnv::self().naluOutputP0() << "             init --  " << " \tavg: " << g_sum[4]/double(nprocs)
                  << " \tmin: " << g_min[4] << " \tmax: " << g_max[4] << std::endl;
  NaluEnv::self().naluOutputP0() << "         assemble --  " << " \tavg: " << g_sum[0]/double(nprocs)
                  << " \tmin: " << g_min[0] << " \tmax: " << g_max[0] << std::endl;
  NaluEnv::self().naluOutputP0() << "    load_complete --  " << " \tavg: " << g_sum[1]/double(nprocs)
                  << " \tmin: " << g_min[1] << " \tmax: " << g_max[1] << std::endl;
  NaluEnv::self().naluOutputP0() << "            solve --  " << " \tavg: " << g_sum[2]/double(nprocs)
                  << " \tmin: " << g_min[2] << " \tmax: " << g_max[2] << std::endl;
  NaluEnv::self().naluOutputP0() << "    precond setup --  " << " \tavg: " << g_sum[5]/double(nprocs)
                  << " \tmin: " << g_min[5] << " \tmax: " << g_max[5] << std::endl;
  NaluEnv::self().naluOutputP0() << "             misc --  " << " \tavg: " << g_sum[3]/double(nprocs)
                  << " \tmin: " << g_min[3] << " \tmax: " << g_max[3] << std::endl;

  if (reportLinearIterations_)
    NaluEnv::self().naluOutputP0() << "linear iterations -- " << " \tavg: " << avgLinearIterations_
                    << " \tmin: " << minLinearIterations_ << " \tmax: "
                    << maxLinearIterations_ << std::endl;

  // reset anytime these are called; 
  // some EquationSystems have no linear system, e.g., LowMach holds .. uvw_p
  timerAssemble_ = 0.0;
  timerLoadComplete_ = 0.0;
  timerMisc_ = 0.0;
  timerSolve_ = 0.0;
  timerInit_ = 0.0;
  timerPrecond_ = 0.0;
  /*
  if ( NULL != linsys_ )
    linsys_->zero_timer_precond();
  */
  avgLinearIterations_ = 0.0;
  minLinearIterations_ = 1.0e10;
  maxLinearIterations_ = 0.0;
  nonLinearIterationCount_ = 0;
}

//--------------------------------------------------------------------------
//-------- update_iteration_statistics -------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::update_iteration_statistics(
  const int & iters)
{
  const double iterations = (double)iters;
  avgLinearIterations_ = (nonLinearIterationCount_*avgLinearIterations_
                          + iterations)/(nonLinearIterationCount_+1);
  maxLinearIterations_ = std::max(maxLinearIterations_,iterations);
  minLinearIterations_ = std::min(minLinearIterations_,iterations);
  nonLinearIterationCount_ += 1;
  reportLinearIterations_ = true;
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::initial_work()
{
  std::vector<Algorithm *>::iterator ii;
  for( ii=copyStateAlg_.begin(); ii!=copyStateAlg_.end(); ++ii )
    (*ii)->execute();
}

//--------------------------------------------------------------------------
//-------- assemble_and_solve ----------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::assemble_and_solve(
  stk::mesh::FieldBase *deltaSolution)
{
  // TBD
}

//--------------------------------------------------------------------------
//-------- bc_data_specified ----------------------------------------------------
//--------------------------------------------------------------------------
bool
EquationSystem::bc_data_specified(
  const UserData &userData, std::string &name)
{
  bool isSpecified = false;
  std::map<std::string, bool>::const_iterator iter
    = userData.bcDataSpecifiedMap_.find(name);
  if ( iter != userData.bcDataSpecifiedMap_.end() ) {
    isSpecified = true;
  }
  return isSpecified;
}

//--------------------------------------------------------------------------
//-------- get_bc_data_type ------------------------------------------------
//--------------------------------------------------------------------------
UserDataType
EquationSystem::get_bc_data_type(
  const UserData &userData, std::string &name)
{
  UserDataType dataType = CONSTANT_UD;
  std::map<std::string, UserDataType>::const_iterator iter
    = userData.bcDataTypeMap_.find(name);
  if ( iter != userData.bcDataTypeMap_.end() ) {
    dataType = (*iter).second;
  }
  return dataType;
}

//--------------------------------------------------------------------------
//-------- get_bc_function_name --------------------------------------------
//--------------------------------------------------------------------------
std::string
EquationSystem::get_bc_function_name(
  const UserData &userData, std::string &name)
{
  std::string fcnName = "no_name";
  std::map<std::string, std::string>::const_iterator iter
    = userData.userFunctionMap_.find(name);
  if ( iter != userData.userFunctionMap_.end() ) {
    fcnName = (*iter).second;
  }
  else {
    throw std::runtime_error("wall function name not found");
  }
  
  return fcnName;
}

//--------------------------------------------------------------------------
//-------- get_bc_function_params ------------------------------------------
//--------------------------------------------------------------------------
std::vector<double>
EquationSystem::get_bc_function_params(
  const UserData &userData, std::string &name)
{
  std::vector<double> theParams;
  std::map<std::string, std::vector<double> >::const_iterator iter
    = userData.functionParams_.find(name);
  if ( iter != userData.functionParams_.end() ) {
    return (*iter).second;
  }
  else {
    return theParams;
  }
}

//--------------------------------------------------------------------------
//-------- get_bc_function_string_params -----------------------------------
//--------------------------------------------------------------------------
std::vector<std::string>
EquationSystem::get_bc_function_string_params(
  const UserData &userData, std::string &name)
{
  std::vector<std::string> theParams;
  std::map<std::string, std::vector<std::string> >::const_iterator iter
    = userData.functionStringParams_.find(name);
  if ( iter != userData.functionStringParams_.end() ) {
    return (*iter).second;
  }
  else {
    return theParams;
  }
}

//--------------------------------------------------------------------------
//-------- evaluate_properties ---------------------------------------------
//--------------------------------------------------------------------------
void
EquationSystem::evaluate_properties()
{
  for ( size_t k = 0; k < propertyAlg_.size(); ++k ) {
    propertyAlg_[k]->execute();
  }
}

void
EquationSystem::pre_iter_work()
{
  for (auto it: preIterAlgDriver_) {
    it->execute();
  }
}

void
EquationSystem::post_iter_work()
{
  for (auto it: postIterAlgDriver_) {
    it->execute();
  }
}

} // namespace nalu
} // namespace Sierra
