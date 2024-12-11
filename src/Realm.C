/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "Realm.h"
#include "Simulation.h"
#include "NaluEnv.h"

#include "AuxFunction.h"
#include "AuxFunctionAlgorithm.h"
#include "ComputeGeometryAlgorithmDriver.h"
#include "ComputeGeometryInteriorAlgorithm.h"
#include "ConstantAuxFunction.h"
#include "Enums.h"
#include "EntityExposedFaceSorter.h"
#include "EquationSystem.h"
#include "EquationSystems.h"
#include "FieldTypeDef.h"
#include "element/Element.h"
#include "MaterialProperty.h"
#include "MaterialPropertys.h"
#include "NaluParsing.h"
#include "OutputInfo.h"
#include "PeriodicManager.h"
#include "Realms.h"
#include "SolutionOptions.h"
#include "TimeIntegrator.h"

// props; algs, evaluators and data
#include "property_evaluator/ConstantPropertyEvaluator.h"
#include "property_evaluator/EnthalpyPropertyEvaluator.h"
#include "property_evaluator/IdealGasPropertyEvaluator.h"
#include "property_evaluator/MaterialPropertyData.h"
#include "property_evaluator/ReferencePropertyData.h"

// transfer
#include "xfer/Transfer.h"

// stk_util
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/perf_util.hpp>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/CreateEdges.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_mesh/base/FieldBLAS.hpp>

// stk_io
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_io/InputFile.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

// yaml for parsing..
#include <yaml-cpp/yaml.h>
#include <NaluParsing.h>
#include <NaluParsingHelper.h>

// basic c++
#include <map>
#include <cmath>
#include <utility>
#include <stdint.h>

#define USE_NALU_PERFORMANCE_TESTING_CALLGRIND 0
#if USE_NALU_PERFORMANCE_TESTING_CALLGRIND
#include "/usr/netpub/valgrind-3.8.1/include/valgrind/callgrind.h"
#endif

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// Realm - do some stuff
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
  Realm::Realm(Realms& realms, const YAML::Node & node)
  : realms_(realms),
    name_("na"),
    type_("multi_physics"),
    inputDBName_("input_unknown"),
    spatialDimension_(3u),  // for convenience; can always get it from meta data
    realmUsesEdges_(false),
    solveFrequency_(1),
    isTurbulent_(false),
    needsEnthalpy_(false),
    l2Scaling_(1.0),
    bulkData_(NULL),
    ioBroker_(NULL),
    resultsFileIndex_(99),
    restartFileIndex_(99),
    computeGeometryAlgDriver_(0),
    numInitialElements_(0),
    timeIntegrator_(0),
    boundaryConditions_(*this),
    initialConditions_(*this),
    materialPropertys_(*this),
    equationSystems_(*this),
    maxCourant_(0.0),
    maxReynolds_(0.0),
    targetCourant_(1.0),
    timeStepChangeFactor_(1.25),
    currentNonlinearIteration_(1),
    solutionOptions_(new SolutionOptions()),
    outputInfo_(new OutputInfo()),
    timerCreateMesh_(0.0),
    timerPopulateMesh_(0.0),
    timerPopulateFieldData_(0.0),
    timerOutputFields_(0.0),
    timerNonconformal_(0.0),
    timerInitializeEqs_(0.0),
    timerPropertyEval_(0.0),
    timerTransferSearch_(0.0),
    timerTransferExecute_(0.0),
    timerSkinMesh_(0.0),
    timerSortExposedFace_(0.0),
    hasNonConformal_(false),
    hasOverset_(false),
    hasMultiPhysicsTransfer_(false),
    hasInitializationTransfer_(false),
    hasIoTransfer_(false),
    hasExternalDataTransfer_(false),
    periodicManager_(NULL),
    hasPeriodic_(false),
    hasFluids_(false),
    globalParameters_(),
    exposedBoundaryPart_(0),
    checkForMissingBcs_(false),
    checkJacobians_(false),
    isothermal_(true),
    uniform_(true),
    gasDynamics_(false),
    provideEntityCount_(false),
    autoDecompType_("None"),
    activateAura_(false),
    activateMemoryDiagnostic_(false),
    supportInconsistentRestart_(false),
    doBalanceNodes_(false),
    wallTimeStart_(stk::wall_time()),
    inputMeshIdx_(-1),
    node_(node),
    usesCVFEM_(true),
    minDualVolume_(1.0)
{
  // deal with specialty options that live off of the realm; 
  // choose to do this now rather than waiting for the load stage
  look_ahead_and_creation(node);
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
Realm::~Realm()
{
  delete ioBroker_;

  if ( NULL != computeGeometryAlgDriver_ )
    delete computeGeometryAlgDriver_;

  // prop algs
  std::vector<Algorithm *>::iterator ii;
  for( ii=initCondAlg_.begin(); ii!=initCondAlg_.end(); ++ii )
    delete *ii;
  for( ii=propertyAlg_.begin(); ii!=propertyAlg_.end(); ++ii ) {
    delete *ii;
  }

  // any bc data
  std::vector<AuxFunctionAlgorithm *>::iterator iaux;
  for( iaux=bcDataAlg_.begin(); iaux!=bcDataAlg_.end(); ++iaux )
    delete *iaux;

  delete solutionOptions_;
  delete outputInfo_;
  
  ElementRepo::clear();
}

void
Realm::breadboard()
{
  if ( usesCVFEM_ )
    computeGeometryAlgDriver_ = new ComputeGeometryAlgorithmDriver(*this);
  equationSystems_.breadboard();
}

bool
Realm::debug() const
{
  return root()->debug_;
}

//--------------------------------------------------------------------------
//-------- get_activate_memory_diagnostic ----------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_activate_memory_diagnostic()
{
  return activateMemoryDiagnostic_;
}

//--------------------------------------------------------------------------
//-------- provide_memory_summary ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_memory_summary() 
{
  size_t now, hwm;
  stk::get_memory_usage(now, hwm);
  // min, max, sum
  size_t global_now[3] = {now,now,now};
  size_t global_hwm[3] = {hwm,hwm,hwm};
  
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceSum<1>( &global_now[2] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMin<1>( &global_now[0] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMax<1>( &global_now[1] ) );
  
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceSum<1>( &global_hwm[2] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMin<1>( &global_hwm[0] ) );
  stk::all_reduce(NaluEnv::self().parallel_comm(), stk::ReduceMax<1>( &global_hwm[1] ) );
  
  NaluEnv::self().naluOutputP0() << "Memory Overview: " << std::endl;
  NaluEnv::self().naluOutputP0() << "nalu memory: total (over all cores) current/high-water mark= "
                                 << std::setw(15) << convert_bytes(global_now[2])
                                 << std::setw(15) << convert_bytes(global_hwm[2])
                                 << std::endl;
  
  NaluEnv::self().naluOutputP0() << "nalu memory:   min (over all cores) current/high-water mark= "
                                 << std::setw(15) << convert_bytes(global_now[0])
                                 << std::setw(15) << convert_bytes(global_hwm[0])
                                 << std::endl;
  
  NaluEnv::self().naluOutputP0() << "nalu memory:   max (over all cores) current/high-water mark= "
                                  << std::setw(15) << convert_bytes(global_now[1])
                                  << std::setw(15) << convert_bytes(global_hwm[1])
                                  << std::endl;
}

//--------------------------------------------------------------------------
//-------- convert_bytes ---------------------------------------------------
//--------------------------------------------------------------------------
std::string 
Realm::convert_bytes(double bytes)
{
  const double K = 1024;
  const double M = K*1024;
  const double G = M*1024;

  std::ostringstream out;
  if (bytes < K) {
    out << bytes << " B";
  } else if (bytes < M) {
    bytes /= K;
    out << bytes << " K";
  } else if (bytes < G) {
    bytes /= M;
    out << bytes << " M";
  } else {
    bytes /= G;
    out << bytes << " G";
  }
  return out.str();
}

//--------------------------------------------------------------------------
//-------- initialize ------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::initialize()
{
  NaluEnv::self().naluOutputP0() << "Realm::initialize() Begin " << std::endl;

  // field registration
  setup_nodal_fields();
  setup_edge_fields();
  setup_element_fields();

  // property maps and evaluation algorithms
  setup_property();

  // interior algorithm creation
  setup_interior_algorithms();

  // create boundary conditions
  setup_bc();
  
  // post processing algorithm creation
  setup_post_processing_algorithms();

  // create initial conditions
  setup_initial_conditions();

  // set global variables that have not yet been set
  initialize_global_variables();

  // Populate_mesh fills in the entities (nodes/elements/etc) and
  // connectivities, but no field-data. Field-data is not allocated yet.
  NaluEnv::self().naluOutputP0() << "Realm::ioBroker_->populate_mesh() Begin" << std::endl;
  double time = -NaluEnv::self().nalu_time();
  ioBroker_->populate_mesh();
  time += NaluEnv::self().nalu_time();
  timerPopulateMesh_ += time;
  NaluEnv::self().naluOutputP0() << "Realm::ioBroker_->populate_mesh() End" << std::endl;

  if (doBalanceNodes_) {
    balance_nodes();
  }

  // If we want to create all internal edges, we want to do it before
  // field-data is allocated because that allows better performance in
  // the create-edges code.
  if (realmUsesEdges_ )
    throw std::runtime_error("add a create_edges");

  // create the nodes for possible data probe

  // output entity counts including max/min
  if ( provideEntityCount_ )
    provide_entity_count();

  // Now the mesh is fully populated, so we're ready to populate
  // field-data including coordinates, and attributes and/or distribution factors
  // if those exist on the input mesh file.
  NaluEnv::self().naluOutputP0() << "Realm::ioBroker_->populate_field_data() Begin" << std::endl;
  time = -NaluEnv::self().nalu_time();
  ioBroker_->populate_field_data();
  time += NaluEnv::self().nalu_time();
  timerPopulateFieldData_ += time;
  NaluEnv::self().naluOutputP0() << "Realm::ioBroker_->populate_field_data() End" << std::endl;

  // manage NaluGlobalId for linear system
  set_global_id();

  // check that all bcs are covering exposed surfaces
  if ( checkForMissingBcs_ )
    enforce_bc_on_exposed_faces();

  // output and restart files
  create_output_mesh();
  create_restart_mesh();

  // sort exposed faces only when using consolidated bc NGP approach
  if ( solutionOptions_->useConsolidatedBcSolverAlg_ ) {
    const double timeSort = NaluEnv::self().nalu_time();
    bulkData_->sort_entities(EntityExposedFaceSorter());
    timerSortExposedFace_ += (NaluEnv::self().nalu_time() - timeSort);
  }
  
  // variables that may come from the initial mesh
  input_variables_from_mesh();

  populate_boundary_data();
 
  if ( has_mesh_deformation() )
    init_current_coordinates();

  if ( hasPeriodic_ )
    periodicManager_->build_constraints();

  compute_geometry();

  compute_l2_scaling();

  equationSystems_.initialize();

  // check job run size after mesh creation, linear system initialization
  check_job(false);

  NaluEnv::self().naluOutputP0() << "Realm::initialize() End " << std::endl;
}

//--------------------------------------------------------------------------
//-------- look_ahead_and_creation -----------------------------------------
//--------------------------------------------------------------------------
void
Realm::look_ahead_and_creation(const YAML::Node & node)
{
  // nothing to do
}
  
//--------------------------------------------------------------------------
//-------- load ------------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::load(const YAML::Node & node)
{

  //======================================
  // realm commands first
  //======================================

  name_ = node["name"].as<std::string>() ;
  inputDBName_ = node["mesh"].as<std::string>() ;
  get_if_present(node, "type", type_, type_);

  // provide a high level banner
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Realm Options Review: " << name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;

  // check for FEM (not yet ready for mixed realms)
  bool usesFEM = false;
  get_if_present(node, "activate_fem", usesFEM, usesFEM);
  usesCVFEM_ = !usesFEM;

  // exposed bc check
  get_if_present(node, "check_for_missing_bcs", checkForMissingBcs_, checkForMissingBcs_);

  // check for bad Jacobians in the mesh
  get_if_present(node, "check_jacobians", checkJacobians_, checkJacobians_);

  // entity count
  get_if_present(node, "provide_entity_count", provideEntityCount_, provideEntityCount_);

  // determine if edges are required and whether or not stk handles this
  get_if_present(node, "use_edges", realmUsesEdges_, realmUsesEdges_);

  // let everyone know about core algorithm
  if ( realmUsesEdges_ ) {
    throw std::runtime_error("Edge-based scheme is not supported");
  }
  else {
    NaluEnv::self().naluOutputP0() <<"Element-based scheme will be activated" << std::endl;
  }

  // how often is the realm solved..
  get_if_present(node, "solve_frequency", solveFrequency_, solveFrequency_);

  // automatic decomposition
  get_if_present(node, "automatic_decomposition_type", autoDecompType_, autoDecompType_);
  if ( "None" != autoDecompType_ ) {
    NaluEnv::self().naluOutputP0() 
      <<"Warning: When using automatic_decomposition_type, one must have a serial file" << std::endl;
  }

  // activate aura
  get_if_present(node, "activate_aura", activateAura_, activateAura_);
  if ( activateAura_ )
    NaluEnv::self().naluOutputP0() << "Nalu will activate aura ghosting" << std::endl;
  else
    NaluEnv::self().naluOutputP0() << "Nalu will deactivate aura ghosting" << std::endl;

  // memory diagnostic
  get_if_present(node, "activate_memory_diagnostic", activateMemoryDiagnostic_, activateMemoryDiagnostic_);
  if ( activateMemoryDiagnostic_ )
    NaluEnv::self().naluOutputP0() << "Nalu will activate detailed memory pulse" << std::endl;
  
  // allow for inconsistent restart (fields are missing)
  get_if_present(node, "support_inconsistent_multi_state_restart", supportInconsistentRestart_, supportInconsistentRestart_);

  // time step control
  const bool dtOptional = true;
  const YAML::Node y_time_step = expect_map(node,"time_step_control", dtOptional);
  if ( y_time_step ) {
    get_if_present(y_time_step, "target_courant", targetCourant_, targetCourant_);
    get_if_present(y_time_step, "time_step_change_factor", timeStepChangeFactor_, timeStepChangeFactor_);
  }

  if (node["balance_nodes_iterations"] || node["balance_nodes_target"] ) {
    throw std::runtime_error("balance nodes not supported");
  }

  //======================================
  // now other commands/actions
  //======================================

  // load output first so we can check for serializing i/o
  outputInfo_->load(node);
  if (root()->serializedIOGroupSize_ == 0)
  {
    // only set from input file if command-line didn't set it
    root()->setSerializedIOGroupSize(outputInfo_->serializedIOGroupSize_);
  }

  // solution options - loaded before create_mesh
  solutionOptions_->load(node);

  // once we know the mesh name, we can open the meta data, and set spatial dimension
  create_mesh();
  spatialDimension_ = meta_data().spatial_dimension();

  // post processing
  //postProcessingInfo_->load(node);

  // boundary, init, material and equation systems "load"
  if ( type_ == "multi_physics" ) {
    NaluEnv::self().naluOutputP0() << std::endl;
    NaluEnv::self().naluOutputP0() << "Boundary Condition Review: " << std::endl;
    NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
    boundaryConditions_.load(node);
    NaluEnv::self().naluOutputP0() << std::endl;
    NaluEnv::self().naluOutputP0() << "Initial Condition Review:  " << std::endl;
    NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
    initialConditions_.load(node);
    NaluEnv::self().naluOutputP0() << std::endl;
    NaluEnv::self().naluOutputP0() << "Material Prop Review:      " << std::endl;
    NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
    materialPropertys_.load(node);
    NaluEnv::self().naluOutputP0() << std::endl;
    NaluEnv::self().naluOutputP0() << "EqSys/options Review:      " << std::endl;
    NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
    equationSystems_.load(node);
  }

  // set number of nodes, check job run size
  check_job(true);
}

Simulation *Realm::root() { return parent()->root(); }
Simulation *Realm::root() const { return parent()->root(); }
Realms *Realm::parent() { return &realms_; }
Realms *Realm::parent() const { return &realms_; }

//--------------------------------------------------------------------------
//-------- setup_nodal_fields ----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_nodal_fields()
{
  // register global id and rank fields on all parts
  const stk::mesh::PartVector parts = meta_data().get_parts();
  for ( size_t ipart = 0; ipart < parts.size(); ++ipart ) {
    naluGlobalId_ = &(meta_data().declare_field<stk::mesh::EntityId>(stk::topology::NODE_RANK, "nalu_global_id"));
    stk::mesh::put_field_on_mesh(*naluGlobalId_, *parts[ipart], nullptr);
  }
  
  // loop over all material props targets and register nodal fields
  std::vector<std::string> targetNames = get_physics_target_names();
  equationSystems_.register_nodal_fields(targetNames);
}

//--------------------------------------------------------------------------
//-------- setup_edge_fields -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_edge_fields()
{
  // removed
}

//--------------------------------------------------------------------------
//-------- setup_element_fields --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_element_fields()
{
  // loop over all material props targets and register element fields
  std::vector<std::string> targetNames = get_physics_target_names();
  equationSystems_.register_element_fields(targetNames);
}

//--------------------------------------------------------------------------
//-------- setup_interior_algorithms ---------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_interior_algorithms()
{
  // loop over all material props targets and register interior algs
  std::vector<std::string> targetNames = get_physics_target_names();
  equationSystems_.register_interior_algorithm(targetNames);
}

//--------------------------------------------------------------------------
//-------- setup_post_processing_algorithms --------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_post_processing_algorithms()
{
  // removed
}

//--------------------------------------------------------------------------
//-------- setup_bc --------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_bc()
{
  // loop over all bcs and register
  for (size_t ibc = 0; ibc < boundaryConditions_.size(); ++ibc) {
    BoundaryCondition& bc = *boundaryConditions_[ibc];
    std::string name = physics_part_name(bc.targetName_);

    switch(bc.theBcType_) {
      case SYMMETRY_BC:
        equationSystems_.register_symmetry_bc(name, *reinterpret_cast<const SymmetryBoundaryConditionData *>(&bc));
        break;
      case PERIODIC_BC:
      {
        STK_ThrowAssert(reinterpret_cast<const PeriodicBoundaryConditionData *>(&bc) != nullptr);
        const auto& pbc = (*reinterpret_cast<const PeriodicBoundaryConditionData *>(&bc));

        std::string monarchName = physics_part_name(pbc.monarchSubject_.monarch_);
        std::string subjectName = physics_part_name(pbc.monarchSubject_.subject_);
        equationSystems_.register_periodic_bc(monarchName, subjectName, pbc);
        break;
      }
      default:
        throw std::runtime_error("unknown bc");
    }
  }
}

//--------------------------------------------------------------------------
//-------- enforce_bc_on_exposed_faces  ------------------------------------
//--------------------------------------------------------------------------
void
Realm::enforce_bc_on_exposed_faces()
{
  double start_time = NaluEnv::self().nalu_time();

  NaluEnv::self().naluOutputP0() << "Realm::skin_mesh(): Begin" << std::endl;

  // first, skin mesh and, therefore, populate
  stk::mesh::Selector activePart = meta_data().locally_owned_part() | meta_data().globally_shared_part();
  stk::mesh::PartVector partVec;
  partVec.push_back(exposedBoundaryPart_);
  stk::mesh::create_exposed_block_boundary_sides(*bulkData_, activePart, partVec);

  stk::mesh::Selector selectRule = stk::mesh::Selector(*exposedBoundaryPart_)
    & !stk::mesh::selectUnion(bcPartVec_);

  stk::mesh::BucketVector const& face_buckets = bulkData_->get_buckets(meta_data().side_rank(), selectRule);

  if (!face_buckets.empty()) {
    NaluEnv::self().naluOutputP0() << "Exposed surfaces found without a boundary condition applied" << std::endl;

    // proceed to show the problem faces
    for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
            ib != face_buckets.end() ; ++ib )
    {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        // extract the face
        stk::mesh::Entity face = b[k];
        
        // report the offending face id
        NaluEnv::self().naluOutput() << "Face Id: " << bulkData_->identifier(face) << " is not properly covered" << std::endl;
      
        // extract face nodes
        const stk::mesh::Entity* face_node_rels = bulkData_->begin_nodes(face); 
        const unsigned numberOfNodes = bulkData_->num_nodes(face);
        NaluEnv::self().naluOutput() << " Number of nodes connected to this face is: " << numberOfNodes << std::endl;
        for ( unsigned n = 0; n < numberOfNodes; ++n ) {
          stk::mesh::Entity node = face_node_rels[n];
          NaluEnv::self().naluOutput() << " attached node Id: " << bulkData_->identifier(node) << std::endl;
        }
      
        // extract the element relations to report to the user and the number of elements connected
        const stk::mesh::Entity* face_elem_rels = bulkData_->begin_elements(face);
        const unsigned numberOfElems = bulkData_->num_elements(face);
        NaluEnv::self().naluOutput() << " Number of elements connected to this face is: " << numberOfElems << std::endl;

        for ( unsigned faceElem = 0; faceElem < numberOfElems; ++faceElem ) {
          stk::mesh::Entity element = face_elem_rels[faceElem];
          NaluEnv::self().naluOutput() << " attached element Id: " << bulkData_->identifier(element) << std::endl;
        }
      }
    }
    throw std::runtime_error("Realm::Error: Please aply bc to problematic exposed surfaces ");
  }

  const double end_time = NaluEnv::self().nalu_time();

  // set mesh reading
  timerSkinMesh_ = (end_time - start_time);

  NaluEnv::self().naluOutputP0() << "Realm::skin_mesh(): End" << std::endl;
}

//--------------------------------------------------------------------------
//-------- setup_initial_conditions ----------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_initial_conditions()
{
  // loop over all ics and register
  for (size_t j_ic = 0; j_ic < initialConditions_.size(); ++j_ic) {
    InitialCondition& initCond = *initialConditions_[j_ic];

    const std::vector<std::string> targetNames = initCond.targetNames_;

    for (size_t itarget=0; itarget < targetNames.size(); ++itarget) {
      const std::string targetName = physics_part_name(targetNames[itarget]);

      // target need not be subsetted since nothing below will depend on topo
      stk::mesh::Part *targetPart = meta_data().get_part(targetName);

      switch(initCond.theIcType_) {

        case CONSTANT_UD:
        {
          const ConstantInitialConditionData& genIC = *reinterpret_cast<const ConstantInitialConditionData *>(&initCond);
          STK_ThrowAssert(genIC.data_.size() == genIC.fieldNames_.size());
          for (size_t ifield = 0; ifield < genIC.fieldNames_.size(); ++ifield) {

            std::vector<double>  genSpec = genIC.data_[ifield];
            stk::mesh::FieldBase *field = stk::mesh::get_field_by_name(genIC.fieldNames_[ifield], meta_data());

            if ( nullptr == field )
              throw std::runtime_error("Realm::setup_initial_conditions: field is null: " + genIC.fieldNames_[ifield]);

            stk::mesh::FieldBase *fieldWithState = ( field->number_of_states() > 1 )
              ? field->field_state(stk::mesh::StateNP1)
              : field->field_state(stk::mesh::StateNone);

            std::vector<double> userGen = genSpec;
            ConstantAuxFunction *theGenFunc = new ConstantAuxFunction(0, genSpec.size(), userGen);
            AuxFunctionAlgorithm *auxGen
              = new AuxFunctionAlgorithm( *this, targetPart,
                                          fieldWithState, theGenFunc, stk::topology::NODE_RANK);
            initCondAlg_.push_back(auxGen);

          }
        }
        break;

        case FUNCTION_UD:
        {
          const UserFunctionInitialConditionData& fcnIC = *reinterpret_cast<const UserFunctionInitialConditionData *>(&initCond);
          equationSystems_.register_initial_condition_fcn(targetPart, fcnIC);
        }
        break;

        case USER_SUB_UD:
          throw std::runtime_error("Realm::setup_initial_conditions: USER_SUB not supported: ");

        case UserDataType_END:
          break;

        default:
          NaluEnv::self().naluOutputP0() << "Realm::setup_initial_conditions: unknown type: " << initCond.theIcType_ << std::endl;
          throw std::runtime_error("Realm::setup_initial_conditions: unknown type:");
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- setup_property --------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::setup_property()
{
  // high level defaults, which can be modified, for [pressure], temperature and R
  double tRef = 298.15;
  double universalR = 8314.4621;
         
  // loop overall material property blocks
  for ( size_t i = 0; i < materialPropertys_.materialPropertyVector_.size(); ++i ) {
    
    MaterialProperty *matPropBlock = materialPropertys_.materialPropertyVector_[i];
    
    // loop over all target names
    for (size_t itarget=0; itarget < matPropBlock->targetNames_.size(); ++itarget) {
      
      // target need not be subsetted since nothing below will depend on topo
      const std::string physicsPartName = physics_part_name(matPropBlock->targetNames_[itarget]);
      stk::mesh::Part *targetPart = meta_data().get_part(physicsPartName);
      
      // loop over propertyMap
      std::map<PropertyIdentifier, ScalarFieldType *>::iterator ii;
      for ( ii=propertyMap_.begin(); ii!=propertyMap_.end(); ++ii ) {
        
        // extract property id and field pointer
        PropertyIdentifier thePropId = (*ii).first;
        ScalarFieldType *thePropField = (*ii).second;
        
        // find the material property data object
        MaterialPropertyData *matData = NULL;
        std::map<PropertyIdentifier, MaterialPropertyData*>::iterator itf =
          matPropBlock->propertyDataMap_.find(thePropId);
        if ( itf == matPropBlock->propertyDataMap_.end() ) {
          // will need to throw
          NaluEnv::self().naluOutputP0() << "issue with property: " << PropertyIdentifierNames[thePropId] << std::endl;
          throw std::runtime_error("Please add property specification ");
        }
        else {
          matData = (*itf).second;
        }
        
        switch( matData->type_) {

        case CONSTANT_MAT:
        { 
          if ( thePropId == GAMMA_ID ) {
            // error checks
            if ( !gasDynamics_ )
              throw std::runtime_error("gas dynamics is the only physics that requires gamma");

            if ( matData->cpConstMap_.size() > 0.0 )
              throw std::runtime_error("gas dynamics does not support species-based Cp");

            if ( !solutionOptions_->accousticallyCompressible_ )
              throw std::runtime_error("gas dynamics requires an accoustically compressible key word");

            // save off constant value
            const double gamma = matData->constValue_;

            // extract reference values
            matPropBlock->extract_universal_constant("reference_temperature", tRef, true);
            matPropBlock->extract_universal_constant("universal_gas_constant", universalR, true);
            
            // extract fields
            ScalarFieldType *cpField 
              = meta_data().get_field<double>(stk::topology::NODE_RANK, "specific_heat");
            ScalarFieldType *cvField 
              = meta_data().get_field<double>(stk::topology::NODE_RANK, "specific_heat_v");

            // Compute Cp and Cv on the fly based on gamma and R/mixtureMW
            std::vector<std::pair<double, double> > mwMassFracVec;
            std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
            for ( itrp = matPropBlock->referencePropertyDataMap_.begin();
                  itrp!= matPropBlock->referencePropertyDataMap_.end(); ++itrp) {
              ReferencePropertyData *propData = (*itrp).second;
              std::pair<double,double> thePair;
              thePair = std::make_pair(propData->mw_,propData->massFraction_);
              mwMassFracVec.push_back(thePair);
            }
                
            double mixtureMW = 0.0;
            for (std::size_t k = 0; k < mwMassFracVec.size(); ++k ){
              mixtureMW += mwMassFracVec[k].second/mwMassFracVec[k].first;
            }
            mixtureMW = 1.0/mixtureMW;

            // compute Cp and Cv with constant gamma
            const double Cp = gamma*universalR/mixtureMW/(gamma - 1.0);
            const double Cv = universalR/mixtureMW/(gamma - 1.0);
            
            // create the nodal population for gamma
            std::vector<double> userGamma(1);
            userGamma[0] = gamma;
            ConstantAuxFunction *auxFuncGamma
              = new ConstantAuxFunction(0, 1, userGamma);
            AuxFunctionAlgorithm *auxAlgGamma
              = new AuxFunctionAlgorithm( *this, targetPart,
                                          thePropField, auxFuncGamma, stk::topology::NODE_RANK);
            propertyAlg_.push_back(auxAlgGamma);
                        
            // create the nodal population for Cp
            std::vector<double> userCp(1);            
            userCp[0] = Cp;
            ConstantAuxFunction *auxFuncCp
              = new ConstantAuxFunction(0, 1, userCp);
            AuxFunctionAlgorithm *auxAlgCp
              = new AuxFunctionAlgorithm( *this, targetPart,
                                          cpField, auxFuncCp, stk::topology::NODE_RANK);
            propertyAlg_.push_back(auxAlgCp);

            // create the nodal population for Cv
            std::vector<double> userCv(1);            
            userCv[0] = Cv;
            ConstantAuxFunction *auxFuncCv
              = new ConstantAuxFunction(0, 1, userCv);
            AuxFunctionAlgorithm *auxAlgCv
              = new AuxFunctionAlgorithm( *this, targetPart,
                                          cvField, auxFuncCv, stk::topology::NODE_RANK);
            propertyAlg_.push_back(auxAlgCv);

            // deal with enthalpy property evaluator
            PropertyEvaluator *theEnthPropEval 
              = new EnthalpyConstSpecHeatPropertyEvaluator(Cp, tRef);
            
            // push to prop eval
            matPropBlock->propertyEvalMap_[ENTHALPY_ID]  = theEnthPropEval;
          }
          else {            
            // all other properties are currently constant
            int theBegin = 0;
            int theEnd = 1;
            
            // create everything
            std::vector<double> userConstData(1);
            userConstData[0] = matData->constValue_;
            ConstantAuxFunction *theAuxFunc
              = new ConstantAuxFunction(theBegin, theEnd, userConstData);
            AuxFunctionAlgorithm *auxAlg
              = new AuxFunctionAlgorithm( *this, targetPart,
					  thePropField, theAuxFunc, stk::topology::NODE_RANK);
            propertyAlg_.push_back(auxAlg);
          }
        }
        break;

        case IDEAL_GAS_MAT:
          {
          if ( DENSITY_ID == thePropId ) {
        
            // require R
            matPropBlock->extract_universal_constant("universal_gas_constant", universalR, true);
            
            // placeholder for the property evaluator and alg
            PropertyEvaluator *rhoPropEval = NULL;
            
            if ( solutionOptions_->accousticallyCompressible_ ) {
              // load mw and reference species
              std::vector<std::pair<double, double> > mwMassFracVec;
              std::map<std::string, ReferencePropertyData*>::const_iterator itrp;
              for ( itrp = matPropBlock->referencePropertyDataMap_.begin();
                    itrp!= matPropBlock->referencePropertyDataMap_.end(); ++itrp) {
                ReferencePropertyData *propData = (*itrp).second;
                std::pair<double,double> thePair;
                thePair = std::make_pair(propData->mw_,propData->massFraction_);
                mwMassFracVec.push_back(thePair);
              }
              // rho = f(p,T,mwRef_)
              rhoPropEval = new IdealGasPTYkrefPropertyEvaluator(universalR, mwMassFracVec, meta_data());
            }
            else {
              throw std::runtime_error("Realm::setup_property: ideal_gas only supported for compressible scheme:");
            }
              
            // push back property evaluator to map
            matPropBlock->propertyEvalMap_[thePropId] = rhoPropEval;
            
          }
          else {
            throw std::runtime_error("Realm::setup_property: ideal_gas_t only supported for density:");
          }
        }
        break;
        
        case GEOMETRIC_MAT:
        {
          throw std::runtime_error("Realm::setup_property: no geometric_mat support:");
        }
        break;
      
        case MaterialPropertyType_END:
          break;
          
        default:
          throw std::runtime_error("Realm::setup_property: unknown type:");
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- initialize_global_variables -------------------------------------
//--------------------------------------------------------------------------
void
Realm::initialize_global_variables()
{
  // other variables created on the fly during Eqs registration
  const bool needInOutput = false;
  const bool needInRestart = true;
  globalParameters_.set_param("timeStepNm1", 1.0, needInOutput, needInRestart);
  globalParameters_.set_param("timeStepCount", 1, needInOutput, needInRestart);
}

//--------------------------------------------------------------------------
//-------- augment_property_map --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_property_map(
  PropertyIdentifier propID,
  ScalarFieldType *theField)
{
  propertyMap_[propID] = theField;
}

//--------------------------------------------------------------------------
//-------- makeSureNodesHaveValidTopology ----------------------------------
//--------------------------------------------------------------------------
void 
Realm::makeSureNodesHaveValidTopology()
{
  //To make sure nodes have valid topology, we have to make sure they are in a part that has NODE topology.
  //So first, let's obtain the node topology part:
  stk::mesh::Part& nodePart = bulkData_->mesh_meta_data().get_topology_root_part(stk::topology::NODE);
  stk::mesh::Selector nodesNotInNodePart = (!nodePart) & bulkData_->mesh_meta_data().locally_owned_part();

  //get all the nodes that are *NOT* in nodePart
  std::vector<stk::mesh::Entity> nodes_vector;
  stk::mesh::get_selected_entities(nodesNotInNodePart, bulkData_->buckets(stk::topology::NODE_RANK), nodes_vector);
  // now we require all nodes are in proper node part
  if (nodes_vector.size())
    std::cout << "nodes_vector= " << nodes_vector.size() << std::endl;
  STK_ThrowRequire(0 == nodes_vector.size());
}

//--------------------------------------------------------------------------
//-------- pre_timestep_work -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::pre_timestep_work()
{
  // check for mesh motion
  if ( solutionOptions_->meshMotion_ ) {
    throw std::runtime_error("Mesh motion (sliding/six-DOF) deactivated");
  }

  // deal with non-topology changes, however, moving mesh
  if ( has_mesh_deformation() ) {
    // extract target parts for this physics
    if ( solutionOptions_->externalMeshDeformation_ ) {
      std::vector<std::string> targetNames = get_physics_target_names();
      for ( size_t itarget = 0; itarget < targetNames.size(); ++itarget ) {
        stk::mesh::Part *targetPart = meta_data().get_part(targetNames[itarget]);
        set_current_coordinates(targetPart);
      }
    }
    compute_geometry();
  }

  // ask the equation system to do some work
  equationSystems_.pre_timestep_work();
}

//--------------------------------------------------------------------------
//-------- evaluate_properties ---------------------------------------------
//--------------------------------------------------------------------------
void
Realm::evaluate_properties()
{
  double start_time = NaluEnv::self().nalu_time();
  for ( size_t k = 0; k < propertyAlg_.size(); ++k ) {
    propertyAlg_[k]->execute();
  }
  equationSystems_.evaluate_properties();
  double end_time = NaluEnv::self().nalu_time();
  timerPropertyEval_ += (end_time - start_time);
}

//--------------------------------------------------------------------------
//-------- advance_time_step -----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::advance_time_step()
{
  // leave if we do not need to solve
  if ( !active_time_step() )
    return;

  NaluEnv::self().naluOutputP0() << name_ << "::advance_time_step() " << std::endl;

  NaluEnv::self().naluOutputP0() << "NLI"
                  << std::setw(8) << std::right << "Name"
                  << std::setw(22) << std::right << "Linear Iter"
                  << std::setw(16) << std::right << "Linear Res"
                  << std::setw(16) << std::right << "NLinear Res"
                  << std::setw(14) << std::right << "Scaled NLR" << std::endl;

  NaluEnv::self().naluOutputP0() << "---"
                  << std::setw(8) << std::right << "----"
                  << std::setw(22) << std::right << "-----------"
                  << std::setw(16) << std::right << "----------"
                  << std::setw(16) << std::right << "-----------"
                  << std::setw(14) << std::right << "----------" << std::endl;

  // evaluate new geometry based on latest mesh motion geometry state (provided that external is active)
  if ( solutionOptions_->externalMeshDeformation_ )
    compute_geometry();

  // evaluate properties based on latest state including boundary and and possible xfer
  evaluate_properties();

  // compute velocity relative to mesh
  compute_vrtm();

  const int numNonLinearIterations = equationSystems_.maxIterations_;
  for ( int i = 0; i < numNonLinearIterations; ++i ) {
    currentNonlinearIteration_ = i+1;
    NaluEnv::self().naluOutputP0()
      << currentNonlinearIteration_
      << "/" << numNonLinearIterations
      << std::setw(29) << std::right << "Equation System Iteration" << std::endl;

    isFinalOuterIter_ = ((i+1) == numNonLinearIterations);

    const bool isConverged = equationSystems_.solve_and_update();

    // evaluate properties based on latest np1 solution
    evaluate_properties();

    if ( isConverged ) {
      NaluEnv::self().naluOutputP0() << "norm convergence criteria met for all equation systems: " << std::endl;
      NaluEnv::self().naluOutputP0() << "max scaled norm is: " << equationSystems_.provide_system_norm() << std::endl;
      break;
    }
  }

}

//--------------------------------------------------------------------------
//-------- active_time_step ------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::active_time_step()
{
  const int timeStepCount = get_time_step_count();
  const bool activeTimeStep = (timeStepCount % solveFrequency_ ) == 0 ? true : false;
  return activeTimeStep;
}

//--------------------------------------------------------------------------
//-------- output_converged_results ----------------------------------------
//--------------------------------------------------------------------------
void
Realm::output_converged_results()
{
  provide_output();
  provide_restart_output();
}

//--------------------------------------------------------------------------
//-------- compute_adaptive_time_step --------------------------------------
//--------------------------------------------------------------------------
double
Realm::compute_adaptive_time_step()
{
  // extract current time
  const double dtN = get_time_step();

  // ratio of how off we are
  const double factorOff = targetCourant_/maxCourant_;

  // scaling for dt and candidate
  const double dtScaling = ( targetCourant_ < maxCourant_ )
    ? std::max(factorOff, 1.0/timeStepChangeFactor_)
    : std::min(factorOff, timeStepChangeFactor_);
  const double candidateDt = dtN*dtScaling;

  return candidateDt;
}

//--------------------------------------------------------------------------
//-------- commit ----------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::commit()
{
  //====================================================
  // Commit the meta data
  //====================================================
  meta_data().commit();
}

//--------------------------------------------------------------------------
//-------- create_mesh() ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::create_mesh()
{
  double start_time = NaluEnv::self().nalu_time();

  NaluEnv::self().naluOutputP0() << "Realm::create_mesh(): Begin" << std::endl;
  stk::ParallelMachine pm = NaluEnv::self().parallel_comm();

  // news for mesh constructs  
  stk::mesh::MeshBuilder builder(pm);
  builder.set_aura_option(activateAura_ ? stk::mesh::BulkData::AUTO_AURA : stk::mesh::BulkData::NO_AUTO_AURA);
  bulkData_ = builder.create();
  bulkData_->mesh_meta_data().use_simple_fields();
 
  ioBroker_ = new stk::io::StkMeshIoBroker(pm);
  ioBroker_->set_bulk_data(bulkData_);
  
  // allow for automatic decomposition
  if (autoDecompType_ != "None") 
    ioBroker_->property_add(Ioss::Property("DECOMPOSITION_METHOD", autoDecompType_));
  
  // Initialize meta data (from exodus file); can possibly be a restart file..
  inputMeshIdx_ = ioBroker_->add_mesh_database( 
   inputDBName_, restarted_simulation() ? stk::io::READ_RESTART : stk::io::READ_MESH );
  ioBroker_->create_input_mesh();

  // declare an exposed part for later bc coverage check
  if ( checkForMissingBcs_ ) {
    exposedBoundaryPart_ = &meta_data().declare_part("exposed_boundary_part",meta_data().side_rank());
  }

  // set mesh creation
  const double end_time = NaluEnv::self().nalu_time();
  timerCreateMesh_ = (end_time - start_time);

  NaluEnv::self().naluOutputP0() << "Realm::create_mesh() End" << std::endl;
}

//--------------------------------------------------------------------------
//-------- create_output_mesh() --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::create_output_mesh()
{
  // exodus output file creation
  if (outputInfo_->hasOutputBlock_ ) {

    double start_time = NaluEnv::self().nalu_time();
    NaluEnv::self().naluOutputP0() << "Realm::create_output_mesh(): Begin" << std::endl;

    if (outputInfo_->outputFreq_ == 0)
      return;

    std::string oname =  outputInfo_->outputDBName_ ;

#ifdef NALU_USES_CATALYST    
    if(!outputInfo_->catalystFileName_.empty()||
       !outputInfo_->paraviewScriptName_.empty()) {

      outputInfo_->outputPropertyManager_->add(Ioss::Property("CATALYST_BLOCK_PARSE_INPUT_DECK_NAME", 
                                                              NaluEnv::self().get_base_name()));
                                               
      if(!outputInfo_->catalystFileName_.empty()) {
        outputInfo_->outputPropertyManager_->add(Ioss::Property("PHACTORI_INPUT_SYNTAX_SCRIPT",
                                                                outputInfo_->catalystFileName_));
      }
      else if(!outputInfo_->paraviewScriptName_.empty()) {
        outputInfo_->outputPropertyManager_->add(Ioss::Property("CATALYST_SCRIPT",
            outputInfo_->paraviewScriptName_.c_str()));
      }
      
      outputInfo_->outputPropertyManager_->add(Ioss::Property("CATALYST_CREATE_SIDE_SETS", 1));
      
      resultsFileIndex_ = ioBroker_->create_output_mesh( oname, stk::io::WRITE_RESULTS,
          *outputInfo_->outputPropertyManager_, "catalyst_exodus" );
    }
    else {
      resultsFileIndex_ = ioBroker_->create_output_mesh( oname, stk::io::WRITE_RESULTS, *outputInfo_->outputPropertyManager_);
    }
#else
    resultsFileIndex_ = ioBroker_->create_output_mesh( oname, stk::io::WRITE_RESULTS, *outputInfo_->outputPropertyManager_);
#endif
    
    // Tell stk_io how to output element block nodal fields:
    // if 'true' passed to function, then output them as nodeset fields;
    // if 'false', then output as nodal fields (on all nodes of the mesh, zero-filled)
    // The option is provided since some post-processing/visualization codes do not
    // correctly handle nodeset fields.
    ioBroker_->use_nodeset_for_part_nodes_fields(resultsFileIndex_, outputInfo_->outputNodeSet_);

    // FIXME: add_field can take user-defined output name, not just varName
    for ( std::set<std::string>::iterator itorSet = outputInfo_->outputFieldNameSet_.begin();
        itorSet != outputInfo_->outputFieldNameSet_.end(); ++itorSet ) {
      std::string varName = *itorSet;
      stk::mesh::FieldBase *theField = stk::mesh::get_field_by_name(varName, meta_data());
      if ( NULL == theField ) {
        NaluEnv::self().naluOutputP0() << " Sorry, no field by the name " << varName << std::endl;
      }
      else {
        // 'varName' is the name that will be written to the database
        // For now, just using the name of the stk field
        ioBroker_->add_field(resultsFileIndex_, *theField, varName);
      }
    }

    // set mesh creation
    const double end_time = NaluEnv::self().nalu_time();
    timerCreateMesh_ = (end_time - start_time);

    NaluEnv::self().naluOutputP0() << "Realm::create_output_mesh() End" << std::endl;
  }
}

//--------------------------------------------------------------------------
//-------- create_restart_mesh() -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::create_restart_mesh()
{
  // exodus restart file creation
  if (outputInfo_->hasRestartBlock_ ) {

    if (outputInfo_->restartFreq_ == 0)
      return;
    
    restartFileIndex_ = ioBroker_->create_output_mesh(outputInfo_->restartDBName_, stk::io::WRITE_RESTART, *outputInfo_->restartPropertyManager_);
    
    // loop over restart variable field names supplied by Eqs
    for ( std::set<std::string>::iterator itorSet = outputInfo_->restartFieldNameSet_.begin();
        itorSet != outputInfo_->restartFieldNameSet_.end(); ++itorSet ) {
      std::string varName = *itorSet;
      stk::mesh::FieldBase *theField = stk::mesh::get_field_by_name(varName,meta_data());
      if ( NULL == theField ) {
        NaluEnv::self().naluOutputP0() << " Sorry, no field by the name " << varName << std::endl;
      }
      else {
        // add the field for a restart output
        ioBroker_->add_field(restartFileIndex_, *theField, varName);
        // if this is a restarted simulation, we will need input
        if ( restarted_simulation() )
          ioBroker_->add_input_field(stk::io::MeshField(*theField, varName));
      }
    }

    // now global params
    stk::util::ParameterMapType::const_iterator i = globalParameters_.begin();
    stk::util::ParameterMapType::const_iterator iend = globalParameters_.end();
    for (; i != iend; ++i) {
      std::string parameterName = (*i).first;
      stk::util::Parameter parameter = (*i).second;
      if(parameter.toRestartFile) {
        ioBroker_->add_global(restartFileIndex_, parameterName, parameter);
      }
    }

    // set max size for restart data base
    ioBroker_->get_output_ioss_region(restartFileIndex_)->get_database()->set_cycle_count(outputInfo_->restartMaxDataBaseStepSize_);
  }

}

//--------------------------------------------------------------------------
//-------- input_variables_from_mesh() -------------------------------------
//--------------------------------------------------------------------------
void
Realm::input_variables_from_mesh()
{
  // no variables from an input mesh if this is a restart
  if ( !restarted_simulation()) {

    // check whether to snap or interpolate data; all fields treated the same
    const stk::io::MeshField::TimeMatchOption fieldInterpOption = solutionOptions_->inputVariablesInterpolateInTime_
        ? stk::io::MeshField::LINEAR_INTERPOLATION
        : stk::io::MeshField::CLOSEST;

    // check for periodic cycling of data based on start time and periodic time; scale time set to unity
    if ( solutionOptions_->inputVariablesPeriodicTime_ > 0.0 ) {
      ioBroker_->get_mesh_database(inputMeshIdx_)
        .set_periodic_time(solutionOptions_->inputVariablesPeriodicTime_, 
                           solutionOptions_->inputVariablesRestorationTime_, 
                           stk::io::InputFile::CYCLIC)
        .set_scale_time(1.0);
    }
    
    std::map<std::string, std::string>::const_iterator iter;
    for ( iter = solutionOptions_->inputVarFromFileMap_.begin();
          iter != solutionOptions_->inputVarFromFileMap_.end(); ++iter) {

      std::string varName = iter->first;
      std::string userName = iter->second;

      stk::mesh::FieldBase *theField = stk::mesh::get_field_by_name(varName,meta_data());
      if ( NULL == theField ) {
        NaluEnv::self().naluOutputP0() << " Sorry, no field by the name " << varName << std::endl;
      }
      else {
        ioBroker_->add_input_field(stk::io::MeshField(*theField, userName, fieldInterpOption));
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- augment_output_variable_list() ----------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_output_variable_list(
    const std::string fieldName)
{
  outputInfo_->outputFieldNameSet_.insert(fieldName);
}

//--------------------------------------------------------------------------
//-------- augment_restart_variable_list -----------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_restart_variable_list(
  std::string restartFieldName)
{
  outputInfo_->restartFieldNameSet_.insert(restartFieldName);
}

//--------------------------------------------------------------------------
//-------- provide_entity_count() ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_entity_count() {

  std::vector<size_t> counts;
  std::vector<size_t> minCounts;
  std::vector<size_t> maxCounts;
  stk::mesh::comm_mesh_counts( *bulkData_ , counts, minCounts, maxCounts);

  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
  NaluEnv::self().naluOutputP0() << "Realm::provide_entity_count:   " << std::endl
		  << "nodes,    " << counts[0] << " min/max: " << minCounts[0] << "/" << maxCounts[0] << std::endl
		  << "edges,    " << counts[1] << " min/max: " << minCounts[1] << "/" << maxCounts[1] << std::endl
		  << "faces,    " << counts[2] << " min/max: " << minCounts[2] << "/" << maxCounts[2] << std::endl
		  << "elements, " << counts[3] << " min/max: " << minCounts[3] << "/" << maxCounts[3] << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;
}

//--------------------------------------------------------------------------
//-------- get_coordinates_name --------------------------------------------
//--------------------------------------------------------------------------
std::string
Realm::get_coordinates_name()
{
  return solutionOptions_->get_coordinates_name();
}

//--------------------------------------------------------------------------
//-------- has_mesh_motion -------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::has_mesh_motion() const
{
  return solutionOptions_->has_mesh_motion();
}

//--------------------------------------------------------------------------
//-------- has_mesh_deformation --------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::has_mesh_deformation() const
{
  return solutionOptions_->has_mesh_deformation();
}

//--------------------------------------------------------------------------
//-------- does_mesh_move --------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::does_mesh_move() const
{
  return solutionOptions_->does_mesh_move();
}

//--------------------------------------------------------------------------
//-------- set_current_coordinates -----------------------------------------
//--------------------------------------------------------------------------
void
Realm::set_current_coordinates(
  stk::mesh::Part *targetPart)
{
  const int nDim = meta_data().spatial_dimension();

  VectorFieldType *modelCoords = meta_data().get_field<double>(stk::topology::NODE_RANK, "coordinates");
  VectorFieldType *currentCoords = meta_data().get_field<double>(stk::topology::NODE_RANK, "current_coordinates");
  VectorFieldType *displacement = meta_data().get_field<double>(stk::topology::NODE_RANK, "mesh_displacement");

  stk::mesh::Selector s_all_nodes = stk::mesh::Selector(*targetPart);

  stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * mCoords = stk::mesh::field_data(*modelCoords, b);
    double * cCoords = stk::mesh::field_data(*currentCoords, b);
    const double * dx = stk::mesh::field_data(*displacement, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        cCoords[offSet+j] = mCoords[offSet+j] + dx[offSet+j];
      }
    }
  }
}


//--------------------------------------------------------------------------
//-------- compute_geometry ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::compute_geometry()
{
  if ( !usesCVFEM_ )
    return;
  
  // interior and boundary
  computeGeometryAlgDriver_->execute();

  // find total volume if the mesh moves at all
  if ( does_mesh_move() ) {
    double totalVolume = 0.0;
    double maxVolume = -1.0e16;
    double minVolume = 1.0e16;

    ScalarFieldType *dualVolume = meta_data().get_field<double>(stk::topology::NODE_RANK, "dual_nodal_volume");

    stk::mesh::Selector s_local_nodes
      = meta_data().locally_owned_part() &stk::mesh::selectField(*dualVolume);

    stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_local_nodes );
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
	  ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      const double * dv = stk::mesh::field_data(*dualVolume, b);
      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const double theVol = dv[k];
        totalVolume += theVol;
        maxVolume = std::max(theVol, maxVolume);
        minVolume = std::min(theVol, minVolume);
      }
    }

    // get min, max and sum over processes
    double g_totalVolume = 0.0, g_minVolume = 0.0, g_maxVolume = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &minVolume, &g_minVolume, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &maxVolume, &g_maxVolume, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &totalVolume, &g_totalVolume, 1);

    NaluEnv::self().naluOutputP0() << " Volume  " << g_totalVolume
		    << " min: " << g_minVolume
		    << " max: " << g_maxVolume << std::endl;
    // set Realm data
    minDualVolume_ = g_minVolume;
  }
}

//--------------------------------------------------------------------------
//-------- compute_vrtm ----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::compute_vrtm()
{
  // compute velocity relative to mesh; must be tied to velocity update...
  if ( hasFluids_ && (solutionOptions_->meshMotion_ || solutionOptions_->externalMeshDeformation_) ) {
    const int nDim = meta_data().spatial_dimension();

    VectorFieldType *velocity = meta_data().get_field<double>(stk::topology::NODE_RANK, "velocity");
    VectorFieldType *meshVelocity = meta_data().get_field<double>(stk::topology::NODE_RANK, "mesh_velocity");
    VectorFieldType *velocityRTM = meta_data().get_field<double>(stk::topology::NODE_RANK, "velocity_rtm");

    stk::mesh::Selector s_all_nodes
       = (meta_data().locally_owned_part() | meta_data().globally_shared_part());

    stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
    for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
	  ib != node_buckets.end() ; ++ib ) {
      stk::mesh::Bucket & b = **ib ;
      const stk::mesh::Bucket::size_type length   = b.size();
      const double * uNp1 = stk::mesh::field_data(*velocity, b);
      const double * vNp1 = stk::mesh::field_data(*meshVelocity, b);
      double * vrtm = stk::mesh::field_data(*velocityRTM, b);

      for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
        const int offSet = k*nDim;
        for ( int j=0; j < nDim; ++j ) {
          vrtm[offSet+j] = uNp1[offSet+j] - vNp1[offSet+j];
        }
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- init_current_coordinates ----------------------------------------
//--------------------------------------------------------------------------
void
Realm::init_current_coordinates()
{

  const int nDim = meta_data().spatial_dimension();

  VectorFieldType *modelCoords = meta_data().get_field<double>(stk::topology::NODE_RANK, "coordinates");
  VectorFieldType *currentCoords = meta_data().get_field<double>(stk::topology::NODE_RANK, "current_coordinates");
  VectorFieldType *displacement = meta_data().get_field<double>(stk::topology::NODE_RANK, "mesh_displacement");

  stk::mesh::Selector s_all_nodes
    = (meta_data().locally_owned_part() | meta_data().globally_shared_part());

  stk::mesh::BucketVector const& node_buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_all_nodes );
  for ( stk::mesh::BucketVector::const_iterator ib = node_buckets.begin() ;
        ib != node_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    const double * mCoords = stk::mesh::field_data(*modelCoords, b);
    double * cCoords = stk::mesh::field_data(*currentCoords, b);
    double * dx = stk::mesh::field_data(*displacement, b);
    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {
      const int offSet = k*nDim;
      for ( int j = 0; j < nDim; ++j ) {
        dx[offSet+j] = 0.0; //RESTART...
        cCoords[offSet+j] = mCoords[offSet+j];
      }
    }
  }
}

//--------------------------------------------------------------------------
//-------- compute_l2_scaling ----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::compute_l2_scaling()
{
  // loop over all material propertys  and save off part vector
  stk::mesh::PartVector partVec;
  const std::vector<std::string> targetNames = get_physics_target_names();
  for (size_t itarget=0; itarget < targetNames.size(); ++itarget) {
    // target need not be subsetted since nothing below will depend on topo
    stk::mesh::Part *targetPart = meta_data().get_part(targetNames[itarget]);
    partVec.push_back(targetPart);
  }

  size_t totalNodes = 0;

  // selector for all locally owned nodes
  stk::mesh::Selector s_locally_owned_union =
    meta_data().locally_owned_part()
    &stk::mesh::selectUnion(partVec);

  stk::mesh::BucketVector const& node_bucket = bulkData_->get_buckets( stk::topology::NODE_RANK, s_locally_owned_union );

  for ( stk::mesh::BucketVector::const_iterator ib = node_bucket.begin() ;
        ib != node_bucket.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();
    totalNodes += length;
  }

  // Parallel assembly of total nodes
  size_t g_totalNodes = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &totalNodes, &g_totalNodes, 1);

  l2Scaling_ = 1.0/std::sqrt(g_totalNodes);

}


//--------------------------------------------------------------------------
//-------- register_nodal_fields -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_nodal_fields(
  stk::mesh::Part *part)
{
  // register high level common fields
  const int nDim = meta_data().spatial_dimension();
  
  // common fields for everyone cvfem-based
  if ( usesCVFEM_ ) {
    ScalarFieldType *dualNodalVolume = &(meta_data().declare_field<double>(stk::topology::NODE_RANK, "dual_nodal_volume"));
    stk::mesh::put_field_on_mesh(*dualNodalVolume, *part, nullptr);
  }

  // mesh motion/deformation is high level
  if ( solutionOptions_->meshMotion_ || solutionOptions_->externalMeshDeformation_) {
    VectorFieldType *displacement = &(meta_data().declare_field<double>(stk::topology::NODE_RANK, "mesh_displacement"));
    stk::mesh::put_field_on_mesh(*displacement, *part, nDim, nullptr);
    stk::io::set_field_output_type(*displacement, stk::io::FieldOutputType::VECTOR_3D);
    VectorFieldType *currentCoords = &(meta_data().declare_field<double>(stk::topology::NODE_RANK, "current_coordinates"));
    stk::mesh::put_field_on_mesh(*currentCoords, *part, nDim, nullptr);
    stk::io::set_field_output_type(*currentCoords, stk::io::FieldOutputType::VECTOR_3D);
    VectorFieldType *meshVelocity = &(meta_data().declare_field<double>(stk::topology::NODE_RANK, "mesh_velocity"));
    stk::mesh::put_field_on_mesh(*meshVelocity, *part, nDim, nullptr);
    stk::io::set_field_output_type(*meshVelocity, stk::io::FieldOutputType::VECTOR_3D);
    VectorFieldType *velocityRTM = &(meta_data().declare_field<double>(stk::topology::NODE_RANK, "velocity_rtm"));
    stk::mesh::put_field_on_mesh(*velocityRTM, *part, nDim, nullptr);
    stk::io::set_field_output_type(*velocityRTM, stk::io::FieldOutputType::VECTOR_3D);
    // only internal mesh motion requires rotation rate
    if ( solutionOptions_->meshMotion_ ) {
      ScalarFieldType *omega = &(meta_data().declare_field<double>(stk::topology::NODE_RANK, "omega"));
      stk::mesh::put_field_on_mesh(*omega, *part, nullptr);
    }
    // only external mesh deformation requires dvi/dxj (for GCL)
    if ( solutionOptions_->externalMeshDeformation_) {
      ScalarFieldType *divV = &(meta_data().declare_field<double>(stk::topology::NODE_RANK, "div_mesh_velocity"));
      stk::mesh::put_field_on_mesh(*divV, *part, nullptr);
    }
  }  
}

//--------------------------------------------------------------------------
//-------- register_interior_algorithm -------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_interior_algorithm(
  stk::mesh::Part *part)
{
  // Track parts that are registered to interior algorithms
  interiorPartVec_.push_back(part);

  //====================================================
  // Register interior algorithms
  //====================================================
  if ( !usesCVFEM_ )
    return;

  const AlgorithmType algType = INTERIOR;
  std::map<AlgorithmType, Algorithm *>::iterator it
    = computeGeometryAlgDriver_->algMap_.find(algType);
  if ( it == computeGeometryAlgDriver_->algMap_.end() ) {
    ComputeGeometryInteriorAlgorithm *theAlg
      = new ComputeGeometryInteriorAlgorithm(*this, part);
    computeGeometryAlgDriver_->algMap_[algType] = theAlg;
  }
  else {
    it->second->partVec_.push_back(part);
  }

}

//--------------------------------------------------------------------------
//-------- register_wall_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_wall_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);
}

//--------------------------------------------------------------------------
//-------- register_inflow_bc ----------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_inflow_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);
}

//--------------------------------------------------------------------------
//-------- register_open_bc ------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_open_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);
}

//--------------------------------------------------------------------------
//-------- register_symmetry_bc --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_symmetry_bc(
  stk::mesh::Part *part,
  const stk::topology &theTopo)
{
  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(part);
}

//--------------------------------------------------------------------------
//-------- register_periodic_bc --------------------------------------------
//--------------------------------------------------------------------------
void
Realm::register_periodic_bc(
  stk::mesh::Part *monarchMeshPart,
  stk::mesh::Part *subjectMeshPart,
  const double &searchTolerance,
  const std::string &searchMethodName)
{
  allPeriodicInteractingParts_.push_back(monarchMeshPart);
  allPeriodicInteractingParts_.push_back(subjectMeshPart);

  // push back the part for book keeping and, later, skin mesh
  bcPartVec_.push_back(monarchMeshPart);
  bcPartVec_.push_back(subjectMeshPart);

  if ( NULL == periodicManager_ ) {
    periodicManager_ = new PeriodicManager(*this);
    hasPeriodic_ = true;
  }

  // add the parts to the manager
  periodicManager_->add_periodic_pair(monarchMeshPart, subjectMeshPart, searchTolerance, searchMethodName);
}

//--------------------------------------------------------------------------
//-------- periodic_field_update -------------------------------------------
//--------------------------------------------------------------------------
void
Realm::periodic_field_update(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField,
  const bool bypassFieldCheck,
  const bool addSubjects,
  const bool setSubjects) const
{
  periodicManager_->apply_constraints(theField, sizeOfField, bypassFieldCheck, addSubjects, setSubjects);
}

//--------------------------------------------------------------------------
//-------- periodic_delta_solution_update ----------------------------------
//--------------------------------------------------------------------------
void
Realm::periodic_delta_solution_update(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField) const
{
  const bool bypassFieldCheck = true;
  const bool addSubjects = false;
  const bool setSubjects = true;
  periodicManager_->apply_constraints(theField, sizeOfField, bypassFieldCheck, addSubjects, setSubjects);
}

//--------------------------------------------------------------------------
//-------- periodic_max_field_update ---------------------------------------
//--------------------------------------------------------------------------
void
Realm::periodic_max_field_update(
  stk::mesh::FieldBase *theField,
  const unsigned &sizeOfField) const
{
  periodicManager_->apply_max_field(theField, sizeOfField);
}

//--------------------------------------------------------------------------
//-------- get_subject_part_vector -----------------------------------------
//--------------------------------------------------------------------------
const stk::mesh::PartVector &
Realm::get_subject_part_vector()
{
  if ( hasPeriodic_)
    return periodicManager_->get_subject_part_vector();
  else
    return emptyPartVector_;
}

//--------------------------------------------------------------------------
//-------- provide_output --------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_output()
{
  stk::diag::TimeBlock mesh_output_timeblock(Simulation::outputTimer());

  if ( outputInfo_->hasOutputBlock_ ) {

    if (outputInfo_->outputFreq_ == 0)
      return;

    const double start_time = NaluEnv::self().nalu_time();

    // process output via io
    const double currentTime = get_current_time();
    const int timeStepCount = get_time_step_count();
    const int modStep = timeStepCount - outputInfo_->outputStart_;

    // check for elapsed WALL time threshold
    bool forcedOutput = false;
    if ( outputInfo_->userWallTimeResults_.first) {
      const double elapsedWallTime = stk::wall_time() - wallTimeStart_;
      // find the max over all core
      double g_elapsedWallTime = 0.0;
      stk::all_reduce_max(NaluEnv::self().parallel_comm(), &elapsedWallTime, &g_elapsedWallTime, 1);
      // convert to hours
      g_elapsedWallTime /= 3600.0;
      // only force output the first time the timer is exceeded
      if ( g_elapsedWallTime > outputInfo_->userWallTimeResults_.second ) {
        forcedOutput = true;
        outputInfo_->userWallTimeResults_.first = false;
        NaluEnv::self().naluOutputP0()
            << "Realm::provide_output()::Forced Result output will be processed at current time: "
            << currentTime << std::endl;
        NaluEnv::self().naluOutputP0()
            <<  " Elapsed (max) WALL time: " << g_elapsedWallTime << " (hours)" << std::endl;
        // provide timer information
        dump_simulation_time();
      }
    }

    const bool isOutput 
      = (timeStepCount >=outputInfo_->outputStart_ && modStep % outputInfo_->outputFreq_ == 0) || forcedOutput;

    if ( isOutput ) {

      NaluEnv::self().naluOutputP0() << "Realm shall provide output files at : currentTime/timeStepCount: "
                                     << currentTime << "/" <<  timeStepCount << " (" << name_ << ")" << std::endl;      

      // not set up for globals
      ioBroker_->process_output_request(resultsFileIndex_, currentTime);
      
      equationSystems_.provide_output();
    }

    const double stop_time = NaluEnv::self().nalu_time();

    // increment time for output
    timerOutputFields_ += (stop_time - start_time);
  }
}

//--------------------------------------------------------------------------
//-------- provide_restart_output ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::provide_restart_output()
{
  stk::diag::TimeBlock mesh_output_timeblock(Simulation::outputTimer());

  if ( outputInfo_->hasRestartBlock_ ) {

    if (outputInfo_->restartFreq_ == 0)
      return;

    const double start_time = NaluEnv::self().nalu_time();

    // process restart via io
    const double currentTime = get_current_time();
    const int timeStepCount = get_time_step_count();
    const int modStep = timeStepCount - outputInfo_->restartStart_;

    // check for elapsed WALL time threshold
    bool forcedOutput = false;
    if ( outputInfo_->userWallTimeRestart_.first) {
      const double elapsedWallTime = stk::wall_time() - wallTimeStart_;
      // find the max over all core
      double g_elapsedWallTime = 0.0;
      stk::all_reduce_max(NaluEnv::self().parallel_comm(), &elapsedWallTime, &g_elapsedWallTime, 1);
      // convert to hours
      g_elapsedWallTime /= 3600.0;
      // only force output the first time the timer is exceeded
      if ( g_elapsedWallTime > outputInfo_->userWallTimeRestart_.second ) {
        forcedOutput = true;
        outputInfo_->userWallTimeRestart_.first = false;
        NaluEnv::self().naluOutputP0()
            << "Realm::provide_restart_output()::Forced Restart output will be processed at current time: "
            << currentTime << std::endl;
        NaluEnv::self().naluOutputP0()
            <<  " Elapsed (max) WALL time: " << g_elapsedWallTime << " (hours)" << std::endl;
      }
    }

    const bool isRestartOutputStep 
      = (timeStepCount >= outputInfo_->restartStart_ && modStep % outputInfo_->restartFreq_ == 0) || forcedOutput;
    
    if ( isRestartOutputStep ) {
      NaluEnv::self().naluOutputP0() << "Realm shall provide restart files at: currentTime/timeStepCount: "
                                     << currentTime << "/" <<  timeStepCount << " (" << name_ << ")" << std::endl;      
      // handle fields
      ioBroker_->begin_output_step(restartFileIndex_, currentTime);
      ioBroker_->write_defined_output_fields(restartFileIndex_);

      // push global variables for time step
      const double timeStepNm1 = timeIntegrator_->get_time_step();
      globalParameters_.set_value("timeStepNm1", timeStepNm1);
      globalParameters_.set_value("timeStepCount", timeStepCount);

      stk::util::ParameterMapType::const_iterator i = globalParameters_.begin();
      stk::util::ParameterMapType::const_iterator iend = globalParameters_.end();
      for (; i != iend; ++i)
      {
        std::string parameterName = (*i).first;
        stk::util::Parameter parameter = (*i).second;
        if ( parameter.toRestartFile ) {
          ioBroker_->write_global(restartFileIndex_, parameterName,  parameter);
        }
      }

      ioBroker_->end_output_step(restartFileIndex_);
    }

    const double stop_time = NaluEnv::self().nalu_time();

    // increment time for output
    timerOutputFields_ += (stop_time - start_time);
  }

}

//--------------------------------------------------------------------------
//-------- swap_states -----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::swap_states()
{
  bulkData_->update_field_data_states();
}

//--------------------------------------------------------------------------
//-------- predict_state ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::predict_state()
{
  equationSystems_.predict_state();
}

//--------------------------------------------------------------------------
//-------- populate_initial_condition --------------------------------------
//--------------------------------------------------------------------------
void
Realm::populate_initial_condition()
{
  for ( size_t k = 0; k < initCondAlg_.size(); ++k ) {
    initCondAlg_[k]->execute();
  }
}

//--------------------------------------------------------------------------
//-------- boundary_data_to_state_data -------------------------------------
//--------------------------------------------------------------------------
void
Realm::boundary_data_to_state_data()
{
  equationSystems_.boundary_data_to_state_data();
}

//--------------------------------------------------------------------------
//-------- populate_restart ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::populate_restart(
  double &timeStepNm1, int &timeStepCount)
{
  double foundRestartTime = get_current_time();
  if ( restarted_simulation() ) {
    // allow restart to skip missed required fields
    const double restartTime = outputInfo_->restartTime_;
    std::vector<stk::io::MeshField> missingFields;
    foundRestartTime = ioBroker_->read_defined_input_fields(restartTime, &missingFields);
    if ( missingFields.size() > 0 ){
      for ( size_t k = 0; k < missingFields.size(); ++k) {
        NaluEnv::self().naluOutputP0() << "WARNING: Restart value for Field "
                                       << missingFields[k].field()->name()
                                       << " is missing; may default to IC specification" << std::endl;
      }
      if ( !supportInconsistentRestart_ ) {
        NaluEnv::self().naluOutputP0() << "The user may desire to set the support_inconsistent_multi_state_restart Realm line command" << std::endl;
        NaluEnv::self().naluOutputP0() << "This is applicable for a BDF2 restart run from a previously run Backward Euler simulation" << std::endl;
      }
    }
    NaluEnv::self().naluOutputP0() << "Realm::populate_restart() candidate restart time: "
        << foundRestartTime << " for Realm: " << name() << std::endl;

    // extract time parameters; okay if they are missing; no need to let the user know
    const bool abortIfNotFound = false;
    ioBroker_->get_global("timeStepNm1", timeStepNm1, abortIfNotFound);
    ioBroker_->get_global("timeStepCount", timeStepCount, abortIfNotFound);

    // allow the user to reset the time; populate from the found data base, however, not reset to a user time
    if ( outputInfo_->restartResetTime_ ) {
      foundRestartTime = outputInfo_->restartResetNewTime_;
      NaluEnv::self().naluOutputP0() << "Realm::populate_restart() candidate restart time: "
                                     << foundRestartTime << " for Realm: " << name() << " (override via reset)" 
                                     << std::endl;
    }
  }
  return foundRestartTime;
}

//--------------------------------------------------------------------------
//-------- populate_variables_from_input -----------------------------------
//--------------------------------------------------------------------------
double
Realm::populate_variables_from_input(const double currentTime)
{
  // no reading fields from mesh if this is a restart
  double foundTime = currentTime;
  if ( !restarted_simulation() && solutionOptions_->inputVarFromFileMap_.size() > 0 ) {
    std::vector<stk::io::MeshField> missingFields;
    foundTime = ioBroker_->read_defined_input_fields(solutionOptions_->inputVariablesRestorationTime_, &missingFields);
    if ( missingFields.size() > 0 ) {
      for ( size_t k = 0; k < missingFields.size(); ++k) {
        NaluEnv::self().naluOutputP0() << "WARNING: Realm::populate_variables_from_input for field "
            << missingFields[k].field()->name()
            << " is missing; will default to IC specification" << std::endl;
      }
    }
    NaluEnv::self().naluOutputP0() << "Realm::populate_variables_form_input() candidate input time: "
        << foundTime << " for Realm: " << name() << std::endl;
  }
  return foundTime;
}

//--------------------------------------------------------------------------
//-------- populate_derived_quantities -------------------------------------
//--------------------------------------------------------------------------
void
Realm::populate_derived_quantities()
{
  equationSystems_.populate_derived_quantities();
}

//--------------------------------------------------------------------------
//-------- initial_work ----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::initial_work()
{
  if ( solutionOptions_->meshMotion_ )
    throw std::runtime_error("Mesh motion (sliding/six-DOF) deactivated");
  equationSystems_.initial_work();
}

//--------------------------------------------------------------------------
//-------- set_global_id ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::set_global_id()
{
  const stk::mesh::Selector s_universal = meta_data().universal_part();
  stk::mesh::BucketVector const& buckets = bulkData_->get_buckets( stk::topology::NODE_RANK, s_universal );

  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
        ib != buckets.end(); ++ib ) {
    const stk::mesh::Bucket & b = **ib;
    const stk::mesh::Bucket::size_type length = b.size();
    stk::mesh::EntityId *naluGlobalIds = stk::mesh::field_data(*naluGlobalId_, b);
    
    for ( stk::mesh::Bucket::size_type k = 0; k < length; ++k ) {
      naluGlobalIds[k] = bulkData_->identifier(b[k]);
    }
  }
}

//--------------------------------------------------------------------------
//-------- populate_boundary_data ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::populate_boundary_data()
{
  // realm first
  for ( size_t k = 0; k < bcDataAlg_.size(); ++k ) {
    bcDataAlg_[k]->execute();
  }
  equationSystems_.populate_boundary_data();
}

//--------------------------------------------------------------------------
//-------- output_banner ---------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::output_banner()
{
  if ( hasFluids_ )
    NaluEnv::self().naluOutputP0() << " Max Courant: " << maxCourant_ << " Max Reynolds: " << maxReynolds_ << " (" << name_ << ")" << std::endl;
}

//--------------------------------------------------------------------------
//-------- check_job -------------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::check_job(bool get_node_count)
{
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "Realm Review:              " << name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "===========================" << std::endl;

  // set number of nodes, check job run size
  size_t l_nodeCount = ioBroker_->get_input_ioss_region()->get_property("node_count").get_int();
  size_t g_nodeCount = 0;
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &l_nodeCount, &g_nodeCount, 1);
  NaluEnv::self().naluOutputP0() << "Node count from meta data = " << g_nodeCount << std::endl;
}

//--------------------------------------------------------------------------
//-------- dump_simulation_time --------------------------------------------
//--------------------------------------------------------------------------
void
  Realm::dump_simulation_time()
{
  NaluEnv::self().naluOutputP0() << std::endl;
  NaluEnv::self().naluOutputP0() << "-------------------------------- " << std::endl;
  NaluEnv::self().naluOutputP0() << "Begin Timer Overview for Realm: " << name_ << std::endl;
  NaluEnv::self().naluOutputP0() << "-------------------------------- " << std::endl;

  // equation system time
  equationSystems_.dump_eq_time();

  const int nprocs = NaluEnv::self().parallel_size();

  // common
  const unsigned ntimers = 6;
  double total_time[ntimers] = {timerCreateMesh_, timerOutputFields_, timerInitializeEqs_, 
                                timerPropertyEval_, timerPopulateMesh_, timerPopulateFieldData_ };
  double g_min_time[ntimers] = {}, g_max_time[ntimers] = {}, g_total_time[ntimers] = {};

  // get min, max and sum over processes
  stk::all_reduce_min(NaluEnv::self().parallel_comm(), &total_time[0], &g_min_time[0], ntimers);
  stk::all_reduce_max(NaluEnv::self().parallel_comm(), &total_time[0], &g_max_time[0], ntimers);
  stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &total_time[0], &g_total_time[0], ntimers);

  NaluEnv::self().naluOutputP0() << "Timing for IO: " << std::endl;
  NaluEnv::self().naluOutputP0() << "   io create mesh --  " << " \tavg: " << g_total_time[0]/double(nprocs)
                  << " \tmin: " << g_min_time[0] << " \tmax: " << g_max_time[0] << std::endl;
  NaluEnv::self().naluOutputP0() << " io output fields --  " << " \tavg: " << g_total_time[1]/double(nprocs)
                  << " \tmin: " << g_min_time[1] << " \tmax: " << g_max_time[1] << std::endl;
  NaluEnv::self().naluOutputP0() << " io populate mesh --  " << " \tavg: " << g_total_time[4]/double(nprocs)
                  << " \tmin: " << g_min_time[4] << " \tmax: " << g_max_time[4] << std::endl;
  NaluEnv::self().naluOutputP0() << " io populate fd   --  " << " \tavg: " << g_total_time[5]/double(nprocs)
                  << " \tmin: " << g_min_time[5] << " \tmax: " << g_max_time[5] << std::endl;
  NaluEnv::self().naluOutputP0() << "Timing for connectivity/finalize lysys: " << std::endl;
  NaluEnv::self().naluOutputP0() << "         eqs init --  " << " \tavg: " << g_total_time[2]/double(nprocs)
                  << " \tmin: " << g_min_time[2] << " \tmax: " << g_max_time[2] << std::endl;

  NaluEnv::self().naluOutputP0() << "Timing for property evaluation:         " << std::endl;
  NaluEnv::self().naluOutputP0() << "            props --  " << " \tavg: " << g_total_time[3]/double(nprocs)
                  << " \tmin: " << g_min_time[3] << " \tmax: " << g_max_time[3] << std::endl;

  // now edge creation; if applicable
  if ( realmUsesEdges_ ) {
    double g_total_edge = 0.0, g_min_edge = 0.0, g_max_edge = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &timerCreateEdges_, &g_min_edge, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &timerCreateEdges_, &g_max_edge, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &timerCreateEdges_, &g_total_edge, 1);

    NaluEnv::self().naluOutputP0() << "Timing for Edge: " << std::endl;
    NaluEnv::self().naluOutputP0() << "    edge creation --  " << " \tavg: " << g_total_edge/double(nprocs)
                    << " \tmin: " << g_min_edge << " \tmax: " << g_max_edge << std::endl;
  }

  // periodic
  if ( hasPeriodic_ ){
    double periodicSearchTime = periodicManager_->get_search_time();
    double g_minPeriodicSearchTime = 0.0, g_maxPeriodicSearchTime = 0.0, g_periodicSearchTime = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &periodicSearchTime, &g_minPeriodicSearchTime, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &periodicSearchTime, &g_maxPeriodicSearchTime, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &periodicSearchTime, &g_periodicSearchTime, 1);

    NaluEnv::self().naluOutputP0() << "Timing for Periodic: " << std::endl;
    NaluEnv::self().naluOutputP0() << "           search --  " << " \tavg: " << g_periodicSearchTime/double(nprocs)
                     << " \tmin: " << g_minPeriodicSearchTime << " \tmax: " << g_maxPeriodicSearchTime << std::endl;
  }

  // transfer
  if ( hasMultiPhysicsTransfer_ || hasInitializationTransfer_ || hasIoTransfer_ || hasExternalDataTransfer_ ) {
    double totalXfer[2] = {timerTransferSearch_, timerTransferExecute_};
    double g_totalXfer[2] = {}, g_minXfer[2] = {}, g_maxXfer[2] = {};
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &totalXfer[0], &g_minXfer[0], 2);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &totalXfer[0], &g_maxXfer[0], 2);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &totalXfer[0], &g_totalXfer[0], 2);

    NaluEnv::self().naluOutputP0() << "Timing for Tranfer (fromRealm):    " << std::endl;
    NaluEnv::self().naluOutputP0() << "           search --  " << " \tavg: " << g_totalXfer[0]/double(nprocs)
                                   << " \tmin: " << g_minXfer[0] << " \tmax: " << g_maxXfer[0] << std::endl;
    NaluEnv::self().naluOutputP0() << "          execute --  " << " \tavg: " << g_totalXfer[1]/double(nprocs)
                                   << " \tmin: " << g_minXfer[1] << " \tmax: " << g_maxXfer[1] << std::endl;
  }

  // skin mesh
  if ( checkForMissingBcs_ || hasOverset_ ) {
    double g_totalSkin = 0.0, g_minSkin= 0.0, g_maxSkin = 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &timerSkinMesh_, &g_minSkin, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &timerSkinMesh_, &g_maxSkin, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &timerSkinMesh_, &g_totalSkin, 1);
    
    NaluEnv::self().naluOutputP0() << "Timing for skin_mesh :    " << std::endl;    
    NaluEnv::self().naluOutputP0() << "        skin_mesh --  " << " \tavg: " << g_totalSkin/double(nprocs)
                                   << " \tmin: " << g_minSkin << " \tmax: " << g_maxSkin << std::endl;
  }

  // consolidated sort
  if (solutionOptions_->useConsolidatedSolverAlg_ ) {
    double g_totalSort= 0.0, g_minSort= 0.0, g_maxSort= 0.0;
    stk::all_reduce_min(NaluEnv::self().parallel_comm(), &timerSortExposedFace_, &g_minSort, 1);
    stk::all_reduce_max(NaluEnv::self().parallel_comm(), &timerSortExposedFace_, &g_maxSort, 1);
    stk::all_reduce_sum(NaluEnv::self().parallel_comm(), &timerSortExposedFace_, &g_totalSort, 1);
    
    NaluEnv::self().naluOutputP0() << "Timing for sort_mesh: " << std::endl;
    NaluEnv::self().naluOutputP0() << "       sort_mesh  -- " << " \tavg: " << g_totalSort/double(nprocs)
                                   << " \tmin: " << g_minSort<< " \tmax: " << g_maxSort<< std::endl;
  }

  NaluEnv::self().naluOutputP0() << std::endl;
}

//--------------------------------------------------------------------------
//-------- provide_mean_norm -----------------------------------------------
//--------------------------------------------------------------------------
double
Realm::provide_mean_norm()
{
  return equationSystems_.provide_mean_system_norm();
}

//--------------------------------------------------------------------------
//-------- get_hybrid_factor -----------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_hybrid_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->hybridDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->hybridMap_.find(dofName);
  if (iter != solutionOptions_->hybridMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_alpha_factor ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_alpha_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->alphaDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->alphaMap_.find(dofName);
  if (iter != solutionOptions_->alphaMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_alpha_upw_factor --------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_alpha_upw_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->alphaUpwDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->alphaUpwMap_.find(dofName);
  if (iter != solutionOptions_->alphaUpwMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_upw_factor --------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_upw_factor(
      const std::string dofName )
{
  double factor = solutionOptions_->upwDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->upwMap_.find(dofName);
  if (iter != solutionOptions_->upwMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- primitive_uses_limiter ------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::primitive_uses_limiter(
  const std::string dofName )
{
  bool usesIt = false;
  std::map<std::string, bool>::const_iterator iter
    = solutionOptions_->limiterMap_.find(dofName);
  if (iter != solutionOptions_->limiterMap_.end()) {
    usesIt = (*iter).second;
  }
  return usesIt;
}

//--------------------------------------------------------------------------
//-------- get_lam_schmidt -------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_lam_schmidt(
  const std::string dofName )
{
  double factor = solutionOptions_->lamScDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->lamScMap_.find(dofName);
  if (iter != solutionOptions_->lamScMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_lam_prandtl -------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_lam_prandtl(
  const std::string dofName, bool &prProvided )
{
  double factor = 1.0;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->lamPrMap_.find(dofName);
  if (iter != solutionOptions_->lamPrMap_.end()) {
    factor = (*iter).second;
    prProvided = true;
  }
  else {
    prProvided = false;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_turb_schmidt ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_turb_schmidt(
  const std::string dofName )
{
  double factor = solutionOptions_->turbScDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->turbScMap_.find(dofName);
  if (iter != solutionOptions_->turbScMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_turb_prandtl ------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_turb_prandtl(
  const std::string dofName )
{
  double factor = solutionOptions_->turbPrDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->turbPrMap_.find(dofName);
  if (iter != solutionOptions_->turbPrMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_noc_usage ---------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_noc_usage(
  const std::string dofName )
{
  bool factor = solutionOptions_->get_noc_usage(dofName);
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_shifted_grad_op ---------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_shifted_grad_op(
  const std::string dofName )
{
  bool factor = solutionOptions_->shiftedGradOpDefault_;
  std::map<std::string, bool>::const_iterator iter
    = solutionOptions_->shiftedGradOpMap_.find(dofName);
  if (iter != solutionOptions_->shiftedGradOpMap_.end()) {
    factor = (*iter).second;
  }
  return factor;
}

//--------------------------------------------------------------------------
//-------- get_skew_symmetric ----------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_skew_symmetric(
  const std::string dofName )
{
  return solutionOptions_->get_skew_symmetric(dofName);
}

//--------------------------------------------------------------------------
//-------- get_tanh_functional_form ----------------------------------------
//--------------------------------------------------------------------------
std::string
Realm::get_tanh_functional_form(
  const std::string dofName )
{
  std::string tanhForm = solutionOptions_->tanhFormDefault_;
  std::map<std::string, std::string>::const_iterator iter
    = solutionOptions_->tanhFormMap_.find(dofName);
  if (iter != solutionOptions_->tanhFormMap_.end()) {
    tanhForm = (*iter).second;
  }
  return tanhForm;
}

//--------------------------------------------------------------------------
//-------- get_tanh_trans --------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_tanh_trans(
  const std::string dofName )
{
  double tanhTrans = solutionOptions_->tanhTransDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->tanhTransMap_.find(dofName);
  if (iter != solutionOptions_->tanhTransMap_.end()) {
    tanhTrans = (*iter).second;
  }
  return tanhTrans;
}

//--------------------------------------------------------------------------
//-------- get_tanh_width --------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_tanh_width(
  const std::string dofName )
{
  double tanhWidth = solutionOptions_->tanhWidthDefault_;
  std::map<std::string, double>::const_iterator iter
    = solutionOptions_->tanhWidthMap_.find(dofName);
  if (iter != solutionOptions_->tanhWidthMap_.end()) {
    tanhWidth = (*iter).second;
  }
  return tanhWidth;
}

//--------------------------------------------------------------------------
//-------- get_consistent_mass_matrix_png ----------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_consistent_mass_matrix_png(
  const std::string dofName )
{
  bool cmmPng = solutionOptions_->consistentMMPngDefault_;
  std::map<std::string, bool>::const_iterator iter
    = solutionOptions_->consistentMassMatrixPngMap_.find(dofName);
  if (iter != solutionOptions_->consistentMassMatrixPngMap_.end()) {
    cmmPng = (*iter).second;
  }
  return cmmPng;
}

//--------------------------------------------------------------------------
//-------- get_divU --------------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_divU()
{
  return solutionOptions_->includeDivU_;
}

//--------------------------------------------------------------------------
//-------- get_cvfem_shifted_mdot ------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_cvfem_shifted_mdot()
{
  return solutionOptions_->cvfemShiftMdot_;
}

//--------------------------------------------------------------------------
//-------- get_cvfem_reduced_sens_poisson ----------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_cvfem_reduced_sens_poisson()
{
  return solutionOptions_->cvfemReducedSensPoisson_;
}

//--------------------------------------------------------------------------
//-------- has_nc_gauss_labatto_quadrature ---------------------------------
//--------------------------------------------------------------------------
bool
Realm::has_nc_gauss_labatto_quadrature()
{
  return solutionOptions_->ncAlgGaussLabatto_;
}

//--------------------------------------------------------------------------
//-------- get_nc_alg_upwind_advection -------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_nc_alg_upwind_advection()
{
  return solutionOptions_->ncAlgUpwindAdvection_;
}

//--------------------------------------------------------------------------
//-------- get_nc_alg_include_pstab ----------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_nc_alg_include_pstab()
{
  return solutionOptions_->ncAlgIncludePstab_;
}

//--------------------------------------------------------------------------
//-------- get_nc_alg_current_normal ---------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_nc_alg_current_normal()
{
  return solutionOptions_->ncAlgCurrentNormal_;
}

//--------------------------------------------------------------------------
//-------- get_material_prop_eval ------------------------------------------
//--------------------------------------------------------------------------
void
Realm::get_material_prop_eval(
  const PropertyIdentifier thePropID,
  std::vector<PropertyEvaluator*> &propEvalVec)
{
  for ( size_t i = 0; i < materialPropertys_.materialPropertyVector_.size(); ++i ) {
    PropertyEvaluator *thePropEval = NULL;
    std::map<PropertyIdentifier, PropertyEvaluator*>::const_iterator iter
      = materialPropertys_.materialPropertyVector_[i]->propertyEvalMap_.find(thePropID);
    if (iter != materialPropertys_.materialPropertyVector_[i]->propertyEvalMap_.end()) {
      thePropEval = (*iter).second;
    }
    propEvalVec.push_back(thePropEval);
  }
}

//--------------------------------------------------------------------------
//-------- is_turbulent ----------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::is_turbulent()
{
  return solutionOptions_->isTurbulent_;
}

//--------------------------------------------------------------------------
//-------- is_turbulent ----------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::is_turbulent( bool isIt )
{
  isTurbulent_ = isIt;
  solutionOptions_->isTurbulent_ = isIt;
}

//--------------------------------------------------------------------------
//-------- needs_enthalpy --------------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::needs_enthalpy()
{
  return needsEnthalpy_;
}

//--------------------------------------------------------------------------
//-------- needs_enthalpy --------------------------------------------------
//--------------------------------------------------------------------------
void
Realm::needs_enthalpy( bool needsEnthalpy )
{
  needsEnthalpy_ = needsEnthalpy;
}

//--------------------------------------------------------------------------
//-------- number_of_states ------------------------------------------------
//--------------------------------------------------------------------------
int
Realm::number_of_states()
{
  const int numStates = (timeIntegrator_->secondOrderTimeAccurate_) ? 3 : 2;
  return numStates;
}

//--------------------------------------------------------------------------
//-------- name ------------------------------------------------------------
//--------------------------------------------------------------------------
std::string
Realm::name()
{
  return name_;
}

//--------------------------------------------------------------------------
//-------- augment_transfer_vector -----------------------------------------
//--------------------------------------------------------------------------
void
Realm::augment_transfer_vector(Transfer *transfer, const std::string transferObjective, Realm *toRealm)
{
  if ( transferObjective == "multi_physics" ) {
    multiPhysicsTransferVec_.push_back(transfer);
    hasMultiPhysicsTransfer_ = true; 
  }
  else if ( transferObjective == "initialization" ) {
    initializationTransferVec_.push_back(transfer);
    hasInitializationTransfer_ = true;
  }
  else if ( transferObjective == "input_output" ) {
    toRealm->ioTransferVec_.push_back(transfer);
    toRealm->hasIoTransfer_ = true;
  }
  else if ( transferObjective == "external_data" ) {
    toRealm->externalDataTransferVec_.push_back(transfer);
    toRealm->hasExternalDataTransfer_ = true;
  }
  else { 
    throw std::runtime_error("Real::augment_transfer_vector: Error, none supported transfer objective: " + transferObjective);
  }
}

//--------------------------------------------------------------------------
//-------- process_multi_physics_transfer ----------------------------------
//--------------------------------------------------------------------------
void
Realm::process_multi_physics_transfer(
  const bool forcedXfer)
{
  if ( !hasMultiPhysicsTransfer_ )
    return;
  
  // only process if an active time step, however, calling class can enforce
  if ( active_time_step() || forcedXfer ) {
    double timeXfer = -NaluEnv::self().nalu_time();
    std::vector<Transfer *>::iterator ii;
    for( ii=multiPhysicsTransferVec_.begin(); ii!=multiPhysicsTransferVec_.end(); ++ii )
      (*ii)->execute();
    timeXfer += NaluEnv::self().nalu_time();
    timerTransferExecute_ += timeXfer;
  }
}

//--------------------------------------------------------------------------
//-------- process_initialization_transfer ---------------------------------
//--------------------------------------------------------------------------
void
Realm::process_initialization_transfer()
{
  if ( !hasInitializationTransfer_ )
    return;

  double timeXfer = -NaluEnv::self().nalu_time();
  std::vector<Transfer *>::iterator ii;
  for( ii=initializationTransferVec_.begin(); ii!=initializationTransferVec_.end(); ++ii ) {
    (*ii)->execute();
  }
  timeXfer += NaluEnv::self().nalu_time();
  timerTransferExecute_ += timeXfer;
}

//--------------------------------------------------------------------------
//-------- process_io_transfer ---------------------------------------------
//--------------------------------------------------------------------------
void
Realm::process_io_transfer()
{
  if ( !hasIoTransfer_ )
    return;

  double timeXfer = -NaluEnv::self().nalu_time();
  // only do at an IO step
  const int timeStepCount = get_time_step_count();
  const bool isOutput = (timeStepCount % outputInfo_->outputFreq_) == 0;
  if ( isOutput ) {
    std::vector<Transfer *>::iterator ii;
    for( ii=ioTransferVec_.begin(); ii!=ioTransferVec_.end(); ++ii )
      (*ii)->execute();
  }
  timeXfer += NaluEnv::self().nalu_time();
  timerTransferExecute_ += timeXfer;
}

//--------------------------------------------------------------------------
//-------- process_external_data_transfer ----------------------------------
//--------------------------------------------------------------------------
void
Realm::process_external_data_transfer()
{
  if ( !hasExternalDataTransfer_ )
    return;

  double timeXfer = -NaluEnv::self().nalu_time();
  std::vector<Transfer *>::iterator ii;
  for( ii=externalDataTransferVec_.begin(); ii!=externalDataTransferVec_.end(); ++ii )
    (*ii)->execute();
  timeXfer += NaluEnv::self().nalu_time();
  timerTransferExecute_ += timeXfer;
}

//--------------------------------------------------------------------------
//-------- post_converged_work ---------------------------------------------
//--------------------------------------------------------------------------
void
Realm::post_converged_work()
{
  equationSystems_.post_converged_work();
}

//--------------------------------------------------------------------------
//-------- part_name(std::string) ------------------------------------------
//--------------------------------------------------------------------------
std::string
Realm::physics_part_name(std::string name) const
{
  return name;
}

std::vector<std::string>
Realm::physics_part_names(std::vector<std::string> names) const
{
  return names;
}

//--------------------------------------------------------------------------
//-------- get_current_time() ----------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_current_time()
{
  return timeIntegrator_->get_current_time();
}

//--------------------------------------------------------------------------
//-------- get_time_step() -------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_time_step()
{
  return timeIntegrator_->get_time_step();
}

double 
Realm::get_time_step_from_file() {
  return timeIntegrator_->get_time_step_from_file();
}

bool 
Realm::get_is_fixed_time_step() {
  return timeIntegrator_->get_is_fixed_time_step();
}

bool 
Realm::get_is_terminate_based_on_time() {
  return timeIntegrator_->get_is_terminate_based_on_time();
}

double 
Realm::get_total_sim_time() {
  return timeIntegrator_->get_total_sim_time();
}

int
Realm::get_max_time_step_count() {
  return timeIntegrator_->get_max_time_step_count();
}

//--------------------------------------------------------------------------
//-------- get_gamma1() ----------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_gamma1()
{
  return timeIntegrator_->get_gamma1();
}

//--------------------------------------------------------------------------
//-------- get_gamma2() ----------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_gamma2()
{
  return timeIntegrator_->get_gamma2();
}

//--------------------------------------------------------------------------
//-------- get_gamma3() ----------------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_gamma3()
{
  return timeIntegrator_->get_gamma3();
}

//--------------------------------------------------------------------------
//-------- get_time_step_count() -------------------------------------------
//--------------------------------------------------------------------------
int
Realm::get_time_step_count() const
{
  return timeIntegrator_->get_time_step_count();
}

//--------------------------------------------------------------------------
//-------- restarted_simulation() ------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::restarted_simulation()
{
  return outputInfo_->activateRestart_ ;
}

//--------------------------------------------------------------------------
//-------- support_inconsistent_restart() ----------------------------------
//--------------------------------------------------------------------------
bool
Realm::support_inconsistent_restart()
{
  return supportInconsistentRestart_ ;
}

//--------------------------------------------------------------------------
//-------- get_stefan_boltzmann() ------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_stefan_boltzmann()
{
  return solutionOptions_->stefanBoltzmann_;
}

//--------------------------------------------------------------------------
//-------- get_turb_model_constant() ---------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_turb_model_constant(
   const TurbulenceModelConstant turbModelEnum)
{
  std::map<TurbulenceModelConstant, double>::iterator it
    = solutionOptions_->turbModelConstantMap_.find(turbModelEnum);
  if ( it != solutionOptions_->turbModelConstantMap_.end() ) {
    return it->second;
  }
  else {
    throw std::runtime_error("unknown (not found) turbulence model constant");
  }
}

//--------------------------------------------------------------------------
//-------- get_buckets() ---------------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::BucketVector const& Realm::get_buckets( 
  stk::mesh::EntityRank rank,
  const stk::mesh::Selector & selector) const
{
  return bulkData_->get_buckets(rank, selector);
}

//--------------------------------------------------------------------------
//-------- bulk_data() -----------------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::BulkData &
Realm::bulk_data()
{
  return *bulkData_;
}

const stk::mesh::BulkData &
Realm::bulk_data() const
{
  return *bulkData_;
}

//--------------------------------------------------------------------------
//-------- meta_data() -----------------------------------------------------
//--------------------------------------------------------------------------
stk::mesh::MetaData &
Realm::meta_data()
{
  return bulkData_->mesh_meta_data();
}

const stk::mesh::MetaData &
Realm::meta_data() const
{
  return bulkData_->mesh_meta_data();
}

//--------------------------------------------------------------------------
//-------- get_activate_aura() ---------------------------------------------
//--------------------------------------------------------------------------
bool
Realm::get_activate_aura()
{
  return activateAura_;
}

//--------------------------------------------------------------------------
//-------- get_inactive_selector() -----------------------------------------
//--------------------------------------------------------------------------
stk::mesh::Selector
Realm::get_inactive_selector()
{
  // accumulate inactive parts relative to the universal part
  if ( hasOverset_ )
    throw std::runtime_error("Realm::get_tanh_blending");

  // find parts that are neither universanl, interior, or bcs
  stk::mesh::Selector nothing;
  stk::mesh::Selector otherInactiveSelector = (
    meta_data().universal_part()
    & !(stk::mesh::selectUnion(interiorPartVec_))
    & !(stk::mesh::selectUnion(bcPartVec_)));

  if (interiorPartVec_.empty() && bcPartVec_.empty()) {
    otherInactiveSelector = nothing;
  }

  return otherInactiveSelector;
}

//--------------------------------------------------------------------------
//-------- push_equation_to_systems() --------------------------------------
//--------------------------------------------------------------------------
void
Realm::push_equation_to_systems(
  EquationSystem *eqSystem)
{
  equationSystems_.equationSystemVector_.push_back(eqSystem);
}

//--------------------------------------------------------------------------
//-------- get_physics_target_names() --------------------------------------
//--------------------------------------------------------------------------
const std::vector<std::string> &
Realm::get_physics_target_names()
{
  // in the future, possibly check for more advanced names;
  // for now, material props holds this'
  return materialPropertys_.targetNames_;
}

//--------------------------------------------------------------------------
//-------- get_tanh_blending() ---------------------------------------------
//--------------------------------------------------------------------------
double
Realm::get_tanh_blending(
 const std::string dofName)
{
  throw std::runtime_error("Realm::get_tanh_blending");
}

//--------------------------------------------------------------------------
//-------- balance_nodes() -------------------------------------------------
//--------------------------------------------------------------------------
void Realm::balance_nodes()
{
  throw std::runtime_error("Realm::balanced_nodes");
}

//--------------------------------------------------------------------------
//-------- get_quad_type() -------------------------------------------------
//--------------------------------------------------------------------------
std::string Realm::get_quad_type() const
{
  STK_ThrowRequire(solutionOptions_ != nullptr);
  return solutionOptions_->quadType_;
}

} // namespace nalu
} // namespace Sierra
