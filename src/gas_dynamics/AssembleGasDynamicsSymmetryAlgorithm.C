/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <gas_dynamics/AssembleGasDynamicsSymmetryAlgorithm.h>
#include <Realm.h>
#include <element/Element.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// stk_util
#include <stk_util/parallel/ParallelReduce.hpp>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleGasDynamicsSymmetryAlgorithm - assembles RHS for gas dynamics; sym
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleGasDynamicsSymmetryAlgorithm::AssembleGasDynamicsSymmetryAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  ScalarFieldType *pressure,
  GenericFieldType *rhsGasDyn)
  : Algorithm(realm, part),
    pressure_(pressure),
    rhsGasDyn_(rhsGasDyn),
    coordinates_(NULL)
{
  // save off mising fields
  stk::mesh::MetaData & metaData = realm_.meta_data();
  coordinates_ = metaData.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleGasDynamicsSymmetryAlgorithm::execute()
{
  stk::mesh::MetaData & metaData = realm_.meta_data();

  // sizes
  const int nDim = metaData.spatial_dimension();

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = metaData.locally_owned_part()
    & stk::mesh::selectUnion(partVec_) 
    & !(realm_.get_inactive_selector());

  //===========================================================
  // assemble edge-based flux operator to the node
  //===========================================================

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( metaData.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const stk::mesh::Bucket::size_type length   = b.size();

    // extract master element
    Element *meFC = sierra::nalu::ElementRepo::get_surface_element(b.topology());

    // extract master element specifics
    const int nodesPerElement = meFC->nodesPerElement_;
    const int numScsIp = meFC->numIntPoints_;
    const int *ipNodeMap = meFC->ipNodeMap();

    // scratch fields to gather
    std::vector<double > ws_coordinates(nodesPerElement*nDim);
    
    // scatch fields for geometry
    std::vector<double > ws_scs_areav(numScsIp*nDim);

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // face node relations for nodal gather
      stk::mesh::Entity const * face_node_rels = b.begin_nodes(k);

      //===============================================
      // gather nodal data; this is how we do it now..
      //===============================================
      int num_nodes = b.num_nodes(k);
      for ( int ni = 0; ni < num_nodes; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        double * coords = stk::mesh::field_data(*coordinates_, node);
        const int offSet = ni*nDim;
        for ( int j=0; j < nDim; ++j ) {
          ws_coordinates[offSet+j] = coords[j];
        }
      }

      // compute scs integration point areavec
      double scs_error = 0.0;
      meFC->determinant(1, &ws_coordinates[0], &ws_scs_areav[0], &scs_error);

      // start assembly
      for ( int ip = 0; ip < numScsIp; ++ip ) {

        const int ipNdim = ip*nDim;
        
        // nearest node
        const int nn = ipNodeMap[ip];
        
        stk::mesh::Entity nodeNN = face_node_rels[nn];

        // pointer to fields to assemble
        double *rhsGasDynNN = stk::mesh::field_data(*rhsGasDyn_, nodeNN);
        
        // suplemental
        double pressureNN = *stk::mesh::field_data(*pressure_, nodeNN);
        
        for ( int i = 0; i < nDim; ++i ) {
          rhsGasDynNN[i] -= pressureNN*ws_scs_areav[ipNdim+i];
        }
      }
    }
  }
}

} // namespace nalu
} // namespace Sierra
