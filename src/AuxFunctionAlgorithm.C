/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include "AuxFunctionAlgorithm.h"
#include "AuxFunction.h"
#include "FieldTypeDef.h"
#include "Realm.h"
#include "Simulation.h"

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/FieldParallel.hpp>

namespace sierra{
namespace nalu{

AuxFunctionAlgorithm::AuxFunctionAlgorithm(
  Realm & realm,
  stk::mesh::Part * part,
  stk::mesh::FieldBase * field,
  AuxFunction * auxFunction,
  stk::mesh::EntityRank entityRank,
  const bool parallelCommunicate,
  const bool periodicCommunicate)
  : Algorithm(realm, part),
    field_(field),
    auxFunction_(auxFunction),
    entityRank_(entityRank),
    parallelCommunicate_(parallelCommunicate),
    periodicCommunicate_(periodicCommunicate)
{
  // does nothing
}

AuxFunctionAlgorithm::~AuxFunctionAlgorithm() {
  // delete Aux
  delete auxFunction_;
}

void
AuxFunctionAlgorithm::execute()
{

  // make sure that partVec_ is size one
  STK_ThrowAssert( partVec_.size() == 1 );

  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const unsigned nDim = meta_data.spatial_dimension();
  const double time = realm_.get_current_time();
  VectorFieldType *coordinates = meta_data.get_field<double>(stk::topology::NODE_RANK, realm_.get_coordinates_name());

  auxFunction_->setup(time);

  stk::mesh::Selector selector = stk::mesh::selectUnion(partVec_) &
    stk::mesh::selectField(*field_);

  stk::mesh::BucketVector const& buckets =
    realm_.get_buckets( entityRank_, selector );

  unsigned fieldSizeSaved = 0;
  for ( stk::mesh::BucketVector::const_iterator ib = buckets.begin();
        ib != buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;
    const unsigned fieldSize = field_bytes_per_entity(*field_, b) / sizeof(double);
    // save off field size for possible nodal field operations
    fieldSizeSaved = fieldSize;

    const stk::mesh::Bucket::size_type length   = b.size();

    // FIXME: because coordinates are only defined at nodes, this actually
    //        only works for nodal fields. Hrmmm.
    const double * coords = stk::mesh::field_data( *coordinates, *b.begin() );
    double * fieldData = (double*) stk::mesh::field_data( *field_, *b.begin() );

    auxFunction_->evaluate(coords, time, nDim, length, fieldData, fieldSize);
  }

  // nodal fields may need to be parallel communicated (random, periodic)
  if ( entityRank_ == stk::topology::NODE_RANK ) {
    if ( parallelCommunicate_ ) {
      std::vector< const stk::mesh::FieldBase *> fieldVec(1, field_);
      stk::mesh::copy_owned_to_shared( realm_.bulk_data(), fieldVec);
      stk::mesh::communicate_field_data(realm_.bulk_data().aura_ghosting(), fieldVec);
    }
    if ( periodicCommunicate_ ) {
      // all nodal fields have the same field size; logics are : bypass, +=, =
      realm_.periodic_field_update(field_, fieldSizeSaved, true, false, true);
    }
  }

}

} // namespace nalu
} // namespace Sierra
