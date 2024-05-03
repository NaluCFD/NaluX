/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef PeriodicManager_h
#define PeriodicManager_h

//==============================================================================
// Includes and forwards
//==============================================================================

#include <FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Ghosting.hpp>

#include <stk_search/BoundingBox.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/IdentProc.hpp>

#include <vector>
#include <list>
#include <map>

namespace sierra {
namespace nalu {

class Realm;

typedef stk::search::IdentProc<stk::mesh::EntityKey,int> theEntityKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef std::pair<Sphere,theEntityKey> sphereBoundingBox;

class PeriodicManager {

 public:

  // constructor and destructor
  PeriodicManager(
    Realm & realm);

  ~PeriodicManager();

  void initialize_error_count();

  void add_periodic_pair(
    stk::mesh::Part* meshPartsMonarch,
    stk::mesh::Part* meshPartsSubject,
    const double &userSearchTolerance,
    const std::string &searchMethodName);

  void build_constraints();

  // holder for monarch += subject; subject = monarch
  void apply_constraints(
    stk::mesh::FieldBase *,
    const unsigned sizeOfField,
    const bool bypassFieldCheck,
    const bool addSubjects,
    const bool setSubjects);

  // find the max
  void apply_max_field(
    stk::mesh::FieldBase *,
    const unsigned &sizeOfField);

  void manage_ghosting_object();

  stk::mesh::Ghosting * get_ghosting_object();

  const stk::mesh::PartVector &get_subject_part_vector();

  double get_search_time();

// private:

  void augment_periodic_selector_pairs();

  void initialize_translation_vector();

  void determine_translation(
     stk::mesh::Selector monarchSelector,
     stk::mesh::Selector subjectSelector,
     std::vector<double> &translationVector,
     std::vector<double> &rotationVector);

  void remove_redundant_subject_nodes();

  void finalize_search();

  void populate_search_key_vec(
    stk::mesh::Selector monarchSelector,
    stk::mesh::Selector subjectSelector,
    std::vector<double> &translationVector,
    const stk::search::SearchMethod searchMethod);

  void error_check();

  void update_global_id_field();

  /* communicate periodicGhosting nodes */
  void
  periodic_parallel_communicate_field(
    stk::mesh::FieldBase *theField);

  /* communicate shared nodes and aura nodes */
  void
  parallel_communicate_field(
    stk::mesh::FieldBase *theField);

  Realm &realm_;

  /* manage tolerances; each block specifies a user tolerance */
  double searchTolerance_;

  stk::mesh::Ghosting *periodicGhosting_;
  const std::string ghostingName_;
  double timerSearch_;

  // algorithm to find/exclude points when M/S does not match
  int errorCount_;
  int maxErrorCount_;
  const double amplificationFactor_;

 public:
  // the data structures to hold monarch/subject information
  typedef std::pair<stk::mesh::Entity, stk::mesh::Entity> EntityPair;
  typedef std::pair<stk::mesh::Selector, stk::mesh::Selector> SelectorPair;
  typedef std::vector<std::pair<theEntityKey,theEntityKey> > SearchKeyVector;

  std::vector<int> ghostCommProcs_;

 private:

  // vector of monarch:subject selector pairs
  std::vector<SelectorPair> periodicSelectorPairs_;

  // vector of subject parts
  stk::mesh::PartVector subjectPartVector_;

  // vector of search types
  std::vector<stk::search::SearchMethod> searchMethodVec_;

  // translation and rotation
  std::vector<std::vector<double> > translationVector_;
  std::vector<std::vector<double> > rotationVector_;

  // vector of monarchEntity:subjectEntity
  std::vector<EntityPair> monarchSubjectCommunicator_;

  // culmination of all searches
  SearchKeyVector searchKeyVector_;

  void add_subject_to_monarch(
    stk::mesh::FieldBase *theField,
    const unsigned &sizeOfField,
    const bool &bypassFieldCheck);

  void set_subject_to_monarch(
    stk::mesh::FieldBase *theField,
    const unsigned &sizeOfField,
    const bool &bypassFieldCheck);

};

} // namespace nalu
} // namespace sierra

#endif
