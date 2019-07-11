// Ourselves:
#include <vertex_module.h>

// Standard library:
#include <iostream>
#include <stdexcept>

// geomtools
#include "bayeux/geomtools/line_3d.h"
#include "bayeux/geomtools/helix_3d.h"
#include "bayeux/geomtools/plane.h"
#include "bayeux/geomtools/blur_spot.h"

// falaise
#include <falaise/snemo/datamodels/data_model.h>
#include <falaise/snemo/datamodels/tracker_trajectory.h>
#include <falaise/snemo/datamodels/tracker_trajectory_solution.h>
#include <falaise/snemo/datamodels/tracker_trajectory_data.h>
#include "falaise/snemo/datamodels/calibrated_calorimeter_hit.h"
#include "falaise/snemo/datamodels/base_trajectory_pattern.h"

// Registration instantiation macro :
DPP_MODULE_REGISTRATION_IMPLEMENT(vertex_module,
				  "vertex_module")


void vertex_module::initialize(const datatools::properties  & setup_,
					       datatools::service_manager   & service_manager_,
					       dpp::module_handle_dict_type & /* module_dict_ */)
{
  DT_THROW_IF (this->is_initialized(),
	       std::logic_error,
	       "Module 'Vertex Extrapolator' is already initialized ! ");
  
  dpp::base_module::_common_initialize(setup_);
 
  // Look for services
  if (service_manager_.has("geometry"))
  {
    const geomtools::geometry_service& GS = service_manager_.get<geomtools::geometry_service> ("geometry");

    // initialize geometry manager
    //    std::cout << "Initialize geo manager " << std::endl;
    geometry_manager_ = &GS.get_geom_manager();
    DT_THROW_IF(!geometry_manager_,
                std::runtime_error,
                "Null pointer to geometry manager return by geometry_service");
  }

  // check the label
  _TTD_label_ = snemo::datamodel::data_info::default_tracker_trajectory_data_label();
  _TCD_label_ = snemo::datamodel::data_info::default_tracker_clustering_data_label();

  eventCounter = 0;
  this->_set_initialized(true);

  return;
}

void vertex_module::reset()
{
  DT_THROW_IF (! this->is_initialized(),
	       std::logic_error,
	       "Module 'Vertex Extrapolator' is not initialized !");
  this->_set_initialized(false);
  
  // clean up
  _TTD_label_.clear();
  _TCD_label_.clear();

  eventCounter = 0;
  std::cout << "Vertex Extrapolator Module finished." << std::endl;
  return;
}

// Constructor :
vertex_module::vertex_module(datatools::logger::priority logging_priority_)
  : dpp::base_module(logging_priority_)
{
}

// Destructor :
vertex_module::~vertex_module()
{
  // MUST reset module at destruction
  if (this->is_initialized()) reset();
}

// Processing :
dpp::base_module::process_status vertex_module::process(datatools::things & data_record_)
{
  DT_THROW_IF (! this->is_initialized(), std::logic_error,
	       "Module 'Vertex Extrapolator' is not initialized !");
  
  ///////////////////////////////////
  // Check tracker clustering data //
  ///////////////////////////////////

  snemo::datamodel::tracker_clustering_data * ptr_cluster_data = 0;
  if (!data_record_.has(_TCD_label_)) {
    std::cerr << "failed to grab TCD bank " << std::endl;
    return dpp::base_module::PROCESS_INVALID;
  }
  else {
    ptr_clustering_data = &(data_record_.grab<snemo::datamodel::tracker_clustering_data>(_TCD_label_));
    snemo::datamodel::tracker_clustering_data & the_clustering_data = *ptr_clustering_data;
    if (! the_clustering_data.has_solutions()) {
      std::cerr << "TCD bank empty, no clustering data." << std::endl;
      return dpp::base_module::PROCESS_INVALID;
    }
  }

  ///////////////////////////////////
  // Check tracker trajectory data // not used
  ///////////////////////////////////

  snemo::datamodel::tracker_trajectory_data * ptr_trajectory_data = 0;
  if (!data_record_.has(_TTD_label_)) {
    std::cerr << "failed to grab TTD bank " << std::endl;
    return dpp::base_module::PROCESS_INVALID;
  }
  else {
    ptr_trajectory_data = &(data_record_.grab<snemo::datamodel::tracker_trajectory_data>(_TTD_label_));
    snemo::datamodel::tracker_trajectory_data & the_trajectory_data = *ptr_trajectory_data;
    if (! the_trajectory_data.has_solutions()) {
      std::cerr << "TTD bank empty, no trajectories." << std::endl;
      return dpp::base_module::PROCESS_INVALID;
    }
  }

  /********************
   * Process the data *
   ********************/
  
  // Main processing method :
  // Process the vertex extrapolator:
  namespace sdm = snemo::datamodel;

  // Process clusters hits for fitting
  std::cout << "In process: event counter = " << eventCounter << std::endl;

  // get all cluster solutions
  const sdm::tracker_clustering_data::solution_col_type& all_cls_solutions = ptr_cluster_data->get_solutions();

  // get all trajectory solutions
  // somehow, making line, etc. _trj collection

  // last cls solution, normally only one entry anyway
  const sdm::tracker_clustering_solution::cluster_col_type &cls_defaults = all_cls_solutions.back().get().get_clusters();
  
  std::vector<VertexInfo> all_info;
  VertexInfo info;
  for (auto entry : cls_defaults) { // get sdm::tracker_cluster handle
    // check all clusters for wire_candidate
    const sdm::calibrated_tracker_hit::collection_type & gg_hits_col = entry.get().get_hits();

    VertexInfo info = check_on_wire(gg_hits_col);
    info.clsid = entry.get().get_cluster_id();
    all_info.push_back(info);
    // into a function to check on wire_candidate with vertex pair<bool, bool> 
    // for <calo extrapolate, foil extrapolate> = if any one pair entry false then wire candidate is true
    // for (auto hit_handle : gg_hits_col) {
    //   const sdm::calibrated_tracker_hit & hit = hit_handle.get();
  }
  // ready to extrapolate, work with loop over all trj containers, check on wire_candidate via clsid
  
  eventCounter++;
  return dpp::base_module::PROCESS_SUCCESS;
}



VertexInfo vertex_module::check_on_wire(const sdm::calibrated_tracker_hit::collection_type & data) {

}
