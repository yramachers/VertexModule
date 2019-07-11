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
#include <falaise/snemo/datamodels/calibrated_data.h>

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
  // Check tracker trajectory data //
  ///////////////////////////////////

  bool preserve_former_output = true; // keep all
  
  // check if some 'tracker_trajectory_data' are available in the data model:
  snemo::datamodel::tracker_trajectory_data * ptr_trajectory_data = 0;
  if (! data_record_.has(_TTD_label_)) {
    ptr_trajectory_data = &(data_record_.add<snemo::datamodel::tracker_trajectory_data>(_TTD_label_));
  } else {
    ptr_trajectory_data = &(data_record_.grab<snemo::datamodel::tracker_trajectory_data>(_TTD_label_));
  }
  snemo::datamodel::tracker_trajectory_data & the_trajectory_data = *ptr_trajectory_data;
  if (the_trajectory_data.has_solutions()) 
    if (! preserve_former_output) 
      the_trajectory_data.reset();
  
  
  /********************
   * Process the data *
   ********************/
  
  // Main processing method :
  // Process the vertex extrapolator:
  namespace sdm = snemo::datamodel;

  // Process events for trajectory consolidation
  // make a trajectory solution
  sdm::tracker_trajectory_solution::handle_type htts(new sdm::tracker_trajectory_solution);
  the_trajectory_data.add_solution(htts, true);
  the_trajectory_data.grab_solutions().back().grab().set_solution_id(the_trajectory_data.get_number_of_solutions() - 1);
  sdm::tracker_trajectory_solution & trajectory_solution = the_trajectory_data.grab_solutions().back().grab(); // maybe store in here a bit

  // Process clusters hits for fitting
  std::cout << "In process: event counter = " << eventCounter << std::endl;

  // get all cluster solutions
  const sdm::tracker_clustering_data::solution_col_type& all_solutions = ptr_cluster_data->get_solutions();

  int nsol = 0;
  int ncl  = 0;
  for (auto entry : all_solutions) { 
    const sdm::tracker_clustering_solution::cluster_col_type &defaults = entry.get().get_clusters();
    std::cout << "cluster solution number: " << ++nsol << std::endl;

    for (auto cl_handle : defaults) {
      const sdm::calibrated_tracker_hit::collection_type & gg_hits_col = cl_handle.get().get_hits();
      ncl = cl_handle.get().get_cluster_id();
      std::cout << "cluster number: " << ncl << std::endl;

      rings.clear();
      std::cout << "number of gg hits: " << gg_hits_col.size() << std::endl;
      for (auto hit_handle : gg_hits_col) {
	// work with geiger hits as members of a given cluster
	const sdm::calibrated_tracker_hit & hit = hit_handle.get();
	ring.rerr   = hit.get_sigma_r();
	ring.zerr   = hit.get_sigma_z();
	ring.radius = hit.get_r();
	ring.wirex  = hit.get_x();
	ring.wirey  = hit.get_y();
	ring.zcoord = hit.get_z();
	mi.hitid  = hit.get_id();
	mi.side   = hit.get_side();
	mi.row    = hit.get_row();
	mi.column = hit.get_layer();
	th.clid = ncl;
	th.mi = mi;
	th.gr = ring;
	rings.push_back(th);
      }

      // ready to extrapolate
    }
  }

  eventCounter++;
  return dpp::base_module::PROCESS_SUCCESS;
}

