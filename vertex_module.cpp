// Ourselves:
#include <vertex_module.h>

// Standard library:
#include <iostream>
#include <stdexcept>
#include <string>

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
#include <falaise/snemo/geometry/calo_locator.h>
#include <falaise/snemo/geometry/gg_locator.h>
#include <falaise/snemo/geometry/gveto_locator.h>
#include <falaise/snemo/geometry/xcalo_locator.h>

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
  if (service_manager_.has("geometry")) {
    const geomtools::geometry_service& GS = service_manager_.get<geomtools::geometry_service> ("geometry");
    
    // initialize geometry manager
    //    std::cout << "Initialize geo manager " << std::endl;
    geometry_manager_ = &GS.get_geom_manager();
    DT_THROW_IF(!geometry_manager_,
		std::runtime_error,
		"Null pointer to geometry manager return by geometry_service");
  }
  
  // get locator plugin
  const geomtools::manager::plugins_dict_type &plugins = geometry_manager.get_plugins();
  for (geomtools::manager::plugins_dict_type::const_iterator ip = plugins.begin();
       ip != plugins.end(); ip++) {
    const std::string &plugin_name = ip->first;
    if (geometry_manager.is_plugin_a<snemo::geometry::locator_plugin>(plugin_name)) {
      DT_LOG_DEBUG(get_logging_priority(), "Find locator plugin with name = " << plugin_name);
      std::string locator_plugin_name = plugin_name;
      break;
    }
  }
  // Access to a given plugin by name and type :
  DT_THROW_IF(!geometry_manager.has_plugin(locator_plugin_name) ||
	      !geometry_manager.is_plugin_a<snemo::geometry::locator_plugin>(locator_plugin_name),
              std::logic_error, "Found no locator plugin named '" << locator_plugin_name << "'");
  locator_plugin_ = &geometry_manager.get_plugin<snemo::geometry::locator_plugin>(locator_plugin_name);


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


  /***************************
   * Calorimeter geometry info
   ***************************/
  
  // Set the calorimeter locators :
  const snemo::geometry::calo_locator  &calo_locator = _locator_plugin_->get_calo_locator();
  const snemo::geometry::xcalo_locator &xcalo_locator = _locator_plugin_->get_xcalo_locator();
  const snemo::geometry::gveto_locator &gveto_locator = _locator_plugin_->get_gveto_locator();

  int whichcalo = 0; // index position in vector
  int side = 0; // back
  std::vector<Plane> planes; // make a total of ten planes to potentially intersect

  // main calo walls
  const double xcalo_bd[2] = {calo_locator.get_wall_window_x(snemo::geometry::utils::SIDE_BACK),
                              calo_locator.get_wall_window_x(snemo::geometry::utils::SIDE_FRONT)};
  ROOT::Math::XYZPoint p(xcalo_bd[0],0.0,0.0); // negative x
  ROOT::Math::XYZVector norm(1.0,0.0,0.0); // looking in pos x direction
  Plane mcalo_negx;
  mcalo_negx.planeid = whichcalo;
  mcalo_negx.side = side;
  mcalo_negx.normal = norm;
  mcalo_negx.point = p;
  planes.push_back(mcalo_negx);
  
  ROOT::Math::XYZPoint p2(xcalo_bd[1],0.0,0.0); // positive x
  ROOT::Math::XYZVector norm2(-1.0,0.0,0.0); // looking in neg x direction
  Plane mcalo_posx;
  mcalo_posx.planeid = ++whichcalo;
  mcalo_posx.side = 1; // positive x
  mcalo_posx.normal = norm2;
  mcalo_posx.point = p2;
  planes.push_back(mcalo_posx);

  // xwalls part 1
  const double ycalo_bd_l[2] = {
    xcalo_locator.get_wall_window_y(side, snemo::geometry::xcalo_locator::WALL_LEFT),
    xcalo_locator.get_wall_window_y(side, snemo::geometry::xcalo_locator::WALL_RIGHT)};
  ROOT::Math::XYZPoint p3(0.0,ycalo_bd_l[0],0.0); // negative y
  ROOT::Math::XYZVector norm3(0.0,1.0,0.0); // looking in pos y direction
  ROOT::Math::XYZPoint p4(0.0,ycalo_bd_l[1],0.0); // positive y
  ROOT::Math::XYZVector norm4(0.0,-1.0,0.0); // looking in neg y direction
  Plane xwall_backl;
  xwall_backl.planeid = ++whichcalo;
  xwall_backl.side = side;
  xwall_backl.normal = norm3;
  xwall_backl.point = p3;

  Plane xwall_backr;
  xwall_backr.planeid = ++whichcalo;
  xwall_backr.side = side;
  xwall_backr.normal = norm4;
  xwall_backr.point = p4;

  planes.push_back(xwall_backl);
  planes.push_back(xwall_backr);

  // gamma veto part 1
  const double zcalo_bd_l[2] = {
    gveto_locator.get_wall_window_z(side, snemo::geometry::gveto_locator::WALL_BOTTOM),
    gveto_locator.get_wall_window_z(side, snemo::geometry::gveto_locator::WALL_TOP)};
  ROOT::Math::XYZPoint p5(0.0,0.0,zcalo_bd_l[0]); // negative z
  ROOT::Math::XYZVector norm5(0.0,0.0,1.0); // looking in pos z direction
  ROOT::Math::XYZPoint p6(0.0,0.0,zcalo_bd_l[1]); // positive z
  ROOT::Math::XYZVector norm6(0.0,0.0,-1.0); // looking in neg z direction
  Plane gv_backb;
  gv_backb.planeid = ++whichcalo;
  gv_backb.side = side;
  gv_backb.normal = norm5;
  gv_backb.point = p5;
  Plane gv_backt;
  gv_backt.planeid = ++whichcalo;
  gv_backt.side = side;
  gv_backt.normal = norm6;
  gv_backt.point = p6;

  planes.push_back(gv_backb);
  planes.push_back(gv_backt);

  // xwalls part 2
  side  = 1; // front
  const double ycalo_bd_r[2] = {
    xcalo_locator.get_wall_window_y(side, snemo::geometry::xcalo_locator::WALL_LEFT),
    xcalo_locator.get_wall_window_y(side, snemo::geometry::xcalo_locator::WALL_RIGHT)};
  ROOT::Math::XYZPoint p7(0.0,ycalo_bd_r[0],0.0); // negative y
  ROOT::Math::XYZVector norm7(0.0,1.0,0.0); // looking in pos y direction
  ROOT::Math::XYZPoint p8(0.0,ycalo_bd_r[1],0.0); // positive y
  ROOT::Math::XYZVector norm8(0.0,-1.0,0.0); // looking in neg y direction
  Plane xwall_frontl;
  xwall_frontl.planeid = ++whichcalo;
  xwall_frontl.side = side;
  xwall_frontl.normal = norm7
  xwall_frontl.point = p7;
  Plane xwall_frontr;
  xwall_frontr.planeid = ++whichcalo;
  xwall_frontr.side = side;
  xwall_frontr.normal = norm8;
  xwall_frontr.point = p8;

  planes.push_back(xwall_frontl);
  planes.push_back(xwall_frontr);

  // gamma veto part 2
  const double zcalo_bd_r[2] = {
    gveto_locator.get_wall_window_z(side, snemo::geometry::gveto_locator::WALL_BOTTOM),
    gveto_locator.get_wall_window_z(side, snemo::geometry::gveto_locator::WALL_TOP)};
  ROOT::Math::XYZPoint p9(0.0,0.0,zcalo_bd_r[0]); // negative z
  ROOT::Math::XYZVector norm9(0.0,0.0,1.0); // looking in pos z direction
  ROOT::Math::XYZPoint p10(0.0,0.0,zcalo_bd_r[1]); // positive z
  ROOT::Math::XYZVector norm10(0.0,0.0,-1.0); // looking in neg z direction
  Plane gv_frontb;
  gv_frontb.planeid = ++whichcalo;
  gv_frontb.side = side;
  gv_frontb.normal = norm9;
  gv_frontb.point = p9;
  Plane gv_frontt;
  gv_frontt.planeid = ++whichcalo;
  gv_frontt.side = side;
  gv_frontt.normal = norm10;
  gv_frontt.point = p10;
  
  planes.push_back(gv_frontb);
  planes.push_back(gv_frontt);
  
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
  for (auto entry : cls_defaults) { // get sdm::tracker_cluster handle
    // check all clusters for wire_candidate
    const sdm::calibrated_tracker_hit::collection_type & gg_hits_col = entry.get().get_hits();

    // into a function to check on wire_candidate with vertex pair<bool, bool> 
    // for <foil direction, calo direction>
    VertexInfo info = check_on_wire(gg_hits_col);
    info.clsid = entry.get().get_cluster_id();
    info.side = entry.get().get_hits().at(0).get().get_side();

    findmaxmin(info, gg_hits_col);
    all_info.push_back(info);
  }

  // ready to extrapolate, work with loop over all trj containers, check on wire_candidate via clsid
  VertexExtrapolator ve(planes); // has all the geometry information now

  const std::vector<std::vector<LineFit> >* ptr_lf = &(data_record__.get<std::vector<std::vector<LineFit> > >("mylinefits"));

  for (auto& onecluster: : *ptr_lf) {
    // loop over all types and get all available intersections
    for (LineFit lf : onecluster) {
      ve.setTrajectory(lf, all_info);
      std::pair<VertexInfo, std:pair<Ellipse, Ellipse> >  vt = ve.fullvertex();
    }
  }
  
  eventCounter++;
  return dpp::base_module::PROCESS_SUCCESS;
}



VertexInfo vertex_module::check_on_wire(const sdm::calibrated_tracker_hit::collection_type & data) {
  // use cluster member hits to check on presence of hits on the tracker perimeter.
  // If not then on either the calo (main and xwall) side or the foil side or both
  // there is a wire vertex candidate. Checking on gamma veto needs the trajectories
  // later, potentially turning a candidate into a proper wire vertex.

  VertexInfo vi;
  bool foilside = false;
  bool caloside = false;
  int foil_perimeter = 0; // layer 0
  int calo_perimeter = 8; // layer 8
  int xwall1 = 0; // row 0
  int xwall2 = 112; // row 112
  for (auto hit_handle : data) { // any cluster hit flips the boolean to true
    const sdm::calibrated_tracker_hit & hit = data.get();
    if (hit.get_layer() == foil_perimeter)
      foilside = true;
    if (hit.get_layer() == calo_perimeter)
      caloside = true;
    else if (hit.get_row() == xwall1)
      caloside = true;
    else if (hit.get_row() == xwall2)
      caloside = true;
  }

  vi.wire_candidate = !(foilside && caloside);
  vi.foilcalo = std::make_pair(foilside, caloside);

  return vi;
}

void vertex_module::findmaxmin(VertexInfo& vi, const sdm::calibrated_tracker_hit::collection_type & data)
{
  // Not knowing what the trajectory is, collect some geiger cell
  // information required for helices to determine valid 
  // plane intersection points and charge from curvature if possible
  MetaInfo mi;
  std::vector<int> dummy
  for (auto hit_handle : data) { 
    const sdm::calibrated_tracker_hit & hit = data.get();
    dummy.push_back(hit.get_layer()); // store x layer [0-8]
  }
  std::vector<int>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  std::vector<int>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  int min = *minit; // search target value
  int max = *maxit; // search target value
  int mymincount = (int) std::count(dummy.begin(), dummy.end(), min);
  int mymaxcount = (int) std::count(dummy.begin(), dummy.end(), max);
  int pos;
  for (int i=0;i<mymincount;i++) { // finds all min layer value entries
    pos = minit - dummy.begin();
    mi.hitid  = data.at(pos).get().get_id();
    mi.side   = data.at(pos).get().get_side();
    mi.row    = data.at(pos).get().get_row();
    mi.column = data.at(pos).get().get_layer();
    mi.wirex  = data.at(pos).get().get_x();
    mi.wirey  = data.at(pos).get().get_y();
    mi.zcoord = data.at(pos).get().get_z();
    vi.minx.push_back(mi);
    ++minit;
    it = std::find(minit, dummy.end(), min);
  }
  for (int i=0;i<mymaxcount;i++) { // finds all max layer value entries
    pos = maxit - dummy.begin();
    mi.hitid  = data.at(pos).get().get_id();
    mi.side   = data.at(pos).get().get_side();
    mi.row    = data.at(pos).get().get_row();
    mi.column = data.at(pos).get().get_layer();
    mi.wirex  = data.at(pos).get().get_x();
    mi.wirey  = data.at(pos).get().get_y();
    mi.zcoord = data.at(pos).get().get_z();
    vi.maxx.push_back(mi);
    ++maxit;
    it = std::find(maxit, dummy.end(), max);
  }

}
