/* Description:
 *
 *   Module for consolidating fitted clusters
 *
 * History:
 *
 */

#ifndef FALAISE_VERTEX_MODULE_H
#define FALAISE_VERTEX_MODULE_H 1

// Third party:
#include <string>

// - Bayeux/dpp:
#include <bayeux/dpp/base_module.h>
#include "bayeux/datatools/service_manager.h"
#include "bayeux/geomtools/manager.h"
#include "bayeux/geomtools/geometry_service.h"

// - Falaise:
#include <falaise/snemo/datamodels/calibrated_tracker_hit.h>
#include <falaise/snemo/geometry/locator_plugin.h>

// This project:
#include <vertex_library.h>

/// \brief Tracker consolidation module for fitted clusters
class vertex_module : public dpp::base_module
{
	
public:

	/// Constructor
	vertex_module(datatools::logger::priority = datatools::logger::PRIO_FATAL);
	
	/// Destructor
	virtual ~vertex_module();
	
	/// Initialization
	virtual void initialize(const datatools::properties  & setup_,
				datatools::service_manager   & service_manager_,
				dpp::module_handle_dict_type & module_dict_);
	
	/// Reset
	virtual void reset();
	
	/// Data record processing
	virtual process_status process(datatools::things & data_);
	

protected:
	VertexInfo check_on_wire(const snemo::datamodel::calibrated_tracker_hit::collection_type & data);

private:
	int eventCounter;
	std::string _TTD_label_; // not used
	std::string _TCD_label_;
	
	// geometry service
	const geomtools::manager* geometry_manager_; //!< The geometry manager
	const snemo::geometry::locator_plugin* locator_plugin_; //!< The SuperNEMO locator plugin

	// Macro to automate the registration of the module :
	DPP_MODULE_REGISTRATION_INTERFACE(vertex_module)
};

#endif
