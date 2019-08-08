// us
#include <vertex_library.h>

// ROOT includes
#include <TMath.h>


// std
#include <iostream>
#include <algorithm>


// *****************************
// ** VertexExtrapolator methods
// *****************************
VertexExtrapolator::VertexExtrapolator()
{
  allInfo.clear();
  allPlanes.clear()
}

// main constructor
VertexExtrapolator::VertexExtrapolator(std::vector<Plane> pl)
{
  allInfo.clear();
  allPlanes = pl;
}


void VertexExtrapolator::intersect(int which) // main working method
{
  switch(which) {
  case 0: // line case
    info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [](const VertexInfo& vi){return vi.clsid == lf.clid}) - allInfo.begin()); // linefit case
    calovertex.clear();
    intersect_line(); // sets all results
    break;

  case 1: // helix case
    info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [](const VertexInfo& vi){return vi.clsid == hf.clid}) - allInfo.begin()); // helixfit case
    calovertex.clear();
    intersect_helix(); // sets all results
    break;

  case 2: // broken line case
    info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [](const VertexInfo& vi){return vi.clsid == blf.clid}) - allInfo.begin()); // brokenlinefit case
    calovertex.clear();
    intersect_brokenline(); // sets all results
    break;
  }
  return;
}

// Line section
// ************
void VertexExtrapolator::intersect_line()
{
  int side = info.side;
  std::vector<ROOT::Math::XYZPoint> lc = linecollection(side); // lf is known

  // check final calo: gvet; finalizes info.foilcalo.second to true or false
  zcheck(lc, side);

  if (info.foilcalo.second) { // only if a calo vertex is required
    for (Plane p : allPlanes) {
      // calo spots - determined by centre fit line intersecting exactly one calo wall
      // and left-overs on other planes, if any
      set_calospot(lc, p, side); 
    }
  }
  if (info.foilcalo.first) { // only if a foil vertex is required
    foilvertex.planeid = 10; // foil id
    foilvertex.side = 1; // irrelevant here
    ROOT::Math::XYZVector pplusy (0.0, lf.ixy+lf.errixy, lf.ixz); // intercepts are
    ROOT::Math::XYZVector pminusy(0.0, lf.ixy-lf.errixy, lf.ixz); // the intersection
    ROOT::Math::XYZVector pplusz (0.0, lf.ixy, lf.ixz+lf.errixz); // points by definition
    ROOT::Math::XYZVector pminusz(0.0, lf.ixy, lf.ixz-lf.errixz); // with foil at x=0 origin.
    
    Interval ybounds(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
    Interval zbounds(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
    // axis1 along y
    if (point_plane_check_x(pplusy, 1)) // side irrelevant here
      foilvertex.axis1.setbound(pplusy.y());
    else 
      foilvertex.axis1.setbound(ybounds.to());
    if (point_plane_check_x(pminusy, 1))
      foilvertex.axis1.setbound(pminusy.y());
    else 
      foilvertex.axis1.setbound(ybounds.from());
    // axis2 along z
    if (point_plane_check_x(pplusz, 1))
      foilvertex.axis2.setbound(pplusz.z());
    else 
      foilvertex.axis2.setbound(zbounds.to());      
    if (point_plane_check_x(pminusz, 1))
      foilvertex.axis2.setbound(pminusz.z());
    else 
      foilvertex.axis2.setbound(zbounds.from());

    double original_area = fabs(pplusy.y() - pminusy.y()) * fabs(pplusz.z() - pminusz.z());
    foilvertex.areafraction = (foilvertex.axis1.width() * foilvertex.axis2.width()) / original_area;
  }
  if ((!info.foilcalo.first) && (!info.foilcalo.second)) {
    // both vertices on wires, empty vertex rectangles
    foilvertex.axis1.clear(); // zero axes
    foilvertex.axis2.clear();
    foilvertex.planeid= -1;
    foilvertex.side   = -1;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
    calovertex.clear(); // empty vector
    return; // end here
  }
}


void VertexExtrapolator::set_calospot(std::vector<ROOT::Math::XYZPoint>& lc, Plane p, int side)
{
  // this sets rectangle axes for calovertex, the calorimeter rectangle container
  // depending on which unique plane is intersected by the best fit line

  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  // centre point from best fit line
  ROOT::Math::XYZVector ubest(1.0, lf.slxy, lf.slxz);
  ROOT::Math::XYZVector pbest(0.0, lf.ixy, lf.ixz);
  Line3D best;
  best.u = ubest;
  best.point = pbest;

  if (side == p.side) { // only planes on the correct side
    ROOT::Math::XYZPoint centre = intersect_line_plane(best, p);
    if (centre.x()>1000.0 &&centre.y()>10000.0 &&centre.z()>10000.0) return; // no intersection signature, end here

    if (p.planeid<2) { // main wall
      if (point_plane_check_x(centre, p.side)) { // best fit hits this plane
	double area = mainwall_check(lc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // xwall
	    double d = xwall_check(lc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    double d = gveto_check(lc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since xwall checked first
	    double d = gveto_check(lc, allPlanes.at(indx2), area);
	} 
      }
    }
    else if ((p.planeid>=2 && p.planeid<4) || (p.planeid>=6 && p.planeid<8)) { // xwall
      if (point_plane_check_y(centre, side)) { // best fit hits this plane
	double area = xwall_check(lc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    double d = mainwall_check(lc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    double d = gveto_check(lc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since mainwall checked first
	    double d = gveto_check(lc, allPlanes.at(indx2), area);
	}
      }
    }
    else { // gveto
      if (point_plane_check_z(centre, side)) { // best fit hits this plane
	double area = gveto_check(lc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    double d = mainwall_check(lc, allPlanes.at(indx1), area);
	  else if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // should be set to xwall then
	    double d = xwall_check(lc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be xwall since mainwall checked first
	    double d = xwall_check(lc, allPlanes.at(indx2), area);
	}
      }
    }
  } // correct side, otherwise do nothing
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_line_plane(Line3d& l, Plane p)
{
  double denom = p.normal.Dot(l.u);
  if (denom>0.0) { // not parallel
    double u = (p.normal.Dot(p.point) - p.normal.Dot(l.point)) / denom;
    ROOT::Math::XYZPoint i(l.point.x() + u * l.u.x(),
			   l.point.y() + u * l.u.y(),
			   l.point.z() + u * l.u.z());
    return i; // just the one intersection point
  }
  else { 
    ROOT::Math::XYZPoint i(1001.0,10001.0,10001.0); // no intersection signature
    return i;
  }
}

std::vector<ROOT::Math::XYZPoint> VertexExtrapolator::linecollection(int side)
{
  ROOT::Math::XYZVector uplusy (1.0, lf.slxy+lf.errslxy, lf.slxz);
  ROOT::Math::XYZVector uminusy(1.0, lf.slxy-lf.errslxy, lf.slxz);
  ROOT::Math::XYZVector uplusz (1.0, lf.slxy, lf.slxz+lf.errslxz);
  ROOT::Math::XYZVector uminusz(1.0, lf.slxy, lf.slxz-lf.errslxz);
  ROOT::Math::XYZVector pplusy (0.0, lf.ixy+lf.errixy, lf.ixz);
  ROOT::Math::XYZVector pminusy(0.0, lf.ixy-lf.errixy, lf.ixz);
  ROOT::Math::XYZVector pplusz (0.0, lf.ixy, lf.ixz+lf.errixz);
  ROOT::Math::XYZVector pminusz(0.0, lf.ixy, lf.ixz-lf.errixz);
  std::vector<ROOT::Math::XYZPoint> collection;
  
  Line3D current;
  if (side ==0) { // Side = 0, negative x
    //sweep in y
    current.u     = uplusy;
    current.point = pplusy;
    collecton.push_back(current);
    current.u     = uminusy;
    current.point = pminusy;
    collecton.push_back(current);
    
    //sweep in z
    current.u     = uplusz;
    current.point = pplusz;
    collecton.push_back(current);
    current.u     = uminusz;
    current.point = pminusz;
    collecton.push_back(current);
    
    return collection;
  }
  else { // Side = 1, positive x
    //sweep in y
    current.u     = uplusy;
    current.point = pminusy;
    collecton.push_back(current);
    current.u     = uminusy;
    current.point = pplusy;
    collecton.push_back(current);
    
    //sweep in z
    current.u     = uplusz;
    current.point = pminusz;
    collecton.push_back(current);
    current.u     = uminusz;
    current.point = pplusz;
    collecton.push_back(current);
    
    return collection;
  }
}


// Broken Line section
// *******************
void VertexExtrapolator::intersect_brokenline()
{
}


// Helix section
// ************
void VertexExtrapolator::intersect_helix()
{
  int side = info.side;
  std::vector<ROOT::Math::XYZPoint> isecs;
  std::vector<int> idcollection;
  Helix3d current;
  ROOT::Math::XYZVector pbest(hf.xc, hf.yc, hf.zc); // centre
  current.point  = pbest;
  current.radius = hf.radius;
  current.pitch  = hf.pitch;
  
  // helix charge from fitting the centre
  double lower_bd = findLowerYBound();
  double upper_bd = findUpperYBound();
  if (side == 0) { // back tracker x<0
    if (pbest.y() > lower_bd)
      current.charge = 1; // positron, right curvature
    else if (pbest.y() < upper_bd)
      current.charge = -1; // electron, left curvature
  }
  if (side == 1) { // front tracker x>0
    if (pbest.y() > lower_bd)
      current.charge = -1;
    else if (pbest.y() < upper_bd)
      current.charge = 1;
  }

  int which = 0;
  for (Plane p : allPlanes) {
    if (p.planeid<2)
      which = 0; // main wall
    else if (p.planeid>=2 && p.planeid<4 || p.planeid>=6 && p.planeid<8)
      which = 1; // xwall
    else
      which = 2; // gamma veto
    std::vector<ROOT::Math::XYZPoint> intersections = intersect_helix_plane(current, p, which);
    // check intersection
    if (point_plane_check(intersection)) {
      isecs.push_back(intersection);
      idcollection.push_back(p.planeid);
    }
  }
}


std::vector<ROOT::Math::XYZPoint> VertexExtrapolator::intersect_helix_plane(Helix3d h, Plane p, int which)
{
  ROOT::Math::XYZPoint isec;
  switch (which) {
  case 0: 
    isec = intersect_helix_mainw(h, p); 
    break;
  case 1:
    isec = intersect_helix_xwall(h, p); // +- y from p.planeid
    break;
  case 2:
    isec = intersect_helix_gveto(h, p); // +- z from p.planeid
    break;
  }
  return isec;
}


std::vector<ROOT::Math::XYZPoint> VertexExtrapolator::intersect_helix_mainw(Helix3d h, Plane p)
{
  std::vector<ROOT::Math::XYZPoint> pointcollection;
  // parameter
  double x0 = p.point.x(); // +-x coordinate of main wall
  double xc = h.point.x(); // helix centre x
  double r  = h.radius;
  double arg = (x0 - xc) / r;
  if (arg>=1)  // not permitted
    return pointcollection; // empty

  double pitch = h.pitch;
  double t  =  TMath::Acos(arg);

  double y  = h.point.y() + TMath::Sqrt(r*r - (x0-xc)*(x0-xc));
  double y2 = h.point.y() - TMath::Sqrt(r*r - (x0-xc)*(x0-xc));
  double z  = h.point.z() + pitch * t  / (2.0 * TMath::Pi());
  double z2 = h.point.z() - pitch * t  / (2.0 * TMath::Pi());
  ROOT::Math::XYZPoint i(x0,y,z);
  pointcollection.push_back(i);
  ROOT::Math::XYZPoint i2(x0,y2,z2);  // two closest plane piercing points
  pointcollection.push_back(i2);
  return pointcollection;
}


std::vector<ROOT::Math::XYZPoint> VertexExtrapolator::intersect_helix_xwall(Helix3d h, Plane p)
{
  std::vector<ROOT::Math::XYZPoint> pointcollection;
  // parameter
  double y0 = p.point.y(); //  coordinate of main wall
  double yc = h.point.y(); // helix centre x
  double r  = h.radius;
  double pitch = h.pitch;
  double arg = (y0 - yc);
  double t  =  TMath::Acos(1/r * TMath::Sqrt(r*r - arg*arg));

  double x  = h.point.x() + TMath::Sqrt(r*r - arg*arg);
  double x2 = h.point.x() - TMath::Sqrt(r*r - arg*arg);
  double z  = h.point.z() + pitch * t / (2.0 * TMath::Pi());
  double z2 = h.point.z() - pitch * t / (2.0 * TMath::Pi());
  ROOT::Math::XYZPoint i(x,y0,z); // two closest plane piercing points
  ROOT::Math::XYZPoint i2(x2,y0,z2);
  pointcollection.push_back(i);
  pointcollection.push_back(i2);
  return pointcollection;
}


std::vector<ROOT::Math::XYZPoint> VertexExtrapolator::intersect_helix_gveto(Helix3d h, Plane p)
{
  std::vector<ROOT::Math::XYZPoint> pointcollection;
  // parameter
  double z0 = p.point.y(); // y coordinate of xwall
  double zc = h.point.y(); // helix centre x
  double r  = h.radius;
  double pitch = h.pitch;
  double t = (z0 - zc) / pitch;
  if (t>=1)  // prevent full or more circle spirals to hit z plane
    return pointcollection; // empty
  double x  = h.point.x() + r * TMath::Cos(2.0 * TMath::Pi() * t);
  double y  = y.point.y() + r * TMath::Sin(2.0 * TMath::Pi() * t);
  ROOT::Math::XYZPoint i(x,y,z0); // just one point piercing z-plane
  pointcollection.push_back(i);
  return pointcollection;
}


// Checking section
// ****************
double VertexExtrapolator::mainwall_check(std::vector<ROOT::Math::XYZPoint>& lc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  ROOT::Math::XYZPoint isec1 = intersect_line_plane(lc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_line_plane(lc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_line_plane(lc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_line_plane(lc.at(3), p);
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  
  // axis 1 of rectangle
  if (point_plane_check_x(isec1, p.side)) // -y error
    spot.axis1.setbound(isec1.y());
  else {
    (p.side==1) ? spot.axis1.setbound(ybound.from()) : spot.axis1.setbound(ybound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 2 : 3, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 2 : 3);
  }
  if (point_plane_check_x(isec2, p.side)) // +y error
    spot.axis1.setbound(isec2.y());
  else {
    (p.side==1) ? spot.axis1.setbound(ybound.to()) : spot.axis1.setbound(ybound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 3 : 2, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 3 : 2);
  }  
  // axis 2 of rectangle
  if (point_plane_check_x(isec3, p.side)) // -z error
    spot.axis2.setbound(isec3.z());
  else {
    (p.side==1) ? spot.axis2.setbound(zbound.from()) : spot.axis2.setbound(zbound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 4 : 5, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 4 : 5);
  }
  if (point_plane_check_x(isec4, p.side)) // +z error
    spot.axis2.setbound(isec4.z());
  else {
    (p.side==1) ? spot.axis2.setbound(zbound.to()) : spot.axis2.setbound(zbound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 5 : 4, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 5 : 4);
  }  
  double original_area;
  double width1 = spot.axis1.width();
  double width2 = spot.axis2.width();
  (area<0.0) ? original_area = fabs(isec2.y() - isec1.y()) * fabs(isec4.z() - isec3.z()) : original_area = area;
  spot.areafraction = (width1 * width2) / original_area;
  calovertex.push_back(spot);
  return original_area;
}


double VertexExtrapolator::xwall_check(std::vector<ROOT::Math::XYZPoint>& lc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  ROOT::Math::XYZPoint isec1 = intersect_line_plane(lc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_line_plane(lc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_line_plane(lc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_line_plane(lc.at(3), p);
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  
  // axis 1 of rectangle
  if (point_plane_check_y(isec1, side)) // -y error
    spot.axis1.setbound(isec1.x());
  else {
    (side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_y(isec2, side)) // +y error
    spot.axis1.setbound(isec2.x());
  else {
    (side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_y(isec3, side)) // -z error
    spot.axis2.setbound(isec3.z());
  else {
    (side==1) ? spot.axis2.setbound(zbound.from()) : spot.axis2.setbound(zbound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 4 : 5, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 4 : 5);
  }

  if (point_plane_check_y(isec4, side)) // +z error
    spot.axis2.setbound(isec4.z());
  else {
    (side==1) ? spot.axis2.setbound(zbound.to()) : spot.axis2.setbound(zbound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 5 : 4, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 5 : 4);
  }
  
  double original_area;
  (area<0.0) ? original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.z() - isec3.z()) : original_area = area;
  spot.areafraction = (spot.axis1.width() * spot.axis2.width()) / original_area;
  calovertex.push_back(spot);
  return original_area;
}


double VertexExtrapolator::gveto_check(std::vector<ROOT::Math::XYZPoint>& lc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  ROOT::Math::XYZPoint isec1 = intersect_line_plane(lc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_line_plane(lc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_line_plane(lc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_line_plane(lc.at(3), p);
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  
  // axis 1 of rectangle
  if (point_plane_check_z(isec1, side)) // -y error
    spot.axis1.setbound(isec1.x());
  else {
    (side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_z(isec2, side)) // +y error
    spot.axis1.setbound(isec2.x());
  else {
    (side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_z(isec3, side)) // -z error
    spot.axis2.setbound(isec3.y());
  else {
    (side==1) ? spot.axis2.setbound(ybound.from()) : spot.axis2.setbound(ybound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 2 : 3, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 2 : 3);
  }

  if (point_plane_check_z(isec4, side)) // +z error
    spot.axis2.setbound(isec4.y());
  else {
    (side==1) ? spot.axis2.setbound(ybound.to()) : spot.axis2.setbound(ybound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 3 : 2, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 3 : 2);
  }
  
  double original_area;
  (area<0.0) ? original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.y() - isec3.y()) : original_area = area;
  spot.areafraction = (spot.axis1.width() * spot.axis2.width()) / original_area;
  calovertex.push_back(spot);
  return original_area;
}


void VertexExtrapolator::zcheck(std::vector<ROOT::Math::XYZPoint>& lc, int side)
{
  Plane pbot = allPlanes.at(4); // same on either side
  Plane ptop = allPlanes.at(5);

  for (Line3d& current : lc) {
    ROOT::Math::XYZPoint intersection_top = intersect_line_plane(current, ptop);
    ROOT::Math::XYZPoint intersection_bot = intersect_line_plane(current, pbot);
    if (intersection_top.x()>1000.0 &&intersection_top.y()>10000.0 &&intersection_top.z()>10000.0) continue; // no intersection signature

    if (point_plane_check_z(intersection_top, side)) // point overlaps gveto
      info.foilcalo.second = true; // not on wire in z
    if (point_plane_check_z(intersection_bot, side)) // point overlaps gveto
      info.foilcalo.second = true; // not on wire in z
  }
}


bool VertexExtrapolator::point_plane_check_x(ROOT::Math::XYZPoint point, int side)
{
  // within bounds of main wall and correct side.
  Interval ybounds(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbounds(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  bool check = false;
  if (side ==0 && ybounds.contains(point.y()) && zbounds.contains(point.z()))
    check = true;
  else if (side ==1 && ybounds.contains(point.y()) && zbounds.contains(point.z()))
    check = true;
  return check;
}

bool VertexExtrapolator::point_plane_check_y(ROOT::Math::XYZPoint point, int side)
{
  // within bounds of xwall and correct side.
  Interval xbounds_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbounds_front(0.0, allPlanes.at(1).point.x());
  Interval zbounds(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  bool check = false;
  if (side ==0 && xbounds_back.contains(point.x()) && zbounds.contains(point.z()))
    check = true;
  else if (side ==1 && xbounds_front.contains(point.x()) && zbounds.contains(point.z()))
    check = true;
  return check;
}

bool VertexExtrapolator::point_plane_check_z(ROOT::Math::XYZPoint point, int side)
{
  // within bounds of gveto and correct side.
  Interval xbounds_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbounds_front(0.0, allPlanes.at(1).point.x());
  Interval ybounds(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());

  bool check = false;
  if (side ==0 && xbounds_back.contains(point.x()) && ybounds.contains(point.y()))
    check = true;
  else if (side ==1 && xbounds_front.contains(point.x()) && ybounds.contains(point.y()))
    check = true;
  return check;
}


double VertexExtrapolator::findLowerYBound() {
  std::vector<double> dummy;
  for (auto val : info.minx)
    dummy.push_back(val.wirey);
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  double min1 = *minit;
  dummy.clear();
  for (auto val : info.maxx)
    dummy.push_back(val.wirey);
  minit = std::min_element(dummy.begin(), dummy.end());
  double min2 = *minit;
  return (min1<=min2) ? min1 : min2;
}

double VertexExtrapolator::findUpperYBound() {
  std::vector<double> dummy;
  for (auto val : info.minx)
    dummy.push_back(val.wirey);
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  double max1 = *maxit;
  dummy.clear();
  for (auto val : info.maxx)
    dummy.push_back(val.wirey);
  maxit = std::max_element(dummy.begin(), dummy.end());
  double max2 = *maxit;
  return (max1>=max2) ? max1 : max2;
}

// *****
// ** Utility Interval methods
// *****
Interval::Interval()
{
  lower = 0.0;
  upper = 0.0;
}


Interval::Interval(double s, double e)
{
  if (s <= e) { // start should be smaller than end
    lower = s;
    upper = e;
  }
  else { // swap order
    lower = e;
    upper = s;
  }
}

bool Interval::overlap(Interval other)
{
  if (lower < other.from())
    return upper > other.from();
  else if (upper > other.to())
    return lower < other.to();
  else return true; // this a subset
}

bool Interval::contains(double value)
{
  return (value <= upper && value>=lower);
}

void Interval::setbound(double value)
{
  // equals need no setting
  if (value < upper)
    lower = value;
  else if (value > lower)
    upper = value;
}

