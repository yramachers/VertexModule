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
  allPlanes.clear();
}

// main constructor
VertexExtrapolator::VertexExtrapolator(std::vector<Plane> pl)
{
  allInfo.clear();
  allPlanes = pl;
}


void VertexExtrapolator::intersect(int which) // main working method
{
  int xlf = lf.clid;
  int xhf = hf.clid;
  int xblf = blf.clid;
  switch(which) {
  case 0: // line case
    if (std::count_if(allInfo.begin(), allInfo.end(), [xlf](const VertexInfo& vi){return vi.clsid == xlf;}) > 0) {
      info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [xlf](const VertexInfo& vi){return vi.clsid == xlf;}) - allInfo.begin()); // linefit case
      calovertex.clear();
      intersect_line(); // sets all results
    }
    break;

  case 1: // helix case
    if (std::count_if(allInfo.begin(), allInfo.end(), [xhf](const VertexInfo& vi){return vi.clsid == xhf;}) > 0) {
      info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [xhf](const VertexInfo& vi){return vi.clsid == xhf;}) - allInfo.begin()); // helixfit case
      calovertex.clear();
      intersect_helix(); // sets all results
    }
    break;

  case 2: // broken line case
    if (std::count_if(allInfo.begin(), allInfo.end(), [xblf](const VertexInfo& vi){return vi.clsid == xblf;}) > 0) {
      info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [xblf](const VertexInfo& vi){return vi.clsid == xblf;}) - allInfo.begin()); // brokenlinefit case
      calovertex.clear();
      intersect_brokenline(); // sets all results
    }
    break;
  }
  return;
}

// Line section
// ************
void VertexExtrapolator::intersect_line()
{
  int side = info.side;
  std::vector<Line3d> lc = linecollection(side); // lf is known

  // check final calo: gvet; finalizes info.foilcalo.second to true or false
  zcheck(lc, side);
  //  std::cout << "after zcheck, foilcalo is (" << info.foilcalo.first << ", " << info.foilcalo.second << ")" << std::endl;

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
    ROOT::Math::XYZPoint pplusy (0.0, lf.ixy+lf.errixy, lf.ixz); // intercepts are
    ROOT::Math::XYZPoint pminusy(0.0, lf.ixy-lf.errixy, lf.ixz); // the intersection
    ROOT::Math::XYZPoint pplusz (0.0, lf.ixy, lf.ixz+lf.errixz); // points by definition
    ROOT::Math::XYZPoint pminusz(0.0, lf.ixy, lf.ixz-lf.errixz); // with foil at x=0 origin.
    
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
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
  }
  if ((!info.foilcalo.first) && (!info.foilcalo.second)) {
    // both vertices on wires, empty vertex rectangles
    foilvertex.axis1.clear(); // zero axes
    foilvertex.axis2.clear();
    foilvertex.planeid= -1;
    foilvertex.side   = -1;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
    calovertex.clear(); // empty vector
  }
  return; // end here
}


void VertexExtrapolator::set_calospot(std::vector<Line3d>& lc, Plane p, int side)
{
  // this sets rectangle axes for calovertex, the calorimeter rectangle container
  // depending on which unique plane is intersected by the best fit line

  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  // centre point from best fit line
  ROOT::Math::XYZVector ubest(1.0, lf.slxy, lf.slxz);
  ROOT::Math::XYZPoint pbest(0.0, lf.ixy, lf.ixz);
  Line3d best;
  best.u = ubest.Unit();
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
  //  std::cout << "wrong plane side, returns." << std::endl;
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_line_plane(Line3d& l, Plane p)
{
  ROOT::Math::XYZPoint i;
  double denom = p.normal.Dot(l.u);
  // std::cout << "Got p normal: (" << p.normal.x() << ", " << p.normal.y() << ", " << p.normal.z() << ")" << std::endl;
  // std::cout << "Got line u: (" << l.u.x() << ", " << l.u.y() << ", " << l.u.z() << ")" << std::endl;
  // std::cout << "denominator = " << denom << std::endl;
  if (fabs(denom)>0.0) { // not parallel
    double u = ((p.point - l.point).Dot(p.normal)) / denom;
    i.SetXYZ(l.point.x() + u * l.u.x(),
	     l.point.y() + u * l.u.y(),
	     l.point.z() + u * l.u.z());
    //    std::cout << "intersection point: (" << i.x() << ", " << i.y() << ", " << i.z() << ")" << std::endl;
    return i; // just the one intersection point
  }
  else { 
    i.SetXYZ(1001.0,10001.0,10001.0); // no intersection signature
    return i;
  }
}

std::vector<Line3d> VertexExtrapolator::linecollection(int side)
{
  ROOT::Math::XYZVector uplusy (1.0, lf.slxy+lf.errslxy, lf.slxz);
  ROOT::Math::XYZVector uminusy(1.0, lf.slxy-lf.errslxy, lf.slxz);
  ROOT::Math::XYZVector uplusz (1.0, lf.slxy, lf.slxz+lf.errslxz);
  ROOT::Math::XYZVector uminusz(1.0, lf.slxy, lf.slxz-lf.errslxz);
  ROOT::Math::XYZPoint pplusy (0.0, lf.ixy+lf.errixy, lf.ixz);
  ROOT::Math::XYZPoint pminusy(0.0, lf.ixy-lf.errixy, lf.ixz);
  ROOT::Math::XYZPoint pplusz (0.0, lf.ixy, lf.ixz+lf.errixz);
  ROOT::Math::XYZPoint pminusz(0.0, lf.ixy, lf.ixz-lf.errixz);
  std::vector<Line3d> collection;
  
  Line3d current;
  if (side ==0) { // Side = 0, negative x
    //sweep in y
    current.u     = uplusy.Unit();
    current.point = pplusy;
    collection.push_back(current);
    current.u     = uminusy.Unit();
    current.point = pminusy;
    collection.push_back(current);
    
    //sweep in z
    current.u     = uplusz.Unit();
    current.point = pplusz;
    collection.push_back(current);
    current.u     = uminusz.Unit();
    current.point = pminusz;
    collection.push_back(current);
    
    return collection;
  }
  else { // Side = 1, positive x
    //sweep in y
    current.u     = uplusy.Unit();
    current.point = pminusy;
    collection.push_back(current);
    current.u     = uminusy.Unit();
    current.point = pplusy;
    collection.push_back(current);
    
    //sweep in z
    current.u     = uplusz.Unit();
    current.point = pminusz;
    collection.push_back(current);
    current.u     = uminusz.Unit();
    current.point = pplusz;
    collection.push_back(current);
    
    return collection;
  }
}


// Broken Line section
// *******************
void VertexExtrapolator::intersect_brokenline()
{
  if (blf.breakpoints.empty()) { // unbroken line
    lf = blf.linefit1;
    intersect_line();
  }
  else {
    // to foil
    lf = blf.linefit1;
    bool previous_info = info.foilcalo.second;
    info.foilcalo.second = false; // enforce no extrapolation to calo for linefit1
    intersect_line();
    info.foilcalo.second = previous_info;
    // to calo
    lf = blf.linefit2;
    previous_info = info.foilcalo.first;
    info.foilcalo.first = false; // enforce no extrapolation to foil for linefit2
    intersect_line();
    info.foilcalo.first = previous_info;
  }
}


// Helix section
// ************
void VertexExtrapolator::intersect_helix()
{
  int side = info.side;
  std::vector<ROOT::Math::XYZPoint> isecs;
  std::vector<int> idcollection;
  Helix3d current;
  ROOT::Math::XYZPoint pbest(hf.xc, hf.yc, hf.zc); // centre
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
    else
      current.charge = 0; // undecided
  }
  if (side == 1) { // front tracker x>0
    if (pbest.y() > lower_bd)
      current.charge = -1;
    else if (pbest.y() < upper_bd)
      current.charge = 1;
    else
      current.charge = 0; // undecided
  }

  // piercing planes from here
  std::vector<Helix3d> hc = helixcollection(current.charge);

  // check final calo: gvet; finalizes info.foilcalo.second to true or false
  zcheck_helix(hc, side);

  if (info.foilcalo.second) { // only if a calo vertex is required
    for (Plane p : allPlanes) {
      // calo spots - determined by centre fit helix intersecting exactly one calo wall
      // and left-overs on other planes, if any
      set_calospot_helix(hc, p, side); 
    }
  }
  if (info.foilcalo.first) { // only if a foil vertex is required
    foilvertex.planeid = 10; // foil id
    foilvertex.side = 1; // irrelevant here
    Interval ybounds(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
    Interval zbounds(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

    ROOT::Math::XYZPoint origin(0.0, 0.0, 0.0);
    ROOT::Math::XYZVector xaxis(1.0, 0.0, 0.0);
    Plane foil;
    foil.planeid = 10;
    foil.side = 1; // irrelevant here
    foil.point = origin;
    foil.normal = xaxis;
    ROOT::Math::XYZPoint isec1 = intersect_helix_plane(hc.at(0), foil);
    ROOT::Math::XYZPoint isec2 = intersect_helix_plane(hc.at(1), foil);
    ROOT::Math::XYZPoint isec3 = intersect_helix_plane(hc.at(2), foil);
    ROOT::Math::XYZPoint isec4 = intersect_helix_plane(hc.at(3), foil);

    // axis1 along y
    if (point_plane_check_x(isec1, 1)) // side irrelevant here
      foilvertex.axis1.setbound(isec1.y());
    else 
      foilvertex.axis1.setbound(ybounds.to());
    if (point_plane_check_x(isec2, 1))
      foilvertex.axis1.setbound(isec2.y());
    else 
      foilvertex.axis1.setbound(ybounds.from());

    // axis2 along z
    if (point_plane_check_x(isec3, 1))
      foilvertex.axis2.setbound(isec3.z());
    else 
      foilvertex.axis2.setbound(zbounds.to());      
    if (point_plane_check_x(isec4, 1))
      foilvertex.axis2.setbound(isec4.z());
    else 
      foilvertex.axis2.setbound(zbounds.from());

    double original_area = fabs(isec2.y() - isec1.y()) * fabs(isec4.z() - isec3.z());
    foilvertex.areafraction = (foilvertex.axis1.width() * foilvertex.axis2.width()) / original_area;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
  }
  if ((!info.foilcalo.first) && (!info.foilcalo.second)) {
    // both vertices on wires, empty vertex rectangles
    foilvertex.axis1.clear(); // zero axes
    foilvertex.axis2.clear();
    foilvertex.planeid= -1;
    foilvertex.side   = -1;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
    calovertex.clear(); // empty vector
  }
  return;
}


std::vector<Helix3d> VertexExtrapolator::helixcollection(int charge)
{
  // vary centre y, pitch to get four intersection points around the centre tracing a rectangle like for lines.
  std::vector<Helix3d> collection;
  
  Helix3d current;
  current.radius = hf.radius; // changes with y-centre
  current.pitch  = hf.pitch; 
  current.charge = charge; // all the same charge

  ROOT::Math::XYZPoint pplusy(hf.xc, hf.yc + hf.erryc, hf.zc); // centre
  ROOT::Math::XYZPoint pminusy(hf.xc, hf.yc - hf.erryc, hf.zc); // centre
  ROOT::Math::XYZPoint pbest(hf.xc, hf.yc, hf.zc); // centre
  double pitchup = hf.pitch + hf.errpitch;
  double pitchdown = hf.pitch - hf.errpitch;


  //sweep in y
  current.point = pplusy;
  collection.push_back(current);
  current.point = pminusy;
  collection.push_back(current);
  
  //sweep in z
  current.point = pbest;
  current.pitch  = pitchup;
  collection.push_back(current);
  current.pitch = pitchdown;
  collection.push_back(current);
  
  return collection;
}


void VertexExtrapolator::set_calospot_helix(std::vector<Helix3d>& hc, Plane p, int side)
{
  // this sets rectangle axes for calovertex, the calorimeter rectangle container
  // depending on which unique plane is intersected by the best fit helix

  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  // centre point from best fit helix
  ROOT::Math::XYZPoint pbest(hf.xc, hf.yc, hf.zc);
  Helix3d best;
  best.point = pbest;
  best.pitch = hf.pitch;
  best.radius = hf.radius;
  best.charge = hc.front().charge; // all the same charge

  if (side == p.side) { // only planes on the correct side
    ROOT::Math::XYZPoint centre = intersect_helix_plane(best, p);
    if (centre.x()>10000.0 || centre.y()>10000.0 || centre.z()>10000.0) return; // no intersection signature, end here

    if (p.planeid<2) { // main wall
      if (point_plane_check_x(centre, p.side)) { // best fit hits this plane
	double area = mainwall_check_helix(hc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // xwall
	    double d = xwall_check_helix(hc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    double d = gveto_check_helix(hc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since xwall checked first
	    double d = gveto_check_helix(hc, allPlanes.at(indx2), area);
	} 
      }
    }
    else if ((p.planeid>=2 && p.planeid<4) || (p.planeid>=6 && p.planeid<8)) { // xwall
      if (point_plane_check_y(centre, side)) { // best fit hits this plane
	double area = xwall_check_helix(hc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    double d = mainwall_check_helix(hc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    double d = gveto_check_helix(hc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since mainwall checked first
	    double d = gveto_check_helix(hc, allPlanes.at(indx2), area);
	}
      }
    }
    else { // gveto
      if (point_plane_check_z(centre, side)) { // best fit hits this plane
	double area = gveto_check_helix(hc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    double d = mainwall_check_helix(hc, allPlanes.at(indx1), area);
	  else if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // should be set to xwall then
	    double d = xwall_check_helix(hc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be xwall since mainwall checked first
	    double d = xwall_check_helix(hc, allPlanes.at(indx2), area);
	}
      }
    }
  } // correct side, otherwise do nothing
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_helix_plane(Helix3d& h, Plane p)
{
  ROOT::Math::XYZPoint isec;
  if (p.planeid<2) // main wall
    isec = intersect_helix_mainw(h, p); 

  else if (p.planeid>=2 && p.planeid<4 || p.planeid>=6 && p.planeid<8) // xwall
    isec = intersect_helix_xwall(h, p); // +- y from p.planeid

  else if (p.planeid == 10) // foil
    isec = intersect_helix_foil(h, p); 

  else // gveto
    isec = intersect_helix_gveto(h, p); // +- z from p.planeid

  return isec;
}


double VertexExtrapolator::Pointdistance(ROOT::Math::XYZPoint p1, ROOT::Math::XYZPoint p2) 
{
  return TMath::Sqrt((p2-p1).Mag2()); // distance
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_helix_mainw(Helix3d& h, Plane p)
{
  // *** sort two points first, only one should come back ***
  std::vector<ROOT::Math::XYZPoint> pointcollection;
  ROOT::Math::XYZPoint imin;
  // parameter
  double x0 = p.point.x(); // +-x coordinate of main wall
  double xc = h.point.x(); // helix centre x
  double r  = h.radius;
  double arg = (x0 - xc) / r;
  if (fabs(arg)>=1)  // not hitting main wall
    return imin.SetXYZ(10001.0, 0.0, 0.0); // empty

  double pitch = h.pitch;
  double t  =  TMath::ACos(arg);

  double y  = h.point.y() + TMath::Sqrt(r*r - (x0-xc)*(x0-xc));
  double y2 = h.point.y() - TMath::Sqrt(r*r - (x0-xc)*(x0-xc));
  double z  = h.point.z() + pitch * t  / (2.0 * TMath::Pi());
  double z2 = h.point.z() - pitch * t  / (2.0 * TMath::Pi());
  ROOT::Math::XYZPoint i(x0,y,z);
  pointcollection.push_back(i);
  ROOT::Math::XYZPoint i2(x0,y2,z);  // four piercing points
  pointcollection.push_back(i2);
  ROOT::Math::XYZPoint i3(x0,y,z2);
  pointcollection.push_back(i3);
  ROOT::Math::XYZPoint i4(x0,y2,z2);
  pointcollection.push_back(i4);
  // reduce to one

  // pick the closest to the collection of calo-side tracker hits
  double minDistance = 10002.0; // dummy big
  std::vector<double> distances;
  for (auto ipoint : pointcollection) {
    for (auto mi : info.maxx) { // has to contain layer 8 by definition
      ROOT::Math::XYZPoint wirep(mi.wirex, mi.wirey, mi.zcoord);
      distances.push_back(Pointdistance(ipoint, wirep));
    }
    std::vector<double>::iterator mit = std::min_element(distances.begin(), distances.end());
    if (*mit < minDistance) {
      minDistance = *mit;
      imin = ipoint;
    }
    distances.clear();
  }

  return imin;
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_helix_foil(Helix3d& h, Plane p)
{
  // *** sort two points first, only one should come back ***
  std::vector<ROOT::Math::XYZPoint> pointcollection;
  ROOT::Math::XYZPoint imin;
  // parameter
  double x0 = p.point.x(); // x coordinate foil
  double xc = h.point.x(); // helix centre x
  double r  = h.radius;
  double arg = (x0 - xc) / r;
  if (fabs(arg)>=1)  // not hitting foil
    return imin.SetXYZ(10001.0, 0.0, 0.0); // empty

  double pitch = h.pitch;
  double t  =  TMath::ACos(arg);

  double y  = h.point.y() + TMath::Sqrt(r*r - (x0-xc)*(x0-xc));
  double y2 = h.point.y() - TMath::Sqrt(r*r - (x0-xc)*(x0-xc));
  double z  = h.point.z() + pitch * t  / (2.0 * TMath::Pi());
  double z2 = h.point.z() - pitch * t  / (2.0 * TMath::Pi());
  ROOT::Math::XYZPoint i(x0,y,z);
  pointcollection.push_back(i);
  ROOT::Math::XYZPoint i2(x0,y2,z);  // four piercing points
  pointcollection.push_back(i2);
  ROOT::Math::XYZPoint i3(x0,y,z2);
  pointcollection.push_back(i3);
  ROOT::Math::XYZPoint i4(x0,y2,z2);
  pointcollection.push_back(i4);
  // reduce to one

  // pick the closest to the collection of calo-side tracker hits
  double minDistance = 10002.0; // dummy big
  std::vector<double> distances;
  for (auto ipoint : pointcollection) {
    for (auto mi : info.minx) { // has to contain layer 0 by definition
      ROOT::Math::XYZPoint wirep(mi.wirex, mi.wirey, mi.zcoord);
      distances.push_back(Pointdistance(ipoint, wirep));
    }
    std::vector<double>::iterator mit = std::min_element(distances.begin(), distances.end());
    if (*mit < minDistance) {
      minDistance = *mit;
      imin = ipoint;
    }
    distances.clear();
  }

  return imin;
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_helix_xwall(Helix3d& h, Plane p)
{
  // *** sort two points first, only one should come back ***
  std::vector<ROOT::Math::XYZPoint> pointcollection;
  ROOT::Math::XYZPoint imin;
  // parameter
  double y0 = p.point.y(); //  coordinate of xwall
  double yc = h.point.y(); // helix centre x
  double r  = h.radius;
  double pitch = h.pitch;
  double arg = (y0 - yc);
  if (arg / r >= 1)  // not hitting xwall
    return imin.SetXYZ(0.0, 10001.0, 0.0); // empty

  double t  =  TMath::ACos(1/r * TMath::Sqrt(r*r - arg*arg));
  double x  = h.point.x() + TMath::Sqrt(r*r - arg*arg);
  double x2 = h.point.x() - TMath::Sqrt(r*r - arg*arg);
  double z  = h.point.z() + pitch * t / (2.0 * TMath::Pi());
  double z2 = h.point.z() - pitch * t / (2.0 * TMath::Pi());
  ROOT::Math::XYZPoint i(x,y0,z); // four piercing points
  ROOT::Math::XYZPoint i2(x2,y0,z);
  ROOT::Math::XYZPoint i3(x,y0,z2);
  ROOT::Math::XYZPoint i4(x2,y0,z2);
  pointcollection.push_back(i);
  pointcollection.push_back(i2);
  pointcollection.push_back(i3);
  pointcollection.push_back(i4);
  // reduce to one
  // pick the closest to the collection of xwall-side tracker hits
  double minDistance = 10002.0; // dummy big
  std::vector<double> distances;
  for (auto ipoint : pointcollection) {
    // check on extreme row presence first - one must be present by definition
    if (std::count_if(info.maxy.begin(), info.maxy.end(), [](const MetaInfo& mi){return mi.row == 112;}) > 0) {
      for (auto mi : info.maxy) {
	ROOT::Math::XYZPoint wirep(mi.wirex, mi.wirey, mi.zcoord);
	distances.push_back(Pointdistance(ipoint, wirep));
      }
    }
    else {
      for (auto mi : info.miny) {
	ROOT::Math::XYZPoint wirep(mi.wirex, mi.wirey, mi.zcoord);
	distances.push_back(Pointdistance(ipoint, wirep));
      }
    }
    std::vector<double>::iterator mit = std::min_element(distances.begin(), distances.end());
    if (*mit < minDistance) {
      minDistance = *mit;
      imin = ipoint;
    }
    distances.clear();
  }
  
  return imin;
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_helix_gveto(Helix3d& h, Plane p)
{
  ROOT::Math::XYZPoint i; // just one point piercing z-plane
  // parameter
  double z0 = p.point.y(); // y coordinate of xwall
  double zc = h.point.y(); // helix centre x
  double r  = h.radius;
  double pitch = h.pitch;
  double t = (z0 - zc) / pitch;
  if (t>=1)  // prevent full or more circle spirals to hit z plane
    return i.SetXYZ(0.0, 0.0, 10001.0); // empty
  double x  = h.point.x() + r * TMath::Cos(2.0 * TMath::Pi() * t);
  double y  = h.point.y() + r * TMath::Sin(2.0 * TMath::Pi() * t);
  
  return i.SetXYZ(x, y, z0);
}


// Checking section
// ****************
// line checks
double VertexExtrapolator::mainwall_check(std::vector<Line3d>& lc, Plane p, double area)
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


double VertexExtrapolator::xwall_check(std::vector<Line3d>& lc, Plane p, double area)
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
  if (point_plane_check_y(isec1, p.side)) // -y error
    spot.axis1.setbound(isec1.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_y(isec2, p.side)) // +y error
    spot.axis1.setbound(isec2.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_y(isec3, p.side)) // -z error
    spot.axis2.setbound(isec3.z());
  else {
    (p.side==1) ? spot.axis2.setbound(zbound.from()) : spot.axis2.setbound(zbound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 4 : 5, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 4 : 5);
  }

  if (point_plane_check_y(isec4, p.side)) // +z error
    spot.axis2.setbound(isec4.z());
  else {
    (p.side==1) ? spot.axis2.setbound(zbound.to()) : spot.axis2.setbound(zbound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 5 : 4, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 5 : 4);
  }
  
  double original_area;
  (area<0.0) ? original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.z() - isec3.z()) : original_area = area;
  spot.areafraction = (spot.axis1.width() * spot.axis2.width()) / original_area;
  calovertex.push_back(spot);
  return original_area;
}


double VertexExtrapolator::gveto_check(std::vector<Line3d>& lc, Plane p, double area)
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
  if (point_plane_check_z(isec1, p.side)) // -y error
    spot.axis1.setbound(isec1.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_z(isec2, p.side)) // +y error
    spot.axis1.setbound(isec2.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_z(isec3, p.side)) // -z error
    spot.axis2.setbound(isec3.y());
  else {
    (p.side==1) ? spot.axis2.setbound(ybound.from()) : spot.axis2.setbound(ybound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 2 : 3, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 2 : 3);
  }

  if (point_plane_check_z(isec4, p.side)) // +z error
    spot.axis2.setbound(isec4.y());
  else {
    (p.side==1) ? spot.axis2.setbound(ybound.to()) : spot.axis2.setbound(ybound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 3 : 2, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 3 : 2);
  }
  
  double original_area;
  (area<0.0) ? original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.y() - isec3.y()) : original_area = area;
  spot.areafraction = (spot.axis1.width() * spot.axis2.width()) / original_area;
  calovertex.push_back(spot);
  return original_area;
}


void VertexExtrapolator::zcheck(std::vector<Line3d>& lc, int side)
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

// Helix checks
double VertexExtrapolator::mainwall_check_helix(std::vector<Helix3d>& hc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  ROOT::Math::XYZPoint isec1 = intersect_helix_plane(hc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_helix_plane(hc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_helix_plane(hc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_helix_plane(hc.at(3), p);
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


double VertexExtrapolator::xwall_check_helix(std::vector<Helix3d>& hc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  ROOT::Math::XYZPoint isec1 = intersect_helix_plane(hc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_helix_plane(hc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_helix_plane(hc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_helix_plane(hc.at(3), p);
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  
  // axis 1 of rectangle
  if (point_plane_check_y(isec1, p.side)) // -y error
    spot.axis1.setbound(isec1.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_y(isec2, p.side)) // +y error
    spot.axis1.setbound(isec2.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_y(isec3, p.side)) // -z error
    spot.axis2.setbound(isec3.z());
  else {
    (p.side==1) ? spot.axis2.setbound(zbound.from()) : spot.axis2.setbound(zbound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 4 : 5, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 4 : 5);
  }

  if (point_plane_check_y(isec4, p.side)) // +z error
    spot.axis2.setbound(isec4.z());
  else {
    (p.side==1) ? spot.axis2.setbound(zbound.to()) : spot.axis2.setbound(zbound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 5 : 4, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 5 : 4);
  }
  
  double original_area;
  (area<0.0) ? original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.z() - isec3.z()) : original_area = area;
  spot.areafraction = (spot.axis1.width() * spot.axis2.width()) / original_area;
  calovertex.push_back(spot);
  return original_area;
}


double VertexExtrapolator::gveto_check_helix(std::vector<Helix3d>& hc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  ROOT::Math::XYZPoint isec1 = intersect_helix_plane(hc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_helix_plane(hc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_helix_plane(hc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_helix_plane(hc.at(3), p);
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  
  // axis 1 of rectangle
  if (point_plane_check_z(isec1, p.side)) // -y error
    spot.axis1.setbound(isec1.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_z(isec2, p.side)) // +y error
    spot.axis1.setbound(isec2.x());
  else {
    (p.side==1) ? spot.axis1.setbound(xbound_front.to()) : spot.axis1.setbound(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_z(isec3, p.side)) // -z error
    spot.axis2.setbound(isec3.y());
  else {
    (p.side==1) ? spot.axis2.setbound(ybound.from()) : spot.axis2.setbound(ybound.to());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 2 : 3, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 2 : 3);
  }

  if (point_plane_check_z(isec4, p.side)) // +z error
    spot.axis2.setbound(isec4.y());
  else {
    (p.side==1) ? spot.axis2.setbound(ybound.to()) : spot.axis2.setbound(ybound.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 3 : 2, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 3 : 2);
  }
  
  double original_area;
  (area<0.0) ? original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.y() - isec3.y()) : original_area = area;
  spot.areafraction = (spot.axis1.width() * spot.axis2.width()) / original_area;
  calovertex.push_back(spot);
  return original_area;
}


void VertexExtrapolator::zcheck_helix(std::vector<Helix3d>& hc, int side)
{
  Plane pbot = allPlanes.at(4); // same on either side
  Plane ptop = allPlanes.at(5);

  for (Helix3d& current : hc) {
    ROOT::Math::XYZPoint intersection_top = intersect_helix_plane(current, ptop);
    ROOT::Math::XYZPoint intersection_bot = intersect_helix_plane(current, pbot);
    if (intersection_top.z() < 10001.0) { // intersection signature
      if (point_plane_check_z(intersection_top, side)) // point overlaps gveto
	info.foilcalo.second = true; // not on wire in z
    }
    else if (intersection_bot.z()<10001.0) { // intersection signature
      if (point_plane_check_z(intersection_bot, side)) // point overlaps gveto
	info.foilcalo.second = true; // not on wire in z
    }
  }
}

// Point checking
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

// for charge orientation
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

