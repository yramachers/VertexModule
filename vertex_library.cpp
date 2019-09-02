// us
#include <vertex_library.h>

// std
#include <iostream>
#include <algorithm>
#include <cmath>


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
  double sum_area = 0.0;
  int side = info.side;
  std::vector<Line3d> lc = linecollection(side); // lf is known

  // check final calo: gvet; finalizes info.foilcalo.second to true or false
  zcheck(lc, side);
  //  std::cout << "after zcheck, foilcalo is (" << info.foilcalo.first << ", " << info.foilcalo.second << ")" << std::endl;

  bool foilside = true;
  if (info.foilcalo.first && info.foilcalo.second) { // vertices on foil and calo
    // vertex on foil
    set_foilspot();

    // calo spots - determined by centre fit line intersecting exactly one calo wall
    // and left-overs on other planes, if any
    for (Plane p : allPlanes) { 
      double d = set_calospot(lc, p, side); 
      if (d>0.0) {
	sum_area += d;
	//	std::cout << "sum area = " << sum_area << std::endl;
      }
    }

    // correct the calo vertex entries to proper area fractions
    for (unsigned int i=0; i<calovertex.size(); i++) 
      calovertex.at(i).areafraction = (calovertex.at(i).axis1.width() * calovertex.at(i).axis2.width()) / sum_area;

    // no wire vertex
    wirevertex.axis1.clear(); // zero axes
    wirevertex.axis2.clear();
    wirevertex.axis3.clear();
    wirevertex.planeid= -1;
    wirevertex.side   = -1;
    wirevertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
  }
  else if (info.foilcalo.first && !info.foilcalo.second) { // only a foil vertex is required
    // vertex on foil
    set_foilspot();

    calovertex.clear(); // empty vector

    set_wirevertex(!foilside);
  }
  else if (!info.foilcalo.first && info.foilcalo.second) { // only a calo vertex is required
    foilvertex.axis1.clear(); // zero axes
    foilvertex.axis2.clear();
    foilvertex.axis3.clear();
    foilvertex.planeid= -1;
    foilvertex.side   = -1;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours

    set_wirevertex(foilside);

    for (Plane p : allPlanes) {
      double d = set_calospot(lc, p, side); 
      if (d>0.0) {
	sum_area += d;
	//	std::cout << "sum area = " << sum_area << std::endl;
      }
    }
    
    // correct the calo vertex entries to proper area fractions
    for (unsigned int i=0; i<calovertex.size(); i++)
      calovertex.at(i).areafraction = (calovertex.at(i).axis1.width() * calovertex.at(i).axis2.width()) / sum_area;
  }
  else {
    // both vertices on wires, empty vertex rectangles
    foilvertex.axis1.clear(); // zero axes
    foilvertex.axis2.clear();
    foilvertex.axis3.clear();
    foilvertex.planeid= -1;
    foilvertex.side   = -1;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours

    calovertex.clear(); // empty vector

    set_wirevertex(foilside);
    set_wirevertex(!foilside);
  }
  return; // end here
}



void VertexExtrapolator::set_foilspot()
{
    foilvertex.planeid = 10; // foil id
    foilvertex.side = 1; // irrelevant here
    ROOT::Math::XYZPoint pplusy (0.0, lf.ixy+lf.errixy, lf.ixz); // intercepts are
    ROOT::Math::XYZPoint pminusy(0.0, lf.ixy-lf.errixy, lf.ixz); // the intersection
    ROOT::Math::XYZPoint pplusz (0.0, lf.ixy, lf.ixz+lf.errixz); // points by definition
    ROOT::Math::XYZPoint pminusz(0.0, lf.ixy, lf.ixz-lf.errixz); // with foil at x=0 origin.
    
    Interval ybounds(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
    Interval zbounds(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
    Interval checky(lf.ixy-lf.errixy, lf.ixy+lf.errixy); // orders low and high values
    Interval checkz(lf.ixz-lf.errixz, lf.ixz+lf.errixz);
    if (ybounds.overlap(checky)) {
      // axis1 along y
      if (point_plane_check_x(pplusy, true)) 
	foilvertex.axis1.sethigh(checky.to());
      else 
	foilvertex.axis1.sethigh(ybounds.to());
      if (point_plane_check_x(pminusy, true))
	foilvertex.axis1.setlow(checky.from());
      else 
	foilvertex.axis1.setlow(ybounds.from());
    }
    else {
      foilvertex.axis1.clear();
    }
    if (zbounds.overlap(checkz)) {
      // axis2 along z
      if (point_plane_check_x(pplusz, false))
	foilvertex.axis2.sethigh(checkz.to());
      else 
	foilvertex.axis2.sethigh(zbounds.to());      
      if (point_plane_check_x(pminusz, false))
	foilvertex.axis2.setlow(checkz.from());
      else 
	foilvertex.axis2.setlow(zbounds.from());
    }
    else {
      foilvertex.axis2.clear();
    }
    double original_area = fabs(pplusy.y() - pminusy.y()) * fabs(pplusz.z() - pminusz.z());
    if (fabs(original_area)>0.0)
      foilvertex.areafraction = (foilvertex.axis1.width() * foilvertex.axis2.width()) / original_area;
    else
      foilvertex.areafraction = 0.0;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
    foilvertex.axis3.clear();
}


double VertexExtrapolator::set_calospot(std::vector<Line3d>& lc, Plane p, int side)
{
  double area;
  double sum_area = 0.0;
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
    std::cout << "hit side: " << p.side << " id = " << p.planeid << std::endl;
    if (centre.x()>1000.0 &&centre.y()>10000.0 &&centre.z()>10000.0) return 0.0; // no intersection signature, end here

    if (p.planeid<2) { // main wall
      //      std::cout << "main calo hit: (" << centre.x() << ", " << centre.y() << ", " << centre.z() << ")" << std::endl;
      if (point_plane_check_x(centre, true) && point_plane_check_x(centre, false)) { // best fit hits this plane in both axes
	area = mainwall_check(lc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // xwall
	    sum_area = xwall_check(lc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    sum_area = gveto_check(lc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since xwall checked first
	    sum_area = gveto_check(lc, allPlanes.at(indx2), sum_area);
	} 
      }
    }
    else if ((p.planeid>=2 && p.planeid<4) || (p.planeid>=6 && p.planeid<8)) { // xwall
      if (point_plane_check_y(centre, side, true) && point_plane_check_y(centre, side, false)) { // best fit hits this plane in both axes
	area = xwall_check(lc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    sum_area = mainwall_check(lc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    sum_area = gveto_check(lc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since mainwall checked first
	    sum_area = gveto_check(lc, allPlanes.at(indx2), sum_area);
	}
      }
    }
    else { // gveto
      if (point_plane_check_z(centre, side, true) && point_plane_check_z(centre, side, false)) { // best fit hits this plane in both axes
	area = gveto_check(lc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    sum_area = mainwall_check(lc, allPlanes.at(indx1), area);
	  else if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // should be set to xwall then
	    sum_area = xwall_check(lc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be xwall since mainwall checked first
	    sum_area = xwall_check(lc, allPlanes.at(indx2), sum_area);
	}
      }
    }

  } // correct side, otherwise do nothing
  //  std::cout << "wrong plane side, returns." << std::endl;
  return sum_area==0.0 ? area : sum_area;
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
  double sum_area = 0.0;
  int side = info.side;
  std::vector<ROOT::Math::XYZPoint> isecs;
  std::vector<int> idcollection;
  Helix3d current;
  ROOT::Math::XYZPoint pbest(hf.xc, hf.yc, hf.zc); // centre
  current.point  = pbest;
  current.radius = hf.radius;
  current.pitch  = hf.pitch;
  
  // helix charge from fitting the centre
  double lower_bd = findLowerYBoundatEnds();
  double upper_bd = findUpperYBoundatEnds();
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
  std::cout << "made charge: " << current.charge << std::endl;
  //  std::cout << "from pbest y: " << pbest.y() << " and lower_bd=" << lower_bd << " and upper_bd=" << upper_bd << std::endl;

  // piercing planes from here
  std::vector<Helix3d> hc = helixcollection(current.charge);

  // check final calo: gvet; finalizes info.foilcalo.second to true or false
  zcheck_helix(hc, side);

  bool foilside = true;
  if (info.foilcalo.first && info.foilcalo.second) {
    // vertex on foil
    set_foilspot_helix(hc);

    // calo spots - determined by centre fit helix intersecting exactly one calo wall
    // and left-overs on other planes, if any
    for (Plane p : allPlanes) {
      double d = set_calospot_helix(hc, p, side); 
      if (d>0.0) 
	sum_area += d;
    }
    // correct the calo vertex entries to proper area fractions
    for (unsigned int i=0; i<calovertex.size(); i++) 
      calovertex.at(i).areafraction = (calovertex.at(i).axis1.width() * calovertex.at(i).axis2.width()) / sum_area;

    // no wire vertex
    wirevertex.axis1.clear(); // zero axes
    wirevertex.axis2.clear();
    wirevertex.axis3.clear();
    wirevertex.planeid= -1;
    wirevertex.side   = -1;
    wirevertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
  }
  else if (info.foilcalo.first && !info.foilcalo.second) {
    // vertex on foil
    set_foilspot_helix(hc);

    calovertex.clear(); // empty vector

    set_wirevertex(!foilside);
  }
  else if (!info.foilcalo.first && info.foilcalo.second) {
    foilvertex.axis1.clear(); // zero axes
    foilvertex.axis2.clear();
    foilvertex.axis3.clear();
    foilvertex.planeid= -1;
    foilvertex.side   = -1;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours

    set_wirevertex(foilside);

    for (Plane p : allPlanes) {
      double d = set_calospot_helix(hc, p, side); 
      if (d>0.0)
	sum_area += d;
    }
    // correct the calo vertex entries to proper area fractions
    for (unsigned int i=0; i<calovertex.size(); i++)
      calovertex.at(i).areafraction = (calovertex.at(i).axis1.width() * calovertex.at(i).axis2.width()) / sum_area;
  }
  else {
    // both vertices on wires, empty vertex rectangles
    foilvertex.axis1.clear(); // zero axes
    foilvertex.axis2.clear();
    foilvertex.axis3.clear();
    foilvertex.planeid= -1;
    foilvertex.side   = -1;
    foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours

    calovertex.clear(); // empty vector

    set_wirevertex(foilside);
    set_wirevertex(!foilside);
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


void VertexExtrapolator::set_foilspot_helix(std::vector<Helix3d>& hc) {
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
  std::cout << "helix foil hit: (" << isec1.x() << ", " << isec1.y() << ", " << isec1.z() << ")" << std::endl;
  std::cout << "helix foil hit: (" << isec2.x() << ", " << isec2.y() << ", " << isec2.z() << ")" << std::endl;
  std::cout << "helix foil hit: (" << isec3.x() << ", " << isec3.y() << ", " << isec3.z() << ")" << std::endl;
  std::cout << "helix foil hit: (" << isec4.x() << ", " << isec4.y() << ", " << isec4.z() << ")" << std::endl;
  
  Interval checky(isec2.y(),isec1.y());
  Interval checkz(isec4.z(),isec3.z());
  if (checky.swapped()) { // swap points
    ROOT::Math::XYZPoint dummy = isec1;
    isec1 = isec2;
    isec2 = dummy;
  }
  if (checkz.swapped()) { // swap points
    ROOT::Math::XYZPoint dummy = isec3;
    isec3 = isec4;
    isec4 = dummy;
  }
  
  // assumes order to points isec1,2 and 3,4
  if (ybounds.overlap(checky)) {
    // axis1 along y
    if (point_plane_check_x(isec1, true)) 
      foilvertex.axis1.sethigh(checky.to());
    else 
      foilvertex.axis1.sethigh(ybounds.to());
    if (point_plane_check_x(isec2, true))
      foilvertex.axis1.setlow(checky.from());
    else 
      foilvertex.axis1.setlow(ybounds.from());
  }
  else {
    foilvertex.axis1.clear();
  }
  if (zbounds.overlap(checkz)) {
    // axis2 along z
    if (point_plane_check_x(isec3, false))
      foilvertex.axis2.sethigh(checkz.to());
    else 
      foilvertex.axis2.sethigh(zbounds.to());      
    if (point_plane_check_x(isec4, false))
      foilvertex.axis2.setlow(checkz.from());
    else 
      foilvertex.axis2.setlow(zbounds.from());
  }
  else {
    foilvertex.axis2.clear();
  }
  double original_area = fabs(isec2.y() - isec1.y()) * fabs(isec4.z() - isec3.z());
  if (fabs(original_area)>0.0)
    foilvertex.areafraction = (foilvertex.axis1.width() * foilvertex.axis2.width()) / original_area;
  else
    foilvertex.areafraction = 0.0;
  foilvertex.neighbourindex = std::make_pair(-1,-1); // no neighbours
  foilvertex.axis3.clear();
}


double VertexExtrapolator::set_calospot_helix(std::vector<Helix3d>& hc, Plane p, int side)
{
  double area;
  double sum_area = 0.0;
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
    std::cout << "hit side: " << p.side << " id = " << p.planeid << std::endl;
    //    std::cout << "centre at (" << centre.x() << ", " << centre.y() << ", " << centre.z() << ")" << std::endl;
    if (centre.x()>10000.0 || centre.y()>10000.0 || centre.z()>10000.0) {
      //      std::cout << "no intersection, return 0" << std::endl;
      return 0.0;; // no intersection signature, end here
    }
    if (p.planeid<2) { // main wall
      if (point_plane_check_x(centre, true) && point_plane_check_x(centre, false)) { // best fit hits this plane in both axes
	area = mainwall_check_helix(hc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // xwall
	    sum_area = xwall_check_helix(hc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    sum_area = gveto_check_helix(hc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since xwall checked first
	    sum_area = gveto_check_helix(hc, allPlanes.at(indx2), sum_area);
	} 
      }
    }
    else if ((p.planeid>=2 && p.planeid<4) || (p.planeid>=6 && p.planeid<8)) { // xwall
      if (point_plane_check_y(centre, side, true) && point_plane_check_y(centre, side, false)) { // best fit hits this plane in both axes
	area = xwall_check_helix(hc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    sum_area = mainwall_check_helix(hc, allPlanes.at(indx1), area);
	  else if ((indx1>=4 && indx1<6) || (indx1>=8 && indx1<10)) // should be set to gveto then
	    sum_area = gveto_check_helix(hc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be gveto since mainwall checked first
	    sum_area = gveto_check_helix(hc, allPlanes.at(indx2), sum_area);
	}
      }
    }
    else { // gveto
      if (point_plane_check_z(centre, side, true) && point_plane_check_z(centre, side, false)) { // best fit hits this plane in both axes
	area = gveto_check_helix(hc, p, -1.0); // deposits main rectangle in vector
	if (calovertex.back().areafraction<1.0) { // overlaps with another wall
	  // check neighbour walls
	  int indx1 = calovertex.back().neighbourindex.first;
	  int indx2 = calovertex.back().neighbourindex.second;
	  if (indx1<2) // mainwall
	    sum_area = mainwall_check_helix(hc, allPlanes.at(indx1), area);
	  else if ((indx1>=2 && indx1<4) || (indx1>=6 && indx1<8)) // should be set to xwall then
	    sum_area = xwall_check_helix(hc, allPlanes.at(indx1), area);
	  if (indx2 > 0) // can only be xwall since mainwall checked first
	    sum_area = xwall_check_helix(hc, allPlanes.at(indx2), sum_area);
	}
      }
    }
  } // correct side, otherwise do nothing
  //  std::cout << "wrong plane side, returns." << std::endl;
  return sum_area==0.0 ? area : sum_area;
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
  return std::sqrt((p2-p1).Mag2()); // distance
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
  if (fabs(arg)>=1) { // not hitting main wall
    imin.SetXYZ(10001.0, 0.0, 0.0); // empty
    return imin;
  }
  double pitch = h.pitch;
  double t  = std::acos(arg);
  double pi = std::acos(-1.0);
  //  std::cout << "z-inputs: x0-xc= " << (x0-xc) << " arg=" << arg << ", t=" << t << " pitch=" << pitch << std::endl;
  double y  = h.point.y() + std::sqrt(r*r - (x0-xc)*(x0-xc));
  double y2 = h.point.y() - std::sqrt(r*r - (x0-xc)*(x0-xc));
  double z  = h.point.z() + pitch * t  / (2.0 * pi);
  double z2 = h.point.z() - pitch * t  / (2.0 * pi);
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
  if (fabs(arg)>=1) { // not hitting foil
    imin.SetXYZ(10001.0, 0.0, 0.0); // empty
    return imin;
  }

  double pitch = h.pitch;
  double t  =  std::acos(arg);
  double pi = std::acos(-1.0);
  //  std::cout << "foil z-inputs: x0-xc= " << (x0-xc) << " arg=" << arg << ", t=" << t << " pitch=" << pitch << std::endl;

  double y  = h.point.y() + std::sqrt(r*r - (x0-xc)*(x0-xc));
  double y2 = h.point.y() - std::sqrt(r*r - (x0-xc)*(x0-xc));
  double z  = h.point.z() + pitch * t  / (2.0 * pi);
  double z2 = h.point.z() - pitch * t  / (2.0 * pi);
  ROOT::Math::XYZPoint i(x0,y,z);
  pointcollection.push_back(i);
  ROOT::Math::XYZPoint i2(x0,y2,z);  // four piercing points
  pointcollection.push_back(i2);
  ROOT::Math::XYZPoint i3(x0,y,z2);
  pointcollection.push_back(i3);
  ROOT::Math::XYZPoint i4(x0,y2,z2);
  pointcollection.push_back(i4);
  // std::cout << "raw foil hit: (" << i.x() << ", " << i.y() << ", " << i.z() << ")" << std::endl;
  // std::cout << "raw foil hit: (" << i2.x() << ", " << i2.y() << ", " << i2.z() << ")" << std::endl;
  // std::cout << "raw foil hit: (" << i3.x() << ", " << i3.y() << ", " << i3.z() << ")" << std::endl;
  // std::cout << "raw foil hit: (" << i4.x() << ", " << i4.y() << ", " << i4.z() << ")" << std::endl;
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
  double yc = h.point.y(); // helix centre y
  double r  = h.radius;
  double pitch = h.pitch;
  double pi = std::acos(-1.0);
  double arg = (y0 - yc);

  if (fabs(arg) / r >= 1) { // not hitting xwall
    imin.SetXYZ(0.0, 10001.0, 0.0); // empty
    return imin;
  }

  double t  =  std::acos(1/r * std::sqrt(r*r - arg*arg));
  double x  = h.point.x() + std::sqrt(r*r - arg*arg);
  double x2 = h.point.x() - std::sqrt(r*r - arg*arg);
  double z  = h.point.z() + pitch * t / (2.0 * pi);
  double z2 = h.point.z() - pitch * t / (2.0 * pi);
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
  double z0 = p.point.z(); // z coordinate of plane
  double zc = h.point.z(); // helix centre z
  double r  = h.radius;
  double pitch = h.pitch;
  double pi = std::acos(-1.0);
  if (fabs(pitch)>0.0) {
    double t = fabs(z0 - zc) / pitch; // units of 2*pi
    if (t>=1) { // prevent full or more circle spirals to hit z plane
      i.SetXYZ(0.0, 0.0, 10001.0); // empty
      return i;
    }
    double x  = h.point.x() + r * std::cos(2.0 * pi * t);
    double y  = h.point.y() + r * std::sin(2.0 * pi * t);
    
    return i.SetXYZ(x, y, z0);
  }
  else { // parallel to z plane
      i.SetXYZ(0.0, 0.0, 10001.0); // empty
      return i;
  }
}


// Checking section
// ****************
// line checks
double VertexExtrapolator::mainwall_check(std::vector<Line3d>& lc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval a1;
  Interval a2;
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  ROOT::Math::XYZPoint isec1 = intersect_line_plane(lc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_line_plane(lc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_line_plane(lc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_line_plane(lc.at(3), p);
  std::cout << "error main calo hit: (" << isec1.x() << ", " << isec1.y() << ", " << isec1.z() << ")" << std::endl;
  std::cout << "error main calo hit: (" << isec2.x() << ", " << isec2.y() << ", " << isec2.z() << ")" << std::endl;
  std::cout << "error main calo hit: (" << isec3.x() << ", " << isec3.y() << ", " << isec3.z() << ")" << std::endl;
  std::cout << "error main calo hit: (" << isec4.x() << ", " << isec4.y() << ", " << isec4.z() << ")" << std::endl;
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  spot.axis3.clear();
  
  // axis 1 of rectangle
  if (point_plane_check_x(isec1, true)) { // -y error
    (p.side==0) ? a1.setlow(isec1.y()) : a1.sethigh(isec1.y());
    //    std::cout << "mwall isec1 check x passed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==0) ? a1.setlow(ybound.from()) : a1.sethigh(ybound.to());
    if (spot.neighbourindex.first<0) { // identify one of four xwalls
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else { 
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
    //    std::cout << "mainwall isec1 check x failed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  //  std::cout << "spot axis1 part 1: (" << a1.from() << ", " << a1.to() << ")" << std::endl;
  if (point_plane_check_x(isec2, true)) { // +y error
    (p.side==0) ? a1.sethigh(isec2.y()) : a1.setlow(isec2.y());
    //    std::cout << "xwall isec2 check x passed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==0) ? a1.sethigh(ybound.to()) : a1.setlow(ybound.from());
    if (spot.neighbourindex.first<0) { // identify one of four xwalls
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else { 
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
    //    std::cout << "mainwall isec2 check x failed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }  
  //  std::cout << "full spot axis1: (" << a1.from() << ", " << a1.to() << ")" << std::endl;
  
  // axis 2 of rectangle
  if (point_plane_check_x(isec3, false)) // -z error
    (p.side==0) ? a2.setlow(isec3.z()) : a2.sethigh(isec3.z());
  else {
    (p.side==0) ? a2.setlow(zbound.from()) : a2.sethigh(zbound.to());
    if (spot.neighbourindex.first<0) { // identify one of four gveto
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else { 
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
  }
  if (point_plane_check_x(isec4, false)) // +z error
    (p.side==0) ? a2.sethigh(isec4.z()) : a2.setlow(isec4.z());
  else {
    (p.side==0) ? a2.sethigh(zbound.to()) : a2.setlow(zbound.from());
    if (spot.neighbourindex.first<0) { // identify one of four gveto
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else { 
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
  }  
  //  std::cout << "full spot axis2: (" << a2.from() << ", " << a2.to() << ")" << std::endl;

  spot.axis1 = Interval(a1.from(), a1.to());
  spot.axis2 = Interval(a2.from(), a2.to());

  double original_area;
  double sumarea;
  double width1 = a1.width();
  double width2 = a2.width();
  if (area<0.0) { // first time calculation
    original_area = fabs(isec2.y() - isec1.y()) * fabs(isec4.z() - isec3.z());
    sumarea = (width1 * width2);
    spot.areafraction = sumarea / original_area; // is there overlap at all -> less than 1
  }
  else { // previous area calculation received
    sumarea = (width1 * width2) + area;
    spot.areafraction = (width1 * width2) / sumarea; // proper fraction of sum projected area
  }
  calovertex.push_back(spot);
  return sumarea; // that projection area
}


double VertexExtrapolator::xwall_check(std::vector<Line3d>& lc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval a1;
  Interval a2;
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  //  std::cout << "xwall hit side: " << p.side << " id = " << p.planeid << std::endl;
  ROOT::Math::XYZPoint isec1 = intersect_line_plane(lc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_line_plane(lc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_line_plane(lc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_line_plane(lc.at(3), p);
  std::cout << "error xwall hit: (" << isec1.x() << ", " << isec1.y() << ", " << isec1.z() << ")" << std::endl;
  std::cout << "error xwall hit: (" << isec2.x() << ", " << isec2.y() << ", " << isec2.z() << ")" << std::endl;
  std::cout << "error xwall hit: (" << isec3.x() << ", " << isec3.y() << ", " << isec3.z() << ")" << std::endl;
  std::cout << "error xwall hit: (" << isec4.x() << ", " << isec4.y() << ", " << isec4.z() << ")" << std::endl;
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  spot.axis3.clear();
  
  // axis 1 of rectangle
  if (point_plane_check_y(isec1, p.side, true)) { // -y error
    (p.side==1) ? a1.setlow(isec1.x()) : a1.sethigh(isec1.x());
    //    std::cout << "xwall isec1 check y passed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==1) ? a1.setlow(xbound_front.to()) : a1.sethigh(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
    //    std::cout << "xwall isec1 check y failed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  //  std::cout << "xwall spot axis1 part 1: (" << a1.from() << ", " << a1.to() << ")" << std::endl;
  if (point_plane_check_y(isec2, p.side, true)) { // +y error
    (p.side==1) ? a1.sethigh(isec2.x()) : a1.setlow(isec2.x());
    //    std::cout << "xwall isec2 check y passed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==1) ? a1.sethigh(xbound_front.to()) : a1.setlow(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  //  std::cout << "xwall full spot axis1: (" << a1.from() << ", " << a1.to() << ")" << std::endl;
  
  // axis 2 of rectangle
  if (point_plane_check_y(isec3, p.side, false)) { // -z error
    (p.side==0) ? a2.setlow(isec3.z()) : a2.sethigh(isec3.z());
    //    std::cout << "xwall isec3 check y passed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==0) ? a2.setlow(zbound.from()) : a2.sethigh(zbound.to());
    if (spot.neighbourindex.first<0) { // identify one of four gveto
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else { 
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
    //    std::cout << "xwall isec3 check y failed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }

  if (point_plane_check_y(isec4, p.side, false)) { // +z error
    (p.side==0) ? a2.sethigh(isec4.z()) : a2.setlow(isec4.z());
    //    std::cout << "xwall isec4 check y passed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==0) ? a2.sethigh(zbound.to()) : a2.setlow(zbound.from());
    if (spot.neighbourindex.first<0) { // identify one of four gveto
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else { 
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
    //    std::cout << "xwall isec4 check y failed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }
  
  spot.axis1 = Interval(a1.from(), a1.to());
  spot.axis2 = Interval(a2.from(), a2.to());

  double original_area;
  double sumarea;
  double width1 = a1.width();
  double width2 = a2.width();
  if (area<0.0) { // first time calculation
    original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.z() - isec3.z());
    sumarea = (width1 * width2);
    spot.areafraction = sumarea / original_area; // is there overlap at all -> less than 1
  }
  else { // previous area calculation received
    sumarea = (width1 * width2) + area;
    spot.areafraction = (width1 * width2) / sumarea; // proper fraction of sum projected area
  }

  calovertex.push_back(spot);
  return sumarea;
}


double VertexExtrapolator::gveto_check(std::vector<Line3d>& lc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval a1;
  Interval a2;
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  ROOT::Math::XYZPoint isec1 = intersect_line_plane(lc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_line_plane(lc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_line_plane(lc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_line_plane(lc.at(3), p);
  std::cout << "error gveto hit: (" << isec1.x() << ", " << isec1.y() << ", " << isec1.z() << ")" << std::endl;
  std::cout << "error gveto hit: (" << isec2.x() << ", " << isec2.y() << ", " << isec2.z() << ")" << std::endl;
  std::cout << "error gveto hit: (" << isec3.x() << ", " << isec3.y() << ", " << isec3.z() << ")" << std::endl;
  std::cout << "error gveto hit: (" << isec4.x() << ", " << isec4.y() << ", " << isec4.z() << ")" << std::endl;
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  spot.axis3.clear();
  
  // axis 1 of rectangle
  if (point_plane_check_z(isec3, p.side, true)) // -y error
    (p.side==1) ? a1.setlow(isec3.x()) : a1.sethigh(isec3.x());
  else {
    (p.side==1) ? a1.setlow(xbound_front.to()) : a1.sethigh(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_z(isec4, p.side, true)) // +y error
    (p.side==1) ? a1.sethigh(isec4.x()) : a1.setlow(isec4.x());
  else {
    (p.side==1) ? a1.sethigh(xbound_front.to()) : a1.setlow(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_z(isec1, p.side, false)) // -z error
    (p.side==0) ? a2.setlow(isec1.y()) : a2.sethigh(isec1.y());
  else {
    (p.side==0) ? a2.setlow(ybound.from()) : a2.sethigh(ybound.to());
    if (spot.neighbourindex.first<0) { // identify one of four xwalls
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else { 
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
  }

  if (point_plane_check_z(isec2, p.side, false)) // +z error
    (p.side==0) ? a2.sethigh(isec2.y()) : a2.setlow(isec2.y());
  else {
    (p.side==0) ? a2.sethigh(ybound.to()) : a2.setlow(ybound.from());
    if (spot.neighbourindex.first<0) { // identify one of four xwalls
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else { 
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
  }
  
  spot.axis1 = Interval(a1.from(), a1.to());
  spot.axis2 = Interval(a2.from(), a2.to());

  double original_area;
  double sumarea;
  double width1 = a1.width();
  double width2 = a2.width();
  if (area<0.0) { // first time calculation
    original_area = fabs(isec2.y() - isec1.y()) * fabs(isec4.x() - isec3.x());
    sumarea = (width1 * width2);
    spot.areafraction = sumarea / original_area; // is there overlap at all -> less than 1
  }
  else { // previous area calculation received
    sumarea = (width1 * width2) + area;
    spot.areafraction = (width1 * width2) / sumarea; // proper fraction of sum projected area
  }

  calovertex.push_back(spot);
  return sumarea;
}


void VertexExtrapolator::zcheck(std::vector<Line3d>& lc, int side)
{
  Plane pbot = allPlanes.at(4); // same on either side
  Plane ptop = allPlanes.at(5);

  for (Line3d& current : lc) {
    ROOT::Math::XYZPoint intersection_top = intersect_line_plane(current, ptop);
    ROOT::Math::XYZPoint intersection_bot = intersect_line_plane(current, pbot);
    if (intersection_top.x()>1000.0 &&intersection_top.y()>10000.0 &&intersection_top.z()>10000.0) continue; // no intersection signature

    if (point_plane_check_z(intersection_top, side, true) && point_plane_check_z(intersection_top, side, false)) // point overlaps gveto, both axes
      info.foilcalo.second = true; // not on wire in z
    if (point_plane_check_z(intersection_bot, side, true) && point_plane_check_z(intersection_bot, side, false)) // point overlaps gveto
      info.foilcalo.second = true; // not on wire in z
  }
}

// Helix checks
double VertexExtrapolator::mainwall_check_helix(std::vector<Helix3d>& hc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval a1;
  Interval a2;
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  ROOT::Math::XYZPoint isec1 = intersect_helix_plane(hc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_helix_plane(hc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_helix_plane(hc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_helix_plane(hc.at(3), p);
  std::cout << "helix mwall hit: (" << isec1.x() << ", " << isec1.y() << ", " << isec1.z() << ")" << std::endl;
  std::cout << "helix mwall hit: (" << isec2.x() << ", " << isec2.y() << ", " << isec2.z() << ")" << std::endl;
  std::cout << "helix mwall hit: (" << isec3.x() << ", " << isec3.y() << ", " << isec3.z() << ")" << std::endl;
  std::cout << "helix mwall hit: (" << isec4.x() << ", " << isec4.y() << ", " << isec4.z() << ")" << std::endl;
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  spot.axis3.clear();

  // axis 1 of rectangle
  if (point_plane_check_x(isec1, true)) {// -y error
    (p.side==0) ? a1.setlow(isec1.y()) : a1.sethigh(isec1.y());
    //    std::cout << "mwall isec1 check x passed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  else {
    //    (p.side==0) ? a1.setlow(ybound.from()) : a1.sethigh(ybound.to());
    (p.side==0) ? a1.setlow(findNearestBound(ybound, isec1.y())) : a1.sethigh(findNearestBound(ybound, isec1.y()));
    if (spot.neighbourindex.first<0) {  // identify one of four xwalls
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else {
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
    //    std::cout << "mainwall isec1 check x failed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  if (point_plane_check_x(isec2, true)) {// +y error
    (p.side==0) ? a1.sethigh(isec2.y()) : a1.setlow(isec2.y());
    //    std::cout << "mwall isec2 check x passed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==0) ? a1.sethigh(findNearestBound(ybound, isec2.y())) : a1.setlow(findNearestBound(ybound, isec2.y()));
    if (spot.neighbourindex.first<0) {
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else {
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
    //    std::cout << "mainwall isec2 check x failed: (" << a1.from() << ", " << a1.to() << ") side = " << p.side << std::endl;
  }  
  //  std::cout << "spot axis1: (" << a1.from() << ", " << a1.to() << ")" << std::endl;

  // axis 2 of rectangle
  if (point_plane_check_x(isec3, false)) { // -z error
    (p.side==0) ? a2.setlow(isec3.z()) : a2.sethigh(isec3.z());
    //    std::cout << "mwall isec3 check x passed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==0) ? a2.setlow(findNearestBound(zbound, isec3.z())) : a2.sethigh(findNearestBound(zbound, isec3.z()));
    if (spot.neighbourindex.first<0) { // identify one of four gveto
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else {
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
    //    std::cout << "mainwall isec3 check x failed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }
  if (point_plane_check_x(isec4, false)) {// +z error
    (p.side==0) ? a2.sethigh(isec4.z()) : a2.setlow(isec4.z());
    //    std::cout << "mwall isec4 check x passed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }
  else {
    (p.side==0) ? a2.sethigh(findNearestBound(zbound, isec4.z())) : a2.setlow(findNearestBound(zbound, isec4.z()));
    if (spot.neighbourindex.first<0) { spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else {
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
    //    std::cout << "mainwall isec4 check x failed: (" << a2.from() << ", " << a2.to() << ") side = " << p.side << std::endl;
  }  
  //  std::cout << "spot axis2: (" << a2.from() << ", " << a2.to() << ")" << std::endl;

  spot.axis1 = Interval(a1.from(), a1.to());
  spot.axis2 = Interval(a2.from(), a2.to());

  double original_area;
  double sumarea;
  double width1 = a1.width();
  double width2 = a2.width();
  if (area<0.0) { // first time calculation
    original_area = fabs(isec2.y() - isec1.y()) * fabs(isec4.z() - isec3.z());
    sumarea = (width1 * width2);
    spot.areafraction = sumarea / original_area; // is there overlap at all -> less than 1
  }
  else { // previous area calculation received
    sumarea = (width1 * width2) + area;
    spot.areafraction = (width1 * width2) / sumarea; // proper fraction of sum projected area
  }

  calovertex.push_back(spot);
  return sumarea;
}


double VertexExtrapolator::xwall_check_helix(std::vector<Helix3d>& hc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval a1;
  Interval a2;
  Interval zbound(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  ROOT::Math::XYZPoint isec1 = intersect_helix_plane(hc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_helix_plane(hc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_helix_plane(hc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_helix_plane(hc.at(3), p);
  std::cout << "helix xwall hit: (" << isec1.x() << ", " << isec1.y() << ", " << isec1.z() << ")" << std::endl;
  std::cout << "helix xwall hit: (" << isec2.x() << ", " << isec2.y() << ", " << isec2.z() << ")" << std::endl;
  std::cout << "helix xwall hit: (" << isec3.x() << ", " << isec3.y() << ", " << isec3.z() << ")" << std::endl;
  std::cout << "helix xwall hit: (" << isec4.x() << ", " << isec4.y() << ", " << isec4.z() << ")" << std::endl;
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  spot.axis3.clear();
  
  // axis 1 of rectangle
  if (point_plane_check_y(isec1, p.side, true)) // -y error
    (p.side==1) ? a1.setlow(isec1.x()) : a1.sethigh(isec1.x());
  else {
    (p.side==1) ? a1.setlow(xbound_front.to()) : a1.sethigh(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  if (point_plane_check_y(isec2, p.side, true)) // +y error
    (p.side==1) ? a1.sethigh(isec2.x()) : a1.setlow(isec2.x());
  else {
    (p.side==1) ? a1.sethigh(xbound_front.to()) : a1.setlow(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  
  // axis 2 of rectangle
  if (point_plane_check_y(isec3, p.side, false)) // -z error
    (p.side==0) ? a2.setlow(isec3.z()) : a2.sethigh(isec3.z());
  else {
    (p.side==0) ? a2.setlow(findNearestBound(zbound, isec3.z())) : a2.sethigh(findNearestBound(zbound, isec3.z()));
    if (spot.neighbourindex.first<0) { // identify one of four gveto
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else {spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
      if (isec3.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
  }

  if (point_plane_check_y(isec4, p.side, false)) // +z error
    (p.side==0) ? a2.sethigh(isec4.z()) : a2.setlow(isec4.z());
  else {
    (p.side==0) ? a2.sethigh(findNearestBound(zbound, isec4.z())) : a2.setlow(findNearestBound(zbound, isec4.z()));
    if (spot.neighbourindex.first<0) {
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 5 : 9, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 4 : 8, spot.neighbourindex.second);
    }
    else {
      if (isec4.z()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 5 : 9);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 4 : 8);
    }
  }

  spot.axis1 = Interval(a1.from(), a1.to());
  spot.axis2 = Interval(a2.from(), a2.to());
  
  double original_area;
  double sumarea;
  double width1 = a1.width();
  double width2 = a2.width();
  if (area<0.0) { // first time calculation
    original_area = fabs(isec2.x() - isec1.x()) * fabs(isec4.z() - isec3.z());
    sumarea = (width1 * width2);
    spot.areafraction = sumarea / original_area; // is there overlap at all -> less than 1
  }
  else { // previous area calculation received
    sumarea = (width1 * width2) + area;
    spot.areafraction = (width1 * width2) / sumarea; // proper fraction of sum projected area
  }

  calovertex.push_back(spot);
  return sumarea;
}


double VertexExtrapolator::gveto_check_helix(std::vector<Helix3d>& hc, Plane p, double area)
{
  Rectangle spot; // into vector calovertex
  Interval a1;
  Interval a2;
  Interval ybound(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval xbound_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbound_front(0.0, allPlanes.at(1).point.x());

  ROOT::Math::XYZPoint isec1 = intersect_helix_plane(hc.at(0), p);
  ROOT::Math::XYZPoint isec2 = intersect_helix_plane(hc.at(1), p);
  ROOT::Math::XYZPoint isec3 = intersect_helix_plane(hc.at(2), p);
  ROOT::Math::XYZPoint isec4 = intersect_helix_plane(hc.at(3), p);
  std::cout << "helix gveto hit: (" << isec1.x() << ", " << isec1.y() << ", " << isec1.z() << ")" << std::endl;
  std::cout << "helix gveto hit: (" << isec2.x() << ", " << isec2.y() << ", " << isec2.z() << ")" << std::endl;
  std::cout << "helix gveto hit: (" << isec3.x() << ", " << isec3.y() << ", " << isec3.z() << ")" << std::endl;
  std::cout << "helix gveto hit: (" << isec4.x() << ", " << isec4.y() << ", " << isec4.z() << ")" << std::endl;
  spot.planeid = p.planeid;
  spot.side = p.side;
  spot.neighbourindex = std::make_pair(-1,-1);
  spot.axis3.clear();
  
  // axis 1 of rectangle
  if (point_plane_check_z(isec3, p.side, true)) // -y error
    (p.side==1) ? a1.setlow(isec3.x()) : a1.sethigh(isec3.x());
  else {
    (p.side==1) ? a1.setlow(xbound_front.to()) : a1.sethigh(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  //  std::cout << "after isec3 check: a1=[" << a1.from() << ", " << a1.to() << "]" << std::endl;
  if (point_plane_check_z(isec4, p.side, true)) // +y error
    (p.side==1) ? a1.sethigh(isec4.x()) : a1.setlow(isec4.x());
  else {
    (p.side==1) ? a1.sethigh(xbound_front.to()) : a1.setlow(xbound_back.from());
    if (spot.neighbourindex.first<0) spot.neighbourindex = std::make_pair((p.side==1) ? 1 : 0, spot.neighbourindex.second);
    else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==1) ? 1 : 0);
  }
  //  std::cout << "after isec4 check: a1=[" << a1.from() << ", " << a1.to() << "]" << std::endl;
  
  // axis 2 of rectangle
  if (point_plane_check_z(isec1, p.side, false)) // -z error
    (p.side==0) ? a2.setlow(isec1.y()) : a2.sethigh(isec1.y());
  else {
    (p.side==0) ? a2.setlow(findNearestBound(ybound, isec1.y())) : a2.sethigh(findNearestBound(ybound, isec1.y()));
    if (spot.neighbourindex.first<0) {
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else { 
      if (isec1.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
  }

  if (point_plane_check_z(isec2, p.side, false)) // +z error
    (p.side==0) ? a2.sethigh(isec2.y()) : a2.setlow(isec2.y());
  else {
    (p.side==0) ? a2.sethigh(findNearestBound(ybound, isec2.y())) : a2.setlow(findNearestBound(ybound, isec2.y()));
    if (spot.neighbourindex.first<0) {
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair((p.side==0) ? 3 : 7, spot.neighbourindex.second);
      else spot.neighbourindex = std::make_pair((p.side==0) ? 2 : 6, spot.neighbourindex.second);
    }
    else { 
      if (isec2.y()>0.0) spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 3 : 7);
      else spot.neighbourindex = std::make_pair(spot.neighbourindex.first, (p.side==0) ? 2 : 6);
    }
  }
      
  spot.axis1 = Interval(a1.from(), a1.to());
  spot.axis2 = Interval(a2.from(), a2.to());
  
  double original_area;
  double sumarea;
  double width1 = a1.width();
  double width2 = a2.width();
  if (area<0.0) { // first time calculation
    original_area = fabs(isec2.y() - isec1.y()) * fabs(isec4.x() - isec3.x());
    sumarea = (width1 * width2);
    spot.areafraction = sumarea / original_area; // is there overlap at all -> less than 1
  }
  else { // previous area calculation received
    sumarea = (width1 * width2) + area;
    spot.areafraction = (width1 * width2) / sumarea; // proper fraction of sum projected area
  }
  calovertex.push_back(spot);
  return sumarea;
}


void VertexExtrapolator::zcheck_helix(std::vector<Helix3d>& hc, int side)
{
  Plane pbot = allPlanes.at(4); // same on either side
  Plane ptop = allPlanes.at(5);

  for (Helix3d& current : hc) {
    ROOT::Math::XYZPoint intersection_top = intersect_helix_plane(current, ptop);
    ROOT::Math::XYZPoint intersection_bot = intersect_helix_plane(current, pbot);
    if (intersection_top.z() < 10001.0) { // intersection signature
      if (point_plane_check_z(intersection_top, side, true) && point_plane_check_z(intersection_top, side, false)) // point overlaps gveto both axes
	info.foilcalo.second = true; // not on wire in z
    }
    else if (intersection_bot.z()<10001.0) { // intersection signature
      if (point_plane_check_z(intersection_bot, side, true) && point_plane_check_z(intersection_bot, side, false)) // point overlaps gveto
	info.foilcalo.second = true; // not on wire in z
    }
  }
}

// Point checking
bool VertexExtrapolator::point_plane_check_x(ROOT::Math::XYZPoint point, bool axis1)
{
  // within bounds of main wall and correct side.
  Interval ybounds(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());
  Interval zbounds(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  bool check = false;
  if (axis1 && ybounds.contains(point.y()))
    check = true;
  if (!axis1 && zbounds.contains(point.z()))
    check = true;

  return check;
}

bool VertexExtrapolator::point_plane_check_y(ROOT::Math::XYZPoint point, int side, bool axis1)
{
  // within bounds of xwall and correct side.
  Interval xbounds_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbounds_front(0.0, allPlanes.at(1).point.x());
  Interval zbounds(allPlanes.at(4).point.z(), allPlanes.at(5).point.z());

  bool check = false;
  if (side ==0) { // checks in order of axis consignments to x, y and z axes
    if (axis1 && xbounds_back.contains(point.x()))
      check = true;
    if (!axis1  && zbounds.contains(point.z()))
      check = true;
  }
  if (side ==1) {
    if (axis1 && xbounds_front.contains(point.x()))
      check = true;
    if (!axis1  && zbounds.contains(point.z()))
      check = true;
  }
  return check;
}


bool VertexExtrapolator::point_plane_check_z(ROOT::Math::XYZPoint point, int side, bool axis1)
{
  // within bounds of gveto and correct side.
  Interval xbounds_back(allPlanes.at(0).point.x(), 0.0);
  Interval xbounds_front(0.0, allPlanes.at(1).point.x());
  Interval ybounds(allPlanes.at(2).point.y(), allPlanes.at(3).point.y());

  bool check = false;
  if (side ==0) { // checks in order of axis consignments to x, y and z axes
    if (axis1 && xbounds_back.contains(point.x()))
      check = true;
    if (!axis1  && ybounds.contains(point.y()))
      check = true;
  }
  if (side ==1) {
    if (axis1 && xbounds_front.contains(point.x()))
      check = true;
    if (!axis1  && ybounds.contains(point.y()))
      check = true;
  }
  return check;
}

// for charge orientation
double VertexExtrapolator::findLowerYBoundatEnds() {
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

double VertexExtrapolator::findUpperYBoundatEnds() {
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

// for wire vertexing
void VertexExtrapolator::set_wirevertex(bool foilside) { // which cell contains the wire vertex
  wirevertex.planeid = -1; // any plane as required = a box
  wirevertex.side = info.side; // where the cluster is
  wirevertex.neighbourindex = std::make_pair(-1,-1); // none
  wirevertex.areafraction = 1.0; // by definition
  wirevertex.axis3.clear();

  if (foilside) {
    MetaInfo mi = findFoilSideCorner();
    wirevertex.axis1.setlow(mi.wirex - 22.0); // x [mm]
    wirevertex.axis1.sethigh(mi.wirex + 22.0); // x
    wirevertex.axis2.setlow(mi.wirey - 22.0); // y
    wirevertex.axis2.sethigh(mi.wirey + 22.0); // y
    wirevertex.axis3.setlow(mi.zcoord - 10.0); // z
    wirevertex.axis3.sethigh(mi.zcoord + 10.0); // z
  }
  else {
    MetaInfo mi = findCaloSideCorner();
    wirevertex.axis1.setlow(mi.wirex - 22.0); // x
    wirevertex.axis1.sethigh(mi.wirex + 22.0); // x
    wirevertex.axis2.setlow(mi.wirey - 22.0); // y
    wirevertex.axis2.sethigh(mi.wirey + 22.0); // y
    wirevertex.axis3.setlow(mi.zcoord - 10.0); // z
    wirevertex.axis3.sethigh(mi.zcoord + 10.0); // z
  }
}


double VertexExtrapolator::findNearestBound(Interval bound, double value) {
  double diff1 = fabs(bound.from() - value);
  double diff2 = fabs(bound.to() - value);
  return (diff1 <= diff2) ? bound.from() : bound.to();
}


MetaInfo VertexExtrapolator::findFoilSideCorner() {
  // only runs on cells not on layer 0 but closest to foil
  if (info.minx.size() == 1) return info.minx.front(); // just the one possibility

  std::vector<double> dummy;
  for (auto val : info.minx) // foil side
    dummy.push_back(val.wirey); // check y coordinate
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  MetaInfo minmi = info.minx.at(minit - dummy.begin());
  MetaInfo maxmi = info.minx.at(maxit - dummy.begin());
  return findCorner(minmi, maxmi);
}


MetaInfo VertexExtrapolator::findCaloSideCorner() {
  // only runs on cells not on layer 8
  if (info.maxx.size() == 1) return info.maxx.front(); // just the one possibility

  std::vector<double> dummy;
  for (auto val : info.maxx) // calo side
    dummy.push_back(val.wirey); // check y coordinate
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  MetaInfo minmi = info.minx.at(minit - dummy.begin());
  MetaInfo maxmi = info.minx.at(maxit - dummy.begin());
  return findCorner(minmi, maxmi);
}


MetaInfo VertexExtrapolator::findCorner(MetaInfo minimum, MetaInfo maximum) {
  // which direction? Straight tendency up or down
  double ybound_up  = findUpperYBoundatEnds(); // box limit in y at ends
  double ybound_low = findLowerYBoundatEnds();

  // bending curves
  std::vector<double> dummy;
  for (auto val : info.miny) // in y for all x, not just end layers
    dummy.push_back(val.wirey); // check y coordinate
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  double miny = *minit;
  dummy.clear();
  for (auto val : info.maxy) // in y for all x, not just end layers
    dummy.push_back(val.wirey); // check y coordinate
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  double maxy = *maxit;

  // bending
  if (miny < ybound_low) return maximum; // track going down then up
  else if (maxy > ybound_up) return minimum; // track going up then down
  // straight tendency
  else if (minimum.wirey <= ybound_low) return minimum; // found the extreme end
  else return maximum; // found the extreme end
}


// *****
// ** Utility Interval methods
// *****
Interval::Interval()
{
  lower = 0.0;
  upper = 0.0;
  swap = false;
}


Interval::Interval(double s, double e)
{
  if (s <= e) { // start should be smaller than end
    lower = s;
    upper = e;
    swap = false;
  }
  else { // swap order
    lower = e;
    upper = s;
    swap = true;
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

