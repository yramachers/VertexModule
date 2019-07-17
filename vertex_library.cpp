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

// results from this method
std::pair<VertexInfo, std:pair<Ellipse, Ellipse> > VertexExtrapolator::fullvertex();
{
  if (allInfo.empty()) {
    std::cout << "VertexExtrapolator: no vertex info supplied - use setTrajectory first" << std::endl;
    return std::make_pair(VertexInfo, std::make_pair(Ellipse, Ellipse)); // all empty return
  }

  return std::make_pair(info, std::make_pair(spot1, spot2)); // should be filled by intersect function
}



bool VertexExtrapolator::zcheck_line()
{
}

bool VertexExtrapolator::zcheck_helix()
{
}

void VertexExtrapolator::intersect(int which) // main working method
{
  switch(which) {
  case 0: // line case
    info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [](const VertexInfo& vi){return vi.clsid == lf.clid}) - allInfo.begin()); // linefit case
    if (info.wire_candidate)
      zcheck_line();
    intersect_line(); // sets all results
    break;

  case 1: // helix case
    info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [](const VertexInfo& vi){return vi.clsid == hf.clid}) - allInfo.begin()); // helixfit case
    if (info.wire_candidate)
      zcheck_helix();
    intersect_helix(); // sets all results
    break;

  case 2: // broken line case
    info = allInfo.at(std::find_if(allInfo.begin(), allInfo.end(), [](const VertexInfo& vi){return vi.clsid == blf.clid}) - allInfo.begin()); // brokenlinefit case
    if (info.wire_candidate)
      zcheck_line();
    intersect_brokenline(); // sets all results
    break;
  }
  return;
}

// Line section
// ************
void VertexExtrapolator::intersect_line()
{
  std::vector<ROOT::Math::XYZPoint> isecs;
  std::vector<int> idcollection;
  Line3d current;

  // central line
  ROOT::Math::XYZVector ubest(1.0, lf.slxy, lf.slxz);
  ROOT::Math::XYZVector pbest(0.0, lf.ixy, lf.ixz);
  current.u = ubest;
  current.point = pbest;
  if (info.foilcalo.second) { // only if a calo vertex is required
    for (Plane p : allPlanes) {
      ROOT::Math::XYZPoint intersection = intersect_line_plane(current, p);
      if (point_plane_check(intersection)) {
	isecs.push_back(intersection);
	idcollection.push_back(p.planeid);
      }
    }
  }
  if (info.foilcalo.first) { // only if a foil vertex is required
    ROOT::Math::XYZVector norm(1.0,0.0,0.0);
    ROOT::Math::XYZVector origin(0.0,0.0,0.0);
    Plane p;
    p.planeid = 10; // foil plane
    p.normal = norm;
    p.point  = origin;
    ROOT::Math::XYZPoint intersection = intersect_line_plane(current, p);
    if (point_plane_check(intersection)) {
      isecs.push_back(intersection);
      idcollection.push_back(p.planeid);
    }
  }
}


ROOT::Math::XYZPoint VertexExtrapolator::intersect_line_plane(Line3d l, Plane p)
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
    ROOT::Math::XYZPoint i(1001.0,0.0,0.0); // no intersection signature x>1000
    return i;
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
  int side = info.minx.front().side; // side from geiger hits (all on the same side)
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
    isec = intersect_helix_mainw(h, p); // side from p.side
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
bool VertexExtrapolator::point_plane_check(ROOT::Math::XYZPoint point)
{
  // within bounds of tracker and correct side.
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


