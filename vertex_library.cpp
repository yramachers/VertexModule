// us
#include <vertex_library.h>

// ROOT includes


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
  for (Plane& p : allPlanes) {
    ROOT::Math::XYZPoint intersection = intersect_line_plane(current, p);
    if (point_plane_check(intersection)) {
      isecs.push_back(intersection);
      idcollection.push_back(p.planeid);
    }
  }
}



void VertexExtrapolator::intersect_helix()
{
}

void VertexExtrapolator::intersect_brokenline()
{
}

ROOT::Math::XYZPoint VertexExtrapolator::intersect_line_plane(Line3d l, Plane p)
{
  double denom = p.normal.Dot(l.u);
  if (denom>0.0) { // not parallel
    double u = (p.normal.Dot(p.point) - p.normal.Dot(l.point)) / denom;
    return ROOT::Math::XYZPoint i(l.point.x() + u * l.u.x(),
				  l.point.y() + u * l.u.y(),
				  l.point.z() + u * l.u.z());
  }
  else 
    return ROOT::Math::XYZPoint i(1001.0,0.0,0.0); // no intersection signature x>1000
}


 bool VertexExtrapolator::point_plane_check(ROOT::Math::XYZPoint point)
{
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


