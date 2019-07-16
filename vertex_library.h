#ifndef YR_SNVertex
#define YR_SNVertex

// std libraries
#include <vector>
#include <utility>

// ROOT includes
#include <Math/Point3D.h>
#include <Math/Vector3D.h>

// geiger data store
struct MetaInfo {
  // geiger hit info
  int hitid;
  int side; 
  int row; 
  int column;
  double wirex;
  double wirey;
  double zcoord;
};

// vertex info storage
struct VertexInfo {
  // how to extrapolate
  std::pair<bool, bool> foilcalo;
  bool wire_candidate;
  // making sense of intersection points for helices
  std::vector<MetaInfo> minx
  std::vector<MetaInfo> maxx
  // id to link to trajectory solution
  int clsid;
  // id to link to plane
  int planeid;
};


// geometric plane struct
struct Plane {
  // consider a plane as normal vector and a point
  ROOT::Math::XYZVector normal;
  ROOT::Math::XYZPoint point; // cartesian coordinates assumed
  int planeid; // from position in vector
  int side; // flags x negative (0) or positive (1)
};


// geometric ellipse struct
struct Ellipse {
  Interval axis1;
  Interval axis2;
  int planeid; // link to plane
  int side; // flags x negative (0) or positive (1)
};


// geometric line3d struct
struct Line3d {
  ROOT::Math::XYZVector u;
  ROOT::Math::XYZPoint point; // cartesian coordinates assumed
};


// geometric helix3d struct
struct Helix3d {
  ROOT::Math::XYZPoint point; // cartesian coordinates assumed
  double radius;
  double pitch;
  int charge;
};


// *** repeat data model used for trajectories
// path point store
struct PathPoint {
  // just a 3D point with errors for simple fitting
  std::pair<int,size_t> pointid;
  double xc;
  double yc;
  double zc;
  double errx;
  double erry;
  double errz;
};

// line fit store, 4 parameter with errors
struct LineFit {
  // par
  double ixy;
  double slxy;
  double ixz;
  double slxz;
  double errixy;
  double errslxy;
  double errixz;
  double errslxz;
  // fit diagnostics
  double chi2;
  double prob;
  int status;
  int clid;
};

// helix fit store, 5 parameter with errors
struct HelixFit {
  // par
  double radius;
  double pitch;
  double xc;
  double yc;
  double zc;
  double raderr;
  double errpitch;
  double errxc;
  double erryc;
  double errzc;
  // fit diagnostics
  double chi2;
  double prob;
  int status;
  int clid;
};

// broken line fit store, 4 parameter with errors and break
struct BrokenLineFit {
  LineFit linefit1; // has its own diagnostics
  LineFit linefit2; // could be empty
  std::vector<int> breakpoints; // first and last of interest
  // fit diagnostics for full BL fit
  std::vector<double> angles;
  std::vector<PathPoint> path;
  double length; // path length in Euclidean metric
  double chi2;
  double prob;
  int status;
  int clid;
};

class Interval
{
  // simple interval mostly for convenience
  // on checking pairs of doubles
  // and interval overlap
private:
  double lower;
  double upper;
  
  
protected:

public:
  
  Interval(); // Default Constructor, not used
  Interval(double s, double e); // Constructor with lower and upper limit
  
  double midinterval() {return 0.5*(lower+upper);} // mean interval value
  double from() {return lower;} // boundary return
  double to() {return upper;} // boundary return
  bool empty() {return lower == upper;} // check for empty interval
  bool overlap(Interval other); // return true if overlap exists
    
};


class VertexExtrapolator
{
  // input some trajectory and retrieve an intersection ellipse
  // in form of two intervals as main axes

 private:
  Ellipse spot1;
  Ellipse spot2;
  LineFit lf;
  HelixFit hf;
  BrokenLineFit blf;
  VertexInfo info;
  std::vector<VertexInfo> allinfo;
  std::vector<Plane> allPlanes;
  void intersect_line();
  void intersect_helix();
  void intersect_brokenline();


protected:
  void zcheck_line(); // consider wire candidate in z to gveto, still a wire candidate?
  void zcheck_helix(); // consider wire candidate in z to gveto, still a wire candidate?
  void intersect(int which); // action, checks on valid plane and trajectory


public:
  
  VertexExtrapolator(); // Default Constructor
  VertexExtrapolator(std::vector<Plane> pl); // main Constructor
  ~VertexExtrapolator() {allPlanes.clear(); allInfo.clear();}
  
  // signal data input and run all intersections
  void setTrajectory(LineFit dummy, std::vector<VertexInfo> vi) {lf = dummy; allInfo = vi; intersect(0);}
  void setTrajectory(HelixFit dummy, std::vector<VertexInfo> vi) {hf = dummy; allInfo = vi; intersect(1);}
  void setTrajectory(BrokenLineFit dummy, std::vector<VertexInfo> vi) {blf = dummy; allInfo = vi; intersect(2);}

  // result of intersection, first foil spot, second calo spot
  std::pair<VertexInfo, std:pair<Ellipse, Ellipse> > fullvertex();
  // could be empty axes for no intersection
};



#endif
