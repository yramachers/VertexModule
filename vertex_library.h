#ifndef YR_SNVertex
#define YR_SNVertex

// std libraries
#include <vector>
#include <utility>

// ROOT includes
#include <Math/Point3D.h>
#include <Math/Vector3D.h>


// vertex info storage
struct VertexInfo {
  // how to extrapolate
  std::pair<bool, bool> foilcalo;
  bool wire_candidate;
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
  Interval axis1;
  Interval axis2;
  LineFit lf;
  HelixFit hf;
  BrokenLineFit blf;
  std::vector<VertexInfo> allinfo;
  std::vector<Plane> all_planes;

protected:
  bool zcheck(); // consider wire candidate in z to gveto, still a wire candidate?
  void intersect(); // action, checks on valid plane and trajectory
  int find_clid(int id);

public:
  
  VertexExtrapolator(); // Default Constructor
  VertexExtrapolator(std::vector<Plane> pl); // main Constructor
  ~VertexExtrapolator() {all_planes.clear();}
  
  // signal data input and run all intersections
  void setTrajectory(LineFit dummy, std::vector<VertexInfo> vi) {lf = dummy; allinfo = vi; intersect();}
  void setTrajectory(HelixFit dummy, std::vector<VertexInfo> vi) {hf = dummy; allinfo = vi; intersect();}
  void setTrajectory(BrokenLineFit dummy, std::vector<VertexInfo> vi) {blf = dummy; allinfo = vi; intersect();}

  // result of intersection
  std::pair<VertexInfo, std:pair<Interval, Interval> > fullvertex();
  // could be empty axes for no intersection
};



#endif
