#ifndef YR_SNVertex
#define YR_SNVertex

// std libraries
#include <vector>
#include <utility>

// ROOT includes
#include <Math/Point3D.h>
#include <Math/Vector3D.h>

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
  double width() {return fabs(upper - lower);} // how wide
  bool empty() {return lower == upper;} // check for empty interval
  bool overlap(Interval other); // return true if overlap exists
  bool contains(double value); // value in interval?
  void clear() {lower=0.0; upper = 0.0;} // set empty
  void setbound(double value);
};


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
  // making sense of intersection points for helices
  std::vector<MetaInfo> minx;
  std::vector<MetaInfo> maxx;
  std::vector<MetaInfo> miny;
  std::vector<MetaInfo> maxy;
  // id to link to trajectory solution
  int clsid;
  int side; // tracker side from cluster
};


// geometric plane struct
struct Plane {
  // consider a plane as normal vector and a point
  ROOT::Math::XYZVector normal;
  ROOT::Math::XYZPoint point; // cartesian coordinates assumed
  int planeid; // from position in vector
  int side; // flags x negative (0) or positive (1)
};


// geometric rectangle struct
struct Rectangle {
  Interval axis1; // through centre
  Interval axis2; // perpendicular
  double areafraction; // on given plane fraction; interval ]0,1]
  int planeid; // link to plane
  int side; // flags x negative (0) or positive (1)
  std::pair<int, int> neighbourindex; // indices for fixed allPlanes vector: -1 = none
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

class VertexExtrapolator
{
  // input some trajectory and retrieve an intersection rectangle
  // in form of two intervals as main axes

 private:
  Rectangle foilvertex;
  std::vector<Rectangle> calovertex;
  LineFit lf;
  HelixFit hf;
  BrokenLineFit blf;
  VertexInfo info;
  std::vector<VertexInfo> allInfo;
  std::vector<Plane> allPlanes;

  void intersect_line();
  void intersect_helix();
  void intersect_brokenline();
  void zcheck(std::vector<Line3d>& lc, int side);
  void zcheck_helix(std::vector<Helix3d>& hc, int side);
  void set_calospot(std::vector<Line3d>& lc, Plane p, int side);
  void set_calospot_helix(std::vector<Helix3d>& hc, Plane p, int side);
  bool point_plane_check_x(ROOT::Math::XYZPoint point, int side);
  bool point_plane_check_y(ROOT::Math::XYZPoint point, int side);
  bool point_plane_check_z(ROOT::Math::XYZPoint point, int side);
  double findLowerYBound();
  double findUpperYBound();
  double Pointdistance(ROOT::Math::XYZPoint p1, ROOT::Math::XYZPoint p2);
  double mainwall_check(std::vector<Line3d>& lc, Plane p, double area);
  double xwall_check(std::vector<Line3d>& lc, Plane p, double area);
  double gveto_check(std::vector<Line3d>& lc, Plane p, double area);
  double mainwall_check_helix(std::vector<Helix3d>& hc, Plane p, double area);
  double xwall_check_helix(std::vector<Helix3d>& hc, Plane p, double area);
  double gveto_check_helix(std::vector<Helix3d>& hc, Plane p, double area);
  ROOT::Math::XYZPoint intersect_helix_mainw(Helix3d& h, Plane p);
  ROOT::Math::XYZPoint intersect_helix_xwall(Helix3d& h, Plane p);
  ROOT::Math::XYZPoint intersect_helix_gveto(Helix3d& h, Plane p);
  ROOT::Math::XYZPoint intersect_line_plane(Line3d& l, Plane p);
  ROOT::Math::XYZPoint intersect_helix_plane(Helix3d& h, Plane p);
  std::vector<Line3d> linecollection(int side);
  std::vector<Helix3d> helixcollection(int charge);


protected:
  void intersect(int which); // action, checks on valid plane and trajectory


public:
  
  VertexExtrapolator(); // Default Constructor
  VertexExtrapolator(std::vector<Plane> pl); // main Constructor
  ~VertexExtrapolator() {allPlanes.clear(); allInfo.clear(); calovertex.clear();}
  
  // signal data input and run all intersections
  void setTrajectory(LineFit dummy, std::vector<VertexInfo> vi) {lf = dummy; allInfo = vi; intersect(0);}
  void setTrajectory(HelixFit dummy, std::vector<VertexInfo> vi) {hf = dummy; allInfo = vi; intersect(1);}
  void setTrajectory(BrokenLineFit dummy, std::vector<VertexInfo> vi) {blf = dummy; allInfo = vi; intersect(2);}

  // Results from here
  Rectangle onfoil() {return foilvertex;} // rectangle struct on foil, empty axes if none
  std::vector<Rectangle> oncalo() {return calovertex;} // up to three vertices possible on calo walls

};



#endif
