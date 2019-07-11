#ifndef YR_SNVertex
#define YR_SNVertex

// std libraries
#include <vector>
#include <utility>

// ROOT includes
#include <Math/Vector3D.h>
#include <Math/Vector2D.h>
#include <TMath.h>
#include <TVector3.h>


// vertex info storage
struct VertexInfo {
  // how to extrapolate
  std::pair<bool, bool> leftright;
  bool wire_candidate;
  // id to link to trajectory solution
  int clsid;
};

// *** repeat data model used for trajectories

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
  
  Interval(); // Default Constructor
  Interval(double s, double e); // Constructor with lower and upper limit
  
  double midinterval() {return 0.5*(lower+upper);} // mean interval value
  double from() {return lower;} // boundary return
  double to() {return upper;} // boundary return
  bool empty() {return lower == upper;} // check for empty interval
  bool overlap(Interval other); // return true if overlap exists
    
};



#endif
