#include "catch.hpp"
#include <vertex_library.h>
#include <vector>
#include <iostream>

std::vector<Plane> make_planes() {
  // hard-coded plane geometry, avoiding falaise locators
  int whichcalo = 0; // index position in vector
  int side = 0; // back
  std::vector<Plane> planes; // make a total of ten planes to potentially intersect

  // main calo walls
  const double xcalo_bd[2] = {-435.0, 435.0}; // [mm]

  ROOT::Math::XYZPoint p(xcalo_bd[0],0.0,0.0); // negative x
  ROOT::Math::XYZVector norm(1.0,0.0,0.0); // looking in pos x direction
  Plane mcalo_negx;
  mcalo_negx.planeid = whichcalo;
  mcalo_negx.side = side;
  mcalo_negx.normal = norm;
  mcalo_negx.point = p;
  planes.push_back(mcalo_negx);
  
  ROOT::Math::XYZPoint p2(xcalo_bd[1],0.0,0.0); // positive x
  ROOT::Math::XYZVector norm2(-1.0,0.0,0.0); // looking in neg x direction
  Plane mcalo_posx;
  mcalo_posx.planeid = ++whichcalo;
  mcalo_posx.side = 1; // positive x
  mcalo_posx.normal = norm2;
  mcalo_posx.point = p2;
  planes.push_back(mcalo_posx);

  // xwalls part 1
  const double ycalo_bd[2] = {-2505.5, 2505.5};
  ROOT::Math::XYZPoint p3(0.0,ycalo_bd[0],0.0); // negative y
  ROOT::Math::XYZVector norm3(0.0,1.0,0.0); // looking in pos y direction
  ROOT::Math::XYZPoint p4(0.0,ycalo_bd[1],0.0); // positive y
  ROOT::Math::XYZVector norm4(0.0,-1.0,0.0); // looking in neg y direction
  Plane xwall_backl;
  xwall_backl.planeid = ++whichcalo;
  xwall_backl.side = side;
  xwall_backl.normal = norm3;
  xwall_backl.point = p3;

  Plane xwall_backr;
  xwall_backr.planeid = ++whichcalo;
  xwall_backr.side = side;
  xwall_backr.normal = norm4;
  xwall_backr.point = p4;

  planes.push_back(xwall_backl);
  planes.push_back(xwall_backr);

  // gamma veto part 1
  const double zcalo_bd[2] = {-1550.0, 1550.0};
  ROOT::Math::XYZPoint p5(0.0,0.0,zcalo_bd[0]); // negative z
  ROOT::Math::XYZVector norm5(0.0,0.0,1.0); // looking in pos z direction
  ROOT::Math::XYZPoint p6(0.0,0.0,zcalo_bd[1]); // positive z
  ROOT::Math::XYZVector norm6(0.0,0.0,-1.0); // looking in neg z direction
  Plane gv_backb;
  gv_backb.planeid = ++whichcalo;
  gv_backb.side = side;
  gv_backb.normal = norm5;
  gv_backb.point = p5;
  Plane gv_backt;
  gv_backt.planeid = ++whichcalo;
  gv_backt.side = side;
  gv_backt.normal = norm6;
  gv_backt.point = p6;

  planes.push_back(gv_backb);
  planes.push_back(gv_backt);

  // xwalls part 2
  side  = 1; // front
  ROOT::Math::XYZPoint p7(0.0,ycalo_bd[0],0.0); // negative y
  ROOT::Math::XYZVector norm7(0.0,1.0,0.0); // looking in pos y direction
  ROOT::Math::XYZPoint p8(0.0,ycalo_bd[1],0.0); // positive y
  ROOT::Math::XYZVector norm8(0.0,-1.0,0.0); // looking in neg y direction
  Plane xwall_frontl;
  xwall_frontl.planeid = ++whichcalo;
  xwall_frontl.side = side;
  xwall_frontl.normal = norm7;
  xwall_frontl.point = p7;
  Plane xwall_frontr;
  xwall_frontr.planeid = ++whichcalo;
  xwall_frontr.side = side;
  xwall_frontr.normal = norm8;
  xwall_frontr.point = p8;

  planes.push_back(xwall_frontl);
  planes.push_back(xwall_frontr);

  // gamma veto part 2
  ROOT::Math::XYZPoint p9(0.0,0.0,zcalo_bd[0]); // negative z
  ROOT::Math::XYZVector norm9(0.0,0.0,1.0); // looking in pos z direction
  ROOT::Math::XYZPoint p10(0.0,0.0,zcalo_bd[1]); // positive z
  ROOT::Math::XYZVector norm10(0.0,0.0,-1.0); // looking in neg z direction
  Plane gv_frontb;
  gv_frontb.planeid = ++whichcalo;
  gv_frontb.side = side;
  gv_frontb.normal = norm9;
  gv_frontb.point = p9;
  Plane gv_frontt;
  gv_frontt.planeid = ++whichcalo;
  gv_frontt.side = side;
  gv_frontt.normal = norm10;
  gv_frontt.point = p10;
  
  planes.push_back(gv_frontb);
  planes.push_back(gv_frontt);
  return planes;
}


void printve(VertexExtrapolator& ve) {
  // Results
  Rectangle fvertex = ve.onfoil();
  std::cout << "Rectangle on foil: 1:[" << fvertex.axis1.from() << ", " << fvertex.axis1.to() << "]; 2: [" << fvertex.axis2.from() << ", " << fvertex.axis2.to() << "];" << std::endl;
  std::cout << "area fraction: " << fvertex.areafraction << " on plane " << fvertex.planeid << " neighbours: (" << fvertex.neighbourindex.first << ", " << fvertex.neighbourindex.second << ")" << std::endl;

  std::vector<Rectangle> wvertex = ve.onwire();
  std::cout << "Rectangle on calo: " << std::endl;
  for (auto& wv : wvertex) {
    std::cout << "Rectangle on wire: 1:[" << wv.axis1.from() << ", " << wv.axis1.to() << "]; 2: [" << wv.axis2.from() << ", " << wv.axis2.to() << "];" << std::endl;
    std::cout << "Rectangle on wire: 3:[" << wv.axis3.from() << ", " << wv.axis3.to() << "]" << std::endl;
    std::cout << "area fraction: " << wv.areafraction << " on plane " << wv.planeid << " neighbours: (" << wv.neighbourindex.first << ", " << wv.neighbourindex.second << ")" << std::endl;
  }

  std::vector<Rectangle> cvertex = ve.oncalo(); // should be just one entry
  std::cout << "Rectangle on calo: " << std::endl;
  for (auto& cv : cvertex) {
    std::cout << "on calo: 1:[" << cv.axis1.from() << ", " << cv.axis1.to() << "]; 2: [" << cv.axis2.from() << ", " << cv.axis2.to() << "];" << std::endl;
    std::cout << "area fraction: " << cv.areafraction << " on plane " << cv.planeid << " neighbours: (" << cv.neighbourindex.first << ", " << cv.neighbourindex.second << ")" << std::endl;
  }
}


void posxhits_cside(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // on side = 1
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 1; // a minimum geiger hit
  mi.wirex = 74.0; // not quite but irrelevant here
  mi.wirey = 0.0; // central
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // on side = 1
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = 412.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


void posxhits_fside(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // on side = 1
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // central
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // on side = 1
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 7; // a maximum geiger hit
  mi.wirex = 368.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


void negxhits_cside(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 1; // a minimum geiger hit
  mi.wirex = -74.0; // not quite but irrelevant here
  mi.wirey = 0.0; // central
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


void negxhits_fside(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // central
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 7; // a maximum geiger hit
  mi.wirex = -368.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


LineFit lineA() { // artificial line fit solution
  LineFit lf;
  lf.ixy = 0.0;
  lf.ixz = 0.0; // intercept origin
  lf.slxy = 1.0;
  lf.slxz = 0.0; // diagonal up in y, flat in z
  lf.errixy = 0.5;
  lf.errixz = 0.5; // half in intercept y, z
  lf.errslxy = 0.2; 
  lf.errslxz = 0.2; // some error in y, z slopes
  lf.status = 0;
  lf.clid = 1;
  return lf;
}


int check_wirevertex(int which) {
  VertexExtrapolator ve(make_planes());
  LineFit lf1 = lineA();

  VertexInfo vi;
  vi.clsid = 1; // same clsid as hook

  if (which==0) {
    vi.side = 1; // positive x
    // set up single wirevertex case
    vi.foilcalo = std::make_pair(false, true); // foil wire vertices
    posxhits_cside(vi);
  }
  else if (which==1) {
    vi.side = 1; // positive x
    // set up single wirevertex case
    vi.foilcalo = std::make_pair(true, false); // calo wire vertices
    posxhits_fside(vi);
  }
  else if (which==2) {
    vi.side = 0; // negative x
    // set up single wirevertex case
    vi.foilcalo = std::make_pair(false, true); // foil wire vertices
    negxhits_cside(vi);
  }
  else if (which==3) {
    vi.side = 0; // negative x
    // set up single wirevertex case
    vi.foilcalo = std::make_pair(true, false); // calo wire vertices
    negxhits_fside(vi);
  }

  // need an info vector
  std::vector<VertexInfo> allinfo;
  allinfo.push_back(vi);

  // start the work
  ve.setTrajectory(lf1, allinfo); // runs intersect(0) for a helixfit struct

  printve(ve);
  return (int)ve.onwire().size(); // empty interval or not
}


int check_caseA(){
  return check_wirevertex(0);
}



TEST_CASE( "Case A", "[falaise][wirecheck1]" ) {
  REQUIRE( check_caseA() == 1 );
}

