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

  std::vector<Rectangle> cvertex = ve.oncalo(); // should be just one entry
  std::cout << "Rectangle on calo: " << std::endl;
  for (auto& cv : cvertex) {
    std::cout << "on calo: 1:[" << cv.axis1.from() << ", " << cv.axis1.to() << "]; 2: [" << cv.axis2.from() << ", " << cv.axis2.to() << "];" << std::endl;
    std::cout << "area fraction: " << cv.areafraction << " on plane " << cv.planeid << " neighbours: (" << cv.neighbourindex.first << ", " << cv.neighbourindex.second << ")" << std::endl;
  }
}


void uphits(VertexInfo& vi) {
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
  mi.column = 8; // a maximum geiger hit
  mi.wirex = 412.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}

void uphits_left(VertexInfo& vi) {
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
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}

void downhits(VertexInfo& vi) {
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
  mi.row = 55; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = 412.0; // not quite but irrelevant here
  mi.wirey = -44.0; // central + 1 down
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}

void downhits_left(VertexInfo& vi) {
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
  mi.row = 55; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // not quite but irrelevant here
  mi.wirey = -44.0; // central + 1 down
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


void xwallhits(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // on side = 1
  mi.row = 111; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = 2438.0; // top right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // on side = 1
  mi.row = 112; // only used for checks at 0 and 112
  mi.column = 4; // a middle geiger hit
  mi.wirex = 184.0; // top right middle
  mi.wirey = 2482.0; // top
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void xwallhits_low(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // on side = 1
  mi.row = 1; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = -2438.0; // bottom right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // on side = 1
  mi.row = 0; // only used for checks at 0 and 112
  mi.column = 4; // a maximum geiger hit
  mi.wirex = 184.0; // middle layer
  mi.wirey = -2482.0; // bottom
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void xwallhits_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 111; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 2438.0; // top right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 112; // only used for checks at 0 and 112
  mi.column = 4; // a middle geiger hit
  mi.wirex = -184.0; // top right middle
  mi.wirey = 2482.0; // top
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void xwallhits_low_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 1; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = -2438.0; // bottom right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 0; // only used for checks at 0 and 112
  mi.column = 4; // a maximum geiger hit
  mi.wirex = -184.0; // middle layer
  mi.wirey = -2482.0; // bottom
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void gvetohits(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // on side = 1
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // centre
  mi.zcoord = 1300.0; // top

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // on side = 1
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 4; // a middle geiger hit
  mi.wirex = 184.0; // right middle
  mi.wirey = 44.0; // bend upwards
  mi.zcoord = 1545.0; // top

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void gvetohits_low(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // on side = 1
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // centre
  mi.zcoord = -1300.0; // bottom

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // on side = 1
  mi.row = 55; // only used for checks at 0 and 112
  mi.column = 4; // a maximum geiger hit
  mi.wirex = 184.0; // middle layer
  mi.wirey = -44.0; // bend downward
  mi.zcoord = -1545.0; // bottom

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void gvetohits_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // top right
  mi.zcoord = 1300.0; // top
  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 4; // a middle geiger hit
  mi.wirex = -184.0; // top right middle
  mi.wirey = 44.0; // bend upwards
  mi.zcoord = 1545.0; // top

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void gvetohits_low_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // bottom right
  mi.zcoord = -1300.0; // bottom

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 55; // only used for checks at 0 and 112
  mi.column = 4; // a maximum geiger hit
  mi.wirex = -184.0; // middle layer
  mi.wirey = -44.0; // bend downwards
  mi.zcoord = -1545.0; // bottom

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void cornerhits_xy(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // positive x
  mi.row = 110; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = 2394.0; // bottom right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // positive x
  mi.row = 112; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = 412.0; // layer 8
  mi.wirey = 2482.0; // bend downwards
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void cornerhits_xy_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // negative x
  mi.row = 110; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 2394.0; // bottom right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // negative x
  mi.row = 112; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // layer 8
  mi.wirey = 2482.0; // bend downwards
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void down_xy(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // positive x
  mi.row = 2; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = -2394.0; // bottom right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // positive x
  mi.row = 0; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = 412.0; // layer 8
  mi.wirey = -2482.0; // bend downwards
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void down_xy_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // negative x
  mi.row = 2; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = -2394.0; // bottom right
  mi.zcoord = 0.0; // central

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // negative x
  mi.row = 0; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // layer 8
  mi.wirey = -2482.0; // bend downwards
  mi.zcoord = 0.0; // central

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void cornerhits_xz_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // central
  mi.zcoord = 1540.0; // top

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = 1540.0; // top

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


void triplehits(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // on side = 0
  mi.row = 110; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = 2394.0; // central
  mi.zcoord = 1540.0; // top

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // on side = 0
  mi.row = 112; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = 412.0; // not quite but irrelevant here
  mi.wirey = 2482.0; // central + 1 up
  mi.zcoord = 1540.0; // top

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


void down_xz_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // on side = 0
  mi.row = 56; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 0.0; // central
  mi.zcoord = -1540.0; // bottom

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // on side = 0
  mi.row = 57; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // not quite but irrelevant here
  mi.wirey = 44.0; // central + 1 up
  mi.zcoord = -1540.0; // bottom

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}


void cornerhits_yz(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 1; // positive x
  mi.row = 110; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = 30.0; // not quite but irrelevant here
  mi.wirey = 2394.0; // bottom right
  mi.zcoord = 1500.0; // top

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 1; // positive x
  mi.row = 112; // only used for checks at 0 and 112
  mi.column = 4; // a maximum geiger hit
  mi.wirex = 184.0; // layer 4
  mi.wirey = 2482.0; // bend downwards
  mi.zcoord = 1548.0; // top

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void cornerhits_yz_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // positive x
  mi.row = 110; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = 2394.0; // bottom right
  mi.zcoord = 1500.0; // top

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // positive x
  mi.row = 112; // only used for checks at 0 and 112
  mi.column = 4; // a maximum geiger hit
  mi.wirex = -184.0; // layer 4
  mi.wirey = 2482.0; // bend downwards
  mi.zcoord = 1548.0; // top

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



void triplehits_left(VertexInfo& vi) {
  // fill the artificial geiger max min information
  MetaInfo mi;
  // min
  mi.hitid = 0;
  mi.side = 0; // positive x
  mi.row = 2; // only used for checks at 0 and 112
  mi.column = 0; // a minimum geiger hit
  mi.wirex = -30.0; // not quite but irrelevant here
  mi.wirey = -2394.0; // bottom right
  mi.zcoord = -1500.0; // top

  vi.minx.push_back(mi);
  vi.miny.push_back(mi);

  // max
  mi.hitid = 1;
  mi.side = 0; // positive x
  mi.row = 0; // only used for checks at 0 and 112
  mi.column = 8; // a maximum geiger hit
  mi.wirex = -412.0; // layer 8
  mi.wirey = -2482.0; // bend downwards
  mi.zcoord = -1548.0; // top

  vi.maxx.push_back(mi);
  vi.maxy.push_back(mi);
}



HelixFit hgveto_up() { // artificial helix fit solution
  HelixFit hf;
  hf.radius = 190.0; // compare to y centre
  hf.pitch = 1000.0; // rising
  hf.xc = -10.0;
  hf.yc = 191.0; 
  hf.zc = 1500.0;
  hf.raderr = 1.0; // small
  hf.errpitch = 10.0;
  hf.errxc = 2.0; // some error in x, y, z centre
  hf.erryc = 2.0;
  hf.errzc = 2.0;
  hf.status = 0;
  hf.clid = 1; // id number
  return hf;
}


HelixFit hgveto_down() { // artificial helix fit solution
  HelixFit hf = hgveto_up();
  hf.zc = -1500.0;
  hf.yc = -191.0; 
  return hf;
}


HelixFit hgveto_up_left() { // artificial helix fit solution
  HelixFit hf = hgveto_up();
  hf.xc = -150.0;
  hf.yc = 191.0;
  hf.zc = 1400.0;
  hf.radius = 190.0; // compare to y centre
  return hf;
}


HelixFit hgveto_down_left() { // artificial helix fit solution
  HelixFit hf = hgveto_up();
  hf.xc = -150.0;
  hf.yc = -191.0; 
  hf.zc = -1360.0;
  hf.radius = 190.0; // compare to y centre
  return hf;
}


HelixFit hxwall_up() { // artificial helix fit solution
  HelixFit hf;
  hf.radius = 500.0; // compare to y centre
  hf.pitch = 0.0; // flat
  hf.xc = -70.0;
  hf.yc = 2930.0; // small helix curve
  hf.zc = 0.0;
  hf.raderr = 1.0; // small
  hf.errpitch = 1.0;  // small
  hf.errxc = 2.0; // some error in x, y, z centre
  hf.erryc = 2.0;
  hf.errzc = 2.0;
  hf.status = 0;
  hf.clid = 1; // id number
  return hf;
}


HelixFit hxwall_down() { // artificial helix fit solution
  HelixFit hf = hxwall_up();
  hf.yc = -2930.0;
  return hf;
}


HelixFit hxwall_up_left() { // artificial helix fit solution
  HelixFit hf = hxwall_up();
  hf.xc = 70.0;
  return hf;
}


HelixFit hxwall_down_left() { // artificial helix fit solution
  HelixFit hf = hxwall_up();
  hf.xc = 70.0;
  hf.yc = -2930.0;
  return hf;
}


HelixFit helixA() { // artificial helix fit solution
  HelixFit hf;
  hf.radius = 5000.0; // compare to y centre
  hf.pitch = 0.0; // flat
  hf.xc = 1.0; // on side = 1
  hf.yc = 5001.0; // central helix curve, y in [1,1001]
  hf.zc = 0.0;
  hf.raderr = 1.0; // small
  hf.errpitch = 1.0;  // small
  hf.errxc = 2.0; // some error in x, y, z centre
  hf.erryc = 2.0;
  hf.errzc = 2.0;
  hf.status = 0;
  hf.clid = 1; // id number
  return hf;
}


HelixFit helixB() { // artificial helix fit solution
  HelixFit hf = helixA();
  hf.xc = -1.0; // side = 0
  return hf;
}


HelixFit helixC() { // artificial helix fit solution
  HelixFit hf = helixA();
  hf.radius = 4000.0; // compare to y centre
  hf.yc = -4001.0; // central helix curve, downward
  return hf;
}


HelixFit helixD() { // artificial helix fit solution
  HelixFit hf = helixA();
  hf.radius = 4000.0; // compare to y centre
  hf.yc = -4001.0; // central helix curve, downward
  hf.xc = -1.0; // side = 0
  return hf;
}



HelixFit hcornerxy() { // corner mainwall xwall
  HelixFit hf;
  hf.radius = 608.0; // compare to y centre
  hf.pitch = 0.0; // flat
  hf.xc = 0.0;
  hf.yc = 2930.0; // corner helix curve
  hf.zc = 0.0;
  hf.raderr = 2.0; // small
  hf.errpitch = 2.0;  // small
  hf.errxc = 10.0; // some error in x, y, z centre
  hf.erryc = 10.0;
  hf.errzc = 10.0;
  hf.status = 0;
  hf.clid = 1; // id number
  return hf;
}


HelixFit hcornerxy_left() { // corner left mainwall, top xwall
  HelixFit hf = hcornerxy();
  return hf;
}


HelixFit hdownxy() { // corner left mainwall, top xwall
  HelixFit hf = hcornerxy();
  hf.yc = -2930.0; // corner helix curve
  return hf;
}


HelixFit hcornerxz() { // corner mainwall gveto
  HelixFit hf;
  hf.radius = 500.0; // compare to y centre
  hf.pitch = 120.0; // tilt
  hf.xc = 0.0;
  hf.yc = 501.0; // corner helix curve
  hf.zc = 1560.0;
  hf.raderr = 2.0; // small
  hf.errpitch = 20.0;  // small
  hf.errxc = 10.0; // some error in x, y, z centre
  hf.erryc = 10.0;
  hf.errzc = 10.0;
  hf.status = 0;
  hf.clid = 1; // id number
  return hf;
}


HelixFit hcornerxz_left() { // corner central mainwall, bottom gveto
  HelixFit hf = hcornerxz();
  hf.zc = 1138.0;
  hf.pitch = 1000.0; // tilt
  return hf;
}


HelixFit hdownxz_left() { // corner central mainwall, bottom gveto
  HelixFit hf = hcornerxz();
  hf.zc = -1138.0;
  hf.pitch = 1000.0; // tilt
  return hf;
}


HelixFit hcorneryz() { // corner xwall gveto
  HelixFit hf;
  hf.radius = 200.0; // compare to y centre
  hf.pitch = 120.0; // flat
  hf.xc = 20.0;
  hf.yc = 2355.0; // corner helix curve
  hf.zc = 1535.0;
  hf.raderr = 2.0; // small
  hf.errpitch = 20.0;  // small
  hf.errxc = 10.0; // some error in x, y, z centre
  hf.erryc = 10.0;
  hf.errzc = 10.0;
  hf.status = 0;
  hf.clid = 1; // id number
  return hf;
}


HelixFit hcorneryz_left() { // corner xwall gveto
  HelixFit hf = hcorneryz();
  hf.radius = 200.0; // compare to y centre
  hf.xc = -500.0;
  hf.yc = 2355.0;
  hf.zc = 1450;
  hf.pitch = 750.0;
  return hf;
}


HelixFit triple() { // corner all surfaces
  HelixFit hf = hcornerxz();
  hf.xc = 5.0;
  hf.yc = 2250.0;
  // hf.zc = 1138.0;
  // hf.pitch = 1000.0; // tilt
  return hf;
}


HelixFit triple_left() { // corner all surfaces
  HelixFit hf = hcornerxz();
  hf.xc = -10.0;
  hf.yc = -2760.0;
  hf.zc = -1138.0;
  hf.pitch = 1000.0; // tilt
  return hf;
}




int check_helix(HelixFit hf, int which) {
  VertexExtrapolator ve(make_planes());
  VertexInfo vi;
  vi.side = 1; // positive x
  vi.foilcalo = std::make_pair(true, true); // no wire vertices
  vi.clsid = 1; // same clsid as hook
  if (which==0)
    uphits(vi);
  else if (which==1)
    downhits(vi);
  else if (which==2)
    xwallhits(vi);
  else if (which==3)
    xwallhits_low(vi);
  else if (which==4)
    gvetohits(vi);
  else if (which==5)
    gvetohits_low(vi);
  else if (which==6)
    cornerhits_xy(vi);
  else if (which==7)
    cornerhits_yz(vi);
  else if (which==8)
    down_xy(vi);
  else if (which==9)
    triplehits(vi);

  // need an info vector
  std::vector<VertexInfo> allinfo;
  allinfo.push_back(vi);
  // start the work
  ve.setTrajectory(hf, allinfo); // runs intersect(0) for a helixfit struct

  printve(ve);
  return (int)ve.oncalo().size();
}

int check_helix_left(HelixFit hf, int which) {
  VertexExtrapolator ve(make_planes());
  VertexInfo vi;
  vi.side = 0; // negative x
  vi.foilcalo = std::make_pair(true, true); // no wire vertices
  vi.clsid = 1; // same clsid as hook
  if (which==0)
    uphits_left(vi);
  else if (which==1)
    downhits_left(vi);
  else if (which==2)
    xwallhits_left(vi);
  else if (which==3)
    xwallhits_low_left(vi);
  else if (which==4)
    gvetohits_left(vi);
  else if (which==5)
    gvetohits_low_left(vi);
  else if (which==6)
    cornerhits_xy_left(vi);
  else if (which==7)
    cornerhits_xz_left(vi);
  else if (which==8)
    down_xz_left(vi);
  else if (which==9)
    down_xy_left(vi);
  else if (which==10)
    cornerhits_yz_left(vi);
  else if (which==11)
    triplehits_left(vi);

  // need an info vector
  std::vector<VertexInfo> allinfo;
  allinfo.push_back(vi);
  // start the work
  ve.setTrajectory(hf, allinfo); // runs intersect(0) for a helixfit struct

  printve(ve);
  return (int)ve.oncalo().size();
}


int check_helixA(){
  HelixFit hf = helixA(); // set vertex info
  return check_helix(hf,0);
}


int check_helixB(){
  HelixFit hf = helixB(); // set vertex info
  return check_helix_left(hf,0);
}


int check_helixC(){
  HelixFit hf = helixC(); // set vertex info
  return check_helix(hf,1);
}


int check_helixD(){
  HelixFit hf = helixD(); // set vertex info
  return check_helix_left(hf,1);
}


int check_xwall_up(){
  HelixFit hf = hxwall_up(); // set vertex info
  return check_helix(hf,2);
}


int check_xwall_down(){
  HelixFit hf = hxwall_down(); // set vertex info
  return check_helix(hf,3);
}


int left_xwall_up(){
  HelixFit hf = hxwall_up_left(); // set vertex info
  return check_helix_left(hf,2);
}


int left_xwall_down(){
  HelixFit hf = hxwall_down_left(); // set vertex info
  return check_helix_left(hf,3);
}


int check_gveto_up(){
  HelixFit hf = hgveto_up(); // set vertex info
  return check_helix(hf,4);
}


int check_gveto_down(){
  HelixFit hf = hgveto_down(); // set vertex info
  return check_helix(hf,5);
}


int left_gveto_up(){
  HelixFit hf = hgveto_up_left(); // set vertex info
  return check_helix_left(hf,4);
}


int left_gveto_down(){
  HelixFit hf = hgveto_down_left(); // set vertex info
  return check_helix_left(hf,5);
}


int check_corner_xy(){
  HelixFit hf = hcornerxy(); // set vertex info
  return check_helix(hf,6);
}


int check_down_xy(){
  HelixFit hf = hdownxy(); // set vertex info
  return check_helix(hf,8);
}


int left_down_xy(){
  HelixFit hf = hdownxy(); // set vertex info
  return check_helix_left(hf,9);
}


int left_corner_xy(){
  HelixFit hf = hcornerxy_left(); // set vertex info
  return check_helix_left(hf,6);
}


int check_corner_xz(){
  HelixFit hf = hcornerxz(); // set vertex info
  return check_helix(hf,0);
}


int left_corner_xz(){
  HelixFit hf = hcornerxz_left(); // set vertex info
  return check_helix_left(hf,7);
}


int left_down_xz(){
  HelixFit hf = hdownxz_left(); // set vertex info
  return check_helix_left(hf,8);
}


int check_corner_yz(){
  HelixFit hf = hcorneryz(); // set vertex info
  return check_helix(hf,7);
}


int left_corner_yz(){
  HelixFit hf = hcorneryz_left(); // set vertex info
  return check_helix_left(hf,10);
}


int left_xyz(){
  HelixFit hf = triple_left(); // set vertex info
  return check_helix_left(hf,11);
}


int xyz(){
  HelixFit hf = triple(); // set vertex info
  return check_helix(hf,9);
}




TEST_CASE( "Helix A", "[falaise][helixcheck1]" ) {
  REQUIRE( check_helixA() == 1 );
}

TEST_CASE( "Helix B", "[falaise][helixcheck2]" ) {
  REQUIRE( check_helixB() == 1 );
}

TEST_CASE( "Helix C", "[falaise][helixcheck3]" ) {
  REQUIRE( check_helixC() == 1 );
}

TEST_CASE( "Helix D", "[falaise][helixcheck4]" ) {
  REQUIRE( check_helixD() == 1 );
}

TEST_CASE( "XWall up", "[falaise][helixcheck5]" ) {
  REQUIRE( check_xwall_up() == 1 );
}

TEST_CASE( "XWall down", "[falaise][helixcheck6]" ) {
  REQUIRE( check_xwall_down() == 1 );
}

TEST_CASE( "left XWall up", "[falaise][helixcheck7]" ) {
  REQUIRE( left_xwall_up() == 1 );
}

TEST_CASE( "left XWall down", "[falaise][helixcheck8]" ) {
  REQUIRE( left_xwall_down() == 1 );
}

TEST_CASE( "gveto up", "[falaise][helixcheck9]" ) {
  REQUIRE( check_gveto_up() == 1 );
}

TEST_CASE( "gveto down", "[falaise][helixcheck10]" ) {
  REQUIRE( check_gveto_down() == 1 );
}

TEST_CASE( "left gveto up", "[falaise][helixcheck11]" ) {
  REQUIRE( left_gveto_up() == 1 );
}

TEST_CASE( "left gveto down", "[falaise][helixcheck12]" ) {
  REQUIRE( left_gveto_down() == 1 );
}

TEST_CASE( "corner xy", "[falaise][helixcheck13]" ) {
  REQUIRE( check_corner_xy() == 2 );
}

TEST_CASE( "corner xz", "[falaise][helixcheck14]" ) {
  REQUIRE( check_corner_xz() == 2 );
}

TEST_CASE( "corner yz", "[falaise][helixcheck15]" ) {
  REQUIRE( check_corner_yz() == 2 );
}

TEST_CASE( "left corner xy", "[falaise][helixcheck16]" ) {
  REQUIRE( left_corner_xy() == 2 );
}

TEST_CASE( "left corner xz", "[falaise][helixcheck17]" ) {
  REQUIRE( left_corner_xz() == 2 );
}

TEST_CASE( "left down xz", "[falaise][helixcheck18]" ) {
  REQUIRE( left_down_xz() == 2 );
}

TEST_CASE( "down xy", "[falaise][helixcheck19]" ) {
  REQUIRE( check_down_xy() == 2 );
}

TEST_CASE( "left down xy", "[falaise][helixcheck20]" ) {
  REQUIRE( left_down_xy() == 2 );
}

TEST_CASE( "left corner yz", "[falaise][helixcheck21]" ) {
  REQUIRE( left_corner_yz() == 2 );
}

TEST_CASE( "left corner xyz", "[falaise][helixcheck22]" ) {
  REQUIRE( left_xyz() == 3 );
}

TEST_CASE( "corner xyz", "[falaise][helixcheck23]" ) {
  REQUIRE( xyz() == 3 );
}

