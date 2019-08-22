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


LineFit lineB() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = 10.0;
  return lf;
}


LineFit lineC() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = 0.0;
  lf.slxz = 10.0; // up in z
  return lf;
}


LineFit lineD() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = -1.0;
  return lf;
}


LineFit lineE() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = -10.0;
  return lf;
}


LineFit lineF() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = 0.0;
  lf.slxz = -10.0;
  return lf;
}


LineFit borderxwall() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = 5.75;
  return lf;
}


LineFit bordermwall() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = 5.9;
  return lf;
}


LineFit borderxwall_low() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = -5.75;
  return lf;
}


LineFit bordermwall_low() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxy = -5.9;
  return lf;
}


LineFit bordergveto() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxz = 3.6; // diagonal up
  return lf;
}



LineFit bordermwallgveto() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxz = 3.5; // diagonal up
  return lf;
}



LineFit bordergvetolow() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxz = -3.6; // diagonal down
  return lf;
}



LineFit bordermwallgvetolow() { // artificial line fit solution
  LineFit lf = lineA();
  lf.slxz = -3.5; // diagonal down
  return lf;
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

int check_line(LineFit lf) {
  VertexExtrapolator ve(make_planes());
  VertexInfo vi;
  vi.side = 1; // positive x
  vi.foilcalo = std::make_pair(true, true); // no wire vertices
  vi.clsid = 1; // same clsid as hook
  // no other settings required on vi for a line

  // need an info vector
  std::vector<VertexInfo> allinfo;
  allinfo.push_back(vi);
  // start the work
  ve.setTrajectory(lf, allinfo); // runs intersect(0) for a linefit struct

  printve(ve);
  return (int)ve.oncalo().size();
}

int check_lineminus(LineFit lf) {
  VertexExtrapolator ve(make_planes());
  VertexInfo vi;
  vi.side = 0; // negative x
  vi.foilcalo = std::make_pair(true, true); // no wire vertices
  vi.clsid = 1; // same clsid as hook
  // no other settings required on vi for a line

  // need an info vector
  std::vector<VertexInfo> allinfo;
  allinfo.push_back(vi);

  // start the work
  ve.setTrajectory(lf, allinfo); // runs intersect(0) for a linefit struct

  printve(ve);
  return (int)ve.oncalo().size();
}


int check_lineA(){
  LineFit lf = lineA(); // set vertex info
  return check_line(lf);
}


int check_lineAminus(){
  LineFit lf = lineA(); // set vertex info
  return check_lineminus(lf);
}


int check_lineB(){
  LineFit lf = lineB(); // set vertex info
  return check_line(lf);
}

int check_lineBminus(){
  LineFit lf = lineB(); // set vertex info
  return check_lineminus(lf);
}


int check_lineC(){
  LineFit lf = lineC(); // set vertex info
  return check_line(lf);
}

int check_lineCminus(){
  LineFit lf = lineC(); // set vertex info
  return check_lineminus(lf);
}


int check_lineD(){
  LineFit lf = lineD(); // set vertex info
  return check_line(lf);
}

int check_lineDminus(){
  LineFit lf = lineD(); // set vertex info
  return check_lineminus(lf);
}


int check_lineE(){
  LineFit lf = lineE(); // set vertex info
  return check_line(lf);
}

int check_lineEminus(){
  LineFit lf = lineE(); // set vertex info
  return check_lineminus(lf);
}


int check_lineF(){
  LineFit lf = lineF(); // set vertex info
  return check_line(lf);
}

int check_lineFminus(){
  LineFit lf = lineF(); // set vertex info
  return check_lineminus(lf);
}


int check_bxwall(){
  LineFit lf = borderxwall(); // set vertex info
  return check_line(lf);
}


int check_bmwall(){
  LineFit lf = bordermwall(); // set vertex info
  return check_line(lf);
}


int check_bxwall_low(){
  LineFit lf = borderxwall_low(); // set vertex info
  return check_line(lf);
}


int check_bmwall_low(){
  LineFit lf = bordermwall_low(); // set vertex info
  return check_line(lf);
}


int minus_bxwall(){
  LineFit lf = borderxwall(); // set vertex info
  return check_lineminus(lf);
}


int minus_bmwall(){
  LineFit lf = bordermwall(); // set vertex info
  return check_lineminus(lf);
}


int minus_bxwall_low(){
  LineFit lf = borderxwall_low(); // set vertex info
  return check_lineminus(lf);
}


int minus_bmwall_low(){
  LineFit lf = bordermwall_low(); // set vertex info
  return check_lineminus(lf);
}


int check_bgveto(){
  LineFit lf = bordergveto(); // set vertex info
  return check_line(lf);
}


int check_bmg(){
  LineFit lf = bordermwallgveto(); // set vertex info
  return check_line(lf);
}


int check_bgveto_low(){
  LineFit lf = bordergvetolow(); // set vertex info
  return check_line(lf);
}


int check_bmg_low(){
  LineFit lf = bordermwallgvetolow(); // set vertex info
  return check_line(lf);
}


int minus_bgveto(){
  LineFit lf = bordergveto(); // set vertex info
  return check_lineminus(lf);
}


int minus_bmg(){
  LineFit lf = bordermwallgveto(); // set vertex info
  return check_lineminus(lf);
}


int minus_bgveto_low(){
  LineFit lf = bordergvetolow(); // set vertex info
  return check_lineminus(lf);
}


int minus_bmg_low(){
  LineFit lf = bordermwallgvetolow(); // set vertex info
  return check_lineminus(lf);
}



TEST_CASE( "Line A", "[falaise][planecheck1]" ) {
  REQUIRE( check_lineA() == 1 );
}

TEST_CASE( "Line A-", "[falaise][planecheck0]" ) {
  REQUIRE( check_lineAminus() == 1 );
}

TEST_CASE( "Line B", "[falaise][planecheck2]" ) {
  REQUIRE( check_lineB() == 1 );
}

TEST_CASE( "Line B-", "[falaise][planecheck3]" ) {
  REQUIRE( check_lineBminus() == 1 );
}

TEST_CASE( "Line C", "[falaise][planecheck4]" ) {
  REQUIRE( check_lineC() == 1 );
}

TEST_CASE( "Line C-", "[falaise][planecheck5]" ) {
  REQUIRE( check_lineCminus() == 1 );
}

TEST_CASE( "Line D", "[falaise][planecheck6]" ) {
  REQUIRE( check_lineD() == 1 );
}

TEST_CASE( "Line D-", "[falaise][planecheck7]" ) {
  REQUIRE( check_lineDminus() == 1 );
}

TEST_CASE( "Line E", "[falaise][planecheck8]" ) {
  REQUIRE( check_lineE() == 1 );
}

TEST_CASE( "Line E-", "[falaise][planecheck9]" ) {
  REQUIRE( check_lineEminus() == 1 );
}

TEST_CASE( "Line F", "[falaise][planecheck10]" ) {
  REQUIRE( check_lineF() == 1 );
}

TEST_CASE( "Line F-", "[falaise][planecheck11]" ) {
  REQUIRE( check_lineFminus() == 1 );
}

TEST_CASE( "Border xwall", "[falaise][planecheck12]" ) {
  REQUIRE( check_bxwall() == 2 );
}

TEST_CASE( "Border mainwall", "[falaise][planecheck13]" ) {
  REQUIRE( check_bmwall() == 2 );
}

TEST_CASE( "Border xwall low", "[falaise][planecheck14]" ) {
  REQUIRE( check_bxwall_low() == 2 );
}

TEST_CASE( "Border mainwall low", "[falaise][planecheck15]" ) {
  REQUIRE( check_bmwall_low() == 2 );
}

TEST_CASE( "Border neg xwall", "[falaise][planecheck16]" ) {
  REQUIRE( minus_bxwall() == 2 );
}

TEST_CASE( "Border neg mainwall", "[falaise][planecheck17]" ) {
  REQUIRE( minus_bmwall() == 2 );
}

TEST_CASE( "Border neg xwall low", "[falaise][planecheck18]" ) {
  REQUIRE( minus_bxwall_low() == 2 );
}

TEST_CASE( "Border neg mainwall low", "[falaise][planecheck19]" ) {
  REQUIRE( minus_bmwall_low() == 2 );
}

TEST_CASE( "Border gveto", "[falaise][planecheck20]" ) {
  REQUIRE( check_bgveto() == 2 );
}

TEST_CASE( "Border mwall gveto", "[falaise][planecheck21]" ) {
  REQUIRE( check_bmg() == 2 );
}

TEST_CASE( "Border gveto low", "[falaise][planecheck22]" ) {
  REQUIRE( check_bgveto_low() == 2 );
}

TEST_CASE( "Border mwall gveto low", "[falaise][planecheck23]" ) {
  REQUIRE( check_bmg_low() == 2 );
}

TEST_CASE( "Border neg gveto", "[falaise][planecheck24]" ) {
  REQUIRE( minus_bgveto() == 2 );
}

TEST_CASE( "Border neg mwall gveto", "[falaise][planecheck25]" ) {
  REQUIRE( minus_bmg() == 2 );
}

TEST_CASE( "Border neg gveto low", "[falaise][planecheck26]" ) {
  REQUIRE( minus_bgveto_low() == 2 );
}

TEST_CASE( "Border neg mwall gveto low", "[falaise][planecheck27]" ) {
  REQUIRE( minus_bmg_low() == 2 );
}

