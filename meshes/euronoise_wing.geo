Mesh.MshFileVersion = 2.2; // choose the correct mesh format

lc=1.5; // caracteristic size of elements
width=105; // physical width of rectangular domain (without PMLs)
height=15; // physical height of rectangular domain (without PMLs)
xleft=0; // x position of left side of computational domain
PMLxb_N = 20; // number of PML elements
PMLxt_N = 20;
PMLyb_N = 0;
PMLyt_N = 32;

// canyon parameters
L = 15; // length of buidlings
H = 12; // height of buidlings
sepRatio = 0.25; // separation ration (H/W)
x1 = 15; // x position of left side of left building
W = H/sepRatio; // canyon width
W = 45;

// physical groups for boundaries
BC_rigid = 1; // rigid / perfectly reflecting BC
BC_impedance = 2; // custom BC
// everything else == nonreflecting BC

// define points
pxbyb = newp; Point(pxbyb) = {xleft, 0, 0, lc}; // bottom left
pxtyb = newp; Point(pxtyb) = {xleft+width, 0, 0, lc}; // bottom right
pxtyt = newp; Point(pxtyt) = {xleft+width, height, 0, lc}; // top right
pxbyt = newp; Point(pxbyt) = {xleft, height, 0, lc}; // top left

// get coordinates, to build PMLs
cxbyb[] = Point{pxbyb};
cxtyb[] = Point{pxtyb};
cxtyt[] = Point{pxtyt};
cxbyt[] = Point{pxbyt};

// define lines from the points
//lbot   = newl; Line(lbot)   = {pxbyb, pxtyb}; // south
lright = newl; Line(lright) = {pxtyb, pxtyt}; // east
ltop   = newl; Line(ltop)   = {pxtyt, pxbyt}; // north
lleft  = newl; Line(lleft)  = {pxbyt, pxbyb}; // west

// manually define the discretization of the outer boundaries
Transfinite Curve {lleft, lright} = height/lc+1 Using Progression 1;
Transfinite Curve {ltop} = width/lc+1 Using Progression 1;
// ground will be defined later

// add PMLs
If (PMLxb_N>0)
  spmlxb[] = Extrude { -PMLxb_N*lc, 0, 0 }{Line{lleft}; Layers{PMLxb_N}; Recombine;};
EndIf
If (PMLxt_N>0)
  spmlxt[] = Extrude { PMLxt_N*lc, 0, 0 }{Line{lright}; Layers{PMLxt_N}; Recombine;};
EndIf
If (PMLyb_N>0)
  spmlyb[] = Extrude { 0, -PMLyb_N*lc, 0 }{Line{lbot}; Layers{PMLyb_N}; Recombine;};
EndIf
If (PMLyt_N>0)
  spmlyt[] = Extrude { 0, PMLyt_N*lc, 0 }{Line{ltop}; Layers{PMLyt_N}; Recombine;};
EndIf
If (PMLxb_N>0 && PMLyt_N>0)
  spmlxbyt[] = Extrude { -PMLxb_N*lc, 0, 0 }{Line{spmlyt[2]}; Layers{PMLxb_N}; Recombine;};
EndIf
If (PMLxt_N>0 && PMLyt_N>0)
  spmlxtyt[] = Extrude { PMLxt_N*lc, 0, 0 }{Line{spmlyt[3]}; Layers{PMLxt_N}; Recombine;};
EndIf
If (PMLxb_N>0 && PMLyb_N>0)
  spmlxbyb[] = Extrude { -PMLxb_N*lc, 0, 0 }{Line{spmlyb[3]}; Layers{PMLxb_N}; Recombine;};
EndIf
If (PMLxt_N>0 && PMLyb_N>0)
  spmlxtyb[] = Extrude { PMLxt_N*lc, 0, 0 }{Line{spmlyb[2]}; Layers{PMLxt_N}; Recombine;};
EndIf

bthick=0.3; // thickness of balconies
blen = 1; // length of balconies
bsep = 3; // separation

// ground/urban geometry
xshift = x1;
p1 = newp; Point(p1) = {xshift+0, 0, 0, lc};
pp = newp; Point(pp) = {xshift+0, H, 0, lc};
pp = newp; Point(pp) = {xshift+L, H, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen, H, 0, lc}; // balc
pp = newp; Point(pp) = {xshift+L+blen, H-bthick, 0, lc};
pp = newp; Point(pp) = {xshift+L, H-bthick, 0, lc}; // end
pp = newp; Point(pp) = {xshift+L, H-bsep+bthick, 0, lc}; // balc
pp = newp; Point(pp) = {xshift+L+blen-bthick, H-bsep+bthick, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen-bthick, H-bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen, H-bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen, H-bsep, 0, lc};
pp = newp; Point(pp) = {xshift+L, H-bsep, 0, lc}; // end
pp = newp; Point(pp) = {xshift+L, H-2*bsep+bthick, 0, lc}; // balc
pp = newp; Point(pp) = {xshift+L+blen-bthick, H-2*bsep+bthick, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen-bthick, H-2*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen, H-2*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen, H-2*bsep, 0, lc};
pp = newp; Point(pp) = {xshift+L, H-2*bsep, 0, lc}; // end
pp = newp; Point(pp) = {xshift+L, H-3*bsep+bthick, 0, lc}; // balc
pp = newp; Point(pp) = {xshift+L+blen-bthick, H-3*bsep+bthick, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen-bthick, H-3*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen, H-3*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+L+blen, H-3*bsep, 0, lc};
pp = newp; Point(pp) = {xshift+L, H-3*bsep, 0, lc}; // end
//pp = newp; Point(pp) = {xshift+L+blen-bthick, H-bsep, 0, lc};
p2 = newp; Point(p2) = {xshift+L, 0, 0, lc};
pgrd1 = {p1:p2}; // list of ground points
Npgrd1 = #pgrd1[]; // number of points

xshift = xshift+L+W;
p3 = newp; Point(p3) = {xshift+0, 0, 0, lc};
pp = newp; Point(pp) = {xshift+0, bsep, 0, lc}; // balc
pp = newp; Point(pp) = {xshift+0-blen, bsep, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen, bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen+bthick, bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen+bthick, bsep+bthick, 0, lc};
pp = newp; Point(pp) = {xshift+0, bsep+bthick, 0, lc}; // end
pp = newp; Point(pp) = {xshift+0, 2*bsep, 0, lc}; // balc
pp = newp; Point(pp) = {xshift+0-blen, 2*bsep, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen, 2*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen+bthick, 2*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen+bthick, 2*bsep+bthick, 0, lc};
pp = newp; Point(pp) = {xshift+0, 2*bsep+bthick, 0, lc}; // end
pp = newp; Point(pp) = {xshift+0, 3*bsep, 0, lc}; // balc
pp = newp; Point(pp) = {xshift+0-blen, 3*bsep, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen, 3*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen+bthick, 3*bsep+blen, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen+bthick, 3*bsep+bthick, 0, lc};
pp = newp; Point(pp) = {xshift+0, 3*bsep+bthick, 0, lc}; // end
pp = newp; Point(pp) = {xshift+0, H-bthick, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen, H-bthick, 0, lc};
pp = newp; Point(pp) = {xshift+0-blen, H, 0, lc}; // end
pp = newp; Point(pp) = {xshift+0, H, 0, lc};
pp = newp; Point(pp) = {xshift+L, H, 0, lc};
p4 = newp; Point(p4) = {xshift+L, 0, 0, lc};
pgrd2 = {p3:p4}; // list of ground points
Npgrd2 = #pgrd2[]; // number of points

// create lines between the points
l1 = newl; Line(l1) = {pxbyb,pgrd1[0]};
For ii In {0:Npgrd1-2}
  ll = newl; Line(ll) = {pgrd1[ii],pgrd1[ii+1]};
EndFor
l2 = newl; Line(l2) = {pgrd1[Npgrd1-1],pgrd2[0]};
For ii In {0:Npgrd2-2}
  ll = newl; Line(ll) = {pgrd2[ii],pgrd2[ii+1]};
EndFor
l3 = newl; Line(l3) = {pgrd2[Npgrd2-1],pxtyb};

// define closed path from lines
llmain = newll; Curve Loop(llmain) = {l1:l3,lright,ltop,lleft};
//llmain = newll; Curve Loop(llmain) = {lbot,lright,ltop,lleft};

// define surface from closed path
smain = news; Plane Surface(smain) = {llmain};

Recombine Surface {smain}; // get quads


/* Define physical groups, for boundary conditions. If you want to add
physical groups, all mesh entites (e.g., lines and surfaces) need to
be assigned a group, otherwise they won't be exported */

//Physical Curve(BC_impedance) = {lbot};
//Physical Curve(BC_impedance) = {spmlxb[2], l1, l2, l3, spmlxt[3]};
Physical Curve(BC_rigid) = {spmlxb[2], l1, l2, l3, spmlxt[3],l1+1:l2-1,l2+1:l3-1};

// for surfaces, you can put any group since we won't use it
Physical Surface(0) = {smain};
If (PMLxb_N>0)
  Physical Surface(100) = {spmlxb[1]};
EndIf
If (PMLxt_N>0)
  Physical Surface(200) = {spmlxt[1]};
EndIf
If (PMLyb_N>0)
  Physical Surface(300) = {spmlyb[1]};
EndIf
If (PMLyt_N>0)
  Physical Surface(400) = {spmlyt[1]};
EndIf
If (PMLxb_N>0 && PMLyt_N>0)
  Physical Surface(500) = {spmlxbyt[1]};
EndIf
If (PMLxt_N>0 && PMLyt_N>0)
  Physical Surface(600) = {spmlxtyt[1]};
EndIf
If (PMLxb_N>0 && PMLyb_N>0)
  Physical Surface(700) = {spmlxbyb[1]};
EndIf
If (PMLxt_N>0 && PMLyb_N>0)
  Physical Surface(800) = {spmlxtyb[1]};
EndIf

Coherence; // remove entities sharing the same coordinates

/*
Physical Curve(0) = {3, 4};
Physical Curve(1) = {2};
Physical Curve(7) = {1};
Physical Surface(smain) = {1};
*/

// Note: you can generate the mesh from the command line with:
//         >> gmsh -2 simple_example.geo

// to display something in the console
//Printf("%f",#pgrd[]);
