
lc = 0.07;
r  = 0.1;
R  = 1.0;
H  = 1.0;

Point(1) = { 0, 0, 0,lc};
Point(2) = { 0, 0, H,lc};
Point(3) = { r, 0, H,lc};

cakePts[0] = newp; Point(cakePts[0]) = {R,      0, H,      lc};
cakePts[1] = newp; Point(cakePts[1]) = {R,      0, H*0.8,  lc};
cakePts[2] = newp; Point(cakePts[2]) = {R*0.75, 0, H*0.35, lc};
cakePts[3] = newp; Point(cakePts[3]) = {R*0.2,  0, 0,      lc};
cakePts[4] = newp; Point(cakePts[4]) = {r,      0, 0,      lc};

BSpline(1) = {cakePts[]};

Line(2) = {cakePts[4],1};
Line(3) = {1,2};
Line(4) = {2,3};
Line(5) = {3,cakePts[0]};
Line(6) = {3,cakePts[4]};

Line Loop(10) = {2,3,4,6};
Line Loop(11) = {1,-6,5};

Plane Surface(20) = {10};
Plane Surface(30) = {11};

inVol1[]  = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Surface{20}; };
inVol2[]  = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Surface{inVol1[0]}; };
inVol3[]  = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Surface{inVol2[0]}; };

outVol1[] = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Surface{30}; };
outVol2[] = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Surface{outVol1[0]}; };
outVol3[] = Extrude{ {0,0,1}, {0,0,0}, 2*Pi/3 } { Surface{outVol2[0]}; };

Physical Volume(1000) = {outVol1[1], outVol2[1], outVol3[1]};
Physical Volume(1001) = {inVol1[1],  inVol2[1],  inVol3[1] };

// hard coded these surfaces... enumeration found from the GUI.
// will change if you add or remove stuff above
Physical Surface(2000) = {88, 105, 96, 113, 122 , 130, 38, 55 , 72, 59, 76, 42};
