// Gmsh project created on Wed May  5 23:32:26 2021
SetFactory("OpenCASCADE");

L = 16.0;
H = 4.0;
r = 1.0;
Lextremo =2.0;
tplaca = 0.4;
textremo = 0.5;
m = 0.4;
mc = 50;

Point(1) = {-L/2, -H/2, 0, m};
Point(2) = {L/2,  -H/2, 0, m};
Point(3) = {L/2,   H/2, 0, m};
Point(4) = {-L/2,  H/2, 0, m};


//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};

Extrude {-Lextremo,0,0} {
	Curve{4};
}

Extrude {Lextremo,0,0} {
	Curve{2};
}


//+
Circle(11) = {0, 0, 0, r, 0, 2*Pi};
//+
Curve Loop(4) = {11};
//+
Plane Surface(4) = {4};


BooleanDifference{ Surface{1}; Delete; }{Surface{4}; Delete;}


//+
Physical Curve("Fixed") = {7};


//+
Physical Curve("Prescribed Traction Stress") = {10};
//+


Transfinite Curve {11} = mc Using Progression 1;
//+

Extrude {0, 0, tplaca/2} {
  Surface{1}; Layers {5}; Recombine;
}
Extrude {0, 0, -tplaca/2} {
  Surface{1}; Layers {5}; Recombine;
}
//+
Extrude {0, 0, textremo/2} {
  Surface{2}; Surface{3}; Layers {5}; Recombine;
}
Extrude {0, 0, -textremo/2} {
  Surface{2}; Surface{3}; Layers {5}; Recombine;
}


