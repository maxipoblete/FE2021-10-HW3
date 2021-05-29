// Gmsh project created on Wed May  5 23:32:26 2021
SetFactory("OpenCASCADE");

L    = 0.16;
H    = 0.04;
r    = 0.01;
Lext = 0.02;
                
m = 1.2;

Point(1) = {-L/2, -H/2, 0, 1.0};
Point(2) = {L/2,  -H/2, 0, 1.0};
Point(3) = {L/2,   H/2, 0, 1.0};
Point(4) = {-L/2,  H/2, 0, 1.0}; 


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {3, 4, 1, 2};
Plane Surface(1) = {1};


//+
Circle(5) = {0, 0, 0, r, 0, 2*Pi};
Curve Loop(2) = {5};

Plane Surface(2) = {2};
Extrude {-Lext, 0, 0} {
  Curve{4}; }

Extrude {Lext, 0, 0} {
  Curve{2}; }

BooleanDifference{ Surface{1}; Delete; }{Surface{2}; Delete;}



Physical Curve("Empotrado") = {8};
Physical Curve("BordeNatural") = {11};
Physical Surface("Placa") = {1};
Physical Surface("Extremos") = {3, 4};


//+
//+
Transfinite Curve {12, 13} = 16*m Using Progression 1;
//+
Transfinite Curve {4, 8, 2, 11} = 4*m Using Progression 1;
//+
Transfinite Curve {6, 7, 9, 10} = 2*m Using Progression 1;
//+
Transfinite Curve {5} = 2*3.14159*m Using Progression 1;
//+
