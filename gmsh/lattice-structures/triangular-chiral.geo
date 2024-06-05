// Gmsh project created on Mon Nov 13 13:27:42 2023
SetFactory("OpenCASCADE");
Geometry.CopyMeshingMethod = 1;

L = 15; // mm
w = 0.1; // ratio of thickness to beam length

// Make hexagonal structure with triangles

// Outer hexagon curve
L_out = 2*L;
Point(1) = {L_out*Cos(0), L_out*Sin(0), 0, 0};
Point(2) = {L_out*Cos(Pi/3), L_out*Sin(Pi/3), 0, 0};
Point(3) = {L_out*Cos(2*Pi/3), L_out*Sin(2*Pi/3), 0, 0};
Point(4) = {L_out*Cos(Pi), L_out*Sin(Pi), 0, 0};
Point(5) = {L_out*Cos(4*Pi/3), L_out*Sin(4*Pi/3), 0, 0};
Point(6) = {L_out*Cos(5*Pi/3), L_out*Sin(5*Pi/3), 0, 0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Curve Loop(1) = {1:6};

// Inner hexagon curve
l = L_out-w*L/Sin(Pi/3);
Point(7) = {l*Cos(0), l*Sin(0), 0, 0};
Point(8) = {l*Cos(Pi/3), l*Sin(Pi/3), 0, 0};
Point(9) = {l*Cos(2*Pi/3), l*Sin(2*Pi/3), 0, 0};
Point(10) = {l*Cos(Pi), l*Sin(Pi), 0, 0};
Point(11) = {l*Cos(4*Pi/3), l*Sin(4*Pi/3), 0, 0};
Point(12) = {l*Cos(5*Pi/3), l*Sin(5*Pi/3), 0, 0};

Line(7) = {7,8};
Line(8) = {8,9};
Line(9) = {9,10};
Line(10) = {10,11};
Line(11) = {11,12};
Line(12) = {12,7};

Curve Loop(2) = {7:12};

Plane Surface(1) = {1,2};

// NOTE: Duplicata automatically assings a tag. Get it from the GUI. Do not assign a tag to avoid errors!

// Strips
Rectangle(2) = {-0.95*L_out, -0.5*w*L, 0, 1.9*L_out, w*L, 0};
Rotate {{0, 0, 1}, {0, 0, 0}, Pi/3} { Duplicata{Surface{2};} }
Rotate {{0, 0, 1}, {0, 0, 0}, -Pi/3} { Duplicata{Surface{2};} }

// Elastic structure
BooleanUnion(5) = { Surface{1}; Delete;}{ Surface{2:4}; Delete;};

//// Magnet slots

// Center
Rectangle(6) = {-1.25*w*L, -1.25*w*L, 0, 2.5*w*L, 2.5*w*L, 0};

// Inner magnets
lm = L-0.5*w*L;
// Surfaces 7-12
Translate {lm*Cos(0), lm*Sin(0), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(Pi/3), lm*Sin(Pi/3), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(2*Pi/3), lm*Sin(2*Pi/3), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(Pi), lm*Sin(Pi), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(4*Pi/3), lm*Sin(4*Pi/3), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(5*Pi/3), lm*Sin(5*Pi/3), 0} { Duplicata{Surface{6};} }

// Outer magnets along the vertices
lm = L_out-w*L;
// Surfaces 13-18
Translate {lm*Cos(0), lm*Sin(0), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(Pi/3), lm*Sin(Pi/3), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(2*Pi/3), lm*Sin(2*Pi/3), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(Pi), lm*Sin(Pi), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(4*Pi/3), lm*Sin(4*Pi/3), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(5*Pi/3), lm*Sin(5*Pi/3), 0} { Duplicata{Surface{6};} }

// Outer magnets along the edges
lm = (L_out-w*L)*Cos(Pi/6);
// Surfaces 19-24
Translate {lm*Cos(Pi/6), lm*Sin(Pi/6), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(3*Pi/2), lm*Sin(3*Pi/6), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(5*Pi/6), lm*Sin(5*Pi/6), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(7*Pi/6), lm*Sin(7*Pi/6), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(9*Pi/6), lm*Sin(9*Pi/6), 0} { Duplicata{Surface{6};} }
Translate {lm*Cos(11*Pi/6), lm*Sin(11*Pi/6), 0} { Duplicata{Surface{6};} }

BooleanUnion(25) = { Surface{5}; Delete;}{ Surface{6:24}; Delete;};

// Add magnets
Rectangle(26) = {-2/3*w*L, -2/3*w*L, 0, 4/3*w*L, 4/3*w*L, 0};

// Inner magnets
lm = L-0.5*w*L;
// Surfaces 27-32
Translate {lm*Cos(0), lm*Sin(0), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(Pi/3), lm*Sin(Pi/3), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(2*Pi/3), lm*Sin(2*Pi/3), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(Pi), lm*Sin(Pi), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(4*Pi/3), lm*Sin(4*Pi/3), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(5*Pi/3), lm*Sin(5*Pi/3), 0} { Duplicata{Surface{26};} }

// Outer magnets along the vertices
lm = L_out-w*L;
// Surfaces 33-38
Translate {lm*Cos(0), lm*Sin(0), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(Pi/3), lm*Sin(Pi/3), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(2*Pi/3), lm*Sin(2*Pi/3), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(Pi), lm*Sin(Pi), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(4*Pi/3), lm*Sin(4*Pi/3), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(5*Pi/3), lm*Sin(5*Pi/3), 0} { Duplicata{Surface{26};} }

// Outer magnets along the edges
lm = (L_out-w*L)*Cos(Pi/6);
// Surfaces 39-44
Translate {lm*Cos(Pi/6), lm*Sin(Pi/6), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(3*Pi/2), lm*Sin(3*Pi/6), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(5*Pi/6), lm*Sin(5*Pi/6), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(7*Pi/6), lm*Sin(7*Pi/6), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(9*Pi/6), lm*Sin(9*Pi/6), 0} { Duplicata{Surface{26};} }
Translate {lm*Cos(11*Pi/6), lm*Sin(11*Pi/6), 0} { Duplicata{Surface{26};} }

BooleanDifference(45) = { Surface{25}; Delete;}{ Surface{26:44}; };

lc = 0.25;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh 2;
