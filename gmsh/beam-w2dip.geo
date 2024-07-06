// Gmsh project created on Mon Nov 13 13:27:42 2023
SetFactory("OpenCASCADE");

// Units are in cm
lm = 1; // magnet diameter, same as thickness 
lb = 1.5; // case cube length

// Beam dimensions (excluding magnets + case)
r = 0.02; // w/L ratio

// (Pi+1)(lb-w) is the minimum length to indent beam with cylinder d = lb-w
w = lb / (1 + 1/(r*(Pi+1))); // width < lb
L = w/r;

// Elastic beam
Rectangle(1) = {-L/2, -w/2, 0, L, w, 0};

// Cylindrical magnets
Disk(2) = {-L/2-lb/2, 0, 0, lm/2};
Disk(3) = {L/2+lb/2, 0, 0, lm/2};

// Magnet case
Rectangle(4) = {-L/2-lb, -lb/2, 0, lb, lb, 0};
Rectangle(5) = {L/2, -lb/2, 0, lb, lb, 0};

// Get elastic part of beam
BooleanDifference{ Surface{1}; Delete; }{ Surface{4,5}; }
BooleanDifference{ Surface{4}; Delete; }{ Surface{2}; }
BooleanDifference{ Surface{5}; Delete; }{ Surface{3}; }

// Make sure nodes are the same for each shared edge
BooleanFragments{ Surface{:}; Delete; }{}
BooleanFragments{ Curve{:}; Delete; }{}
BooleanFragments{ Point{:}; Delete; }{}
Coherence;
lc = 0.4;

Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFactor = lc;
Mesh.MeshSizeMin = lc/16;
Mesh.MeshSizeMax = lc;
Mesh.Algorithm = 5;

// Refine mesh on elastic segments
Field[1] = Box;
Field[1].VIn = lc/8;
Field[1].VOut = lc;
Field[1].XMin = -L/2;
Field[1].XMax = L/2;
Field[1].YMin = -w/2;
Field[1].YMax = w/2;
Field[1].Thickness = 0.5;

Background Field = 1;
Mesh 2;

Physical Surface("ela", 1) = {1};
Physical Surface("M1", 2) = {2};
Physical Surface("M2", 3) = {3};
Physical Surface("C1", 4) = {4};
Physical Surface("C2", 6) = {5};
