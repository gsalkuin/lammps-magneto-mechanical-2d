// Gmsh project created on Mon Nov 13 13:27:42 2023
SetFactory("OpenCASCADE");
Geometry.CopyMeshingMethod = 1;

L = 15; // mm
w = 0.1; // ratio of thickness to beam length

// NOTE: For achiral, double the length
// Unit square
Rectangle(1) = {0, 0, 0, 2*L, 2*L, 0};
Rectangle(2) = {w*L, w*L, 0, 2*L-2*w*L, 2*L-2*w*L, 0};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

// Add magnet slots
// Surfaces 4-11
Rectangle(4) = {-w*L, -w*L, 0, 3*w*L, 3*w*L, 0};
Translate {L*(1-w/2), 0, 0} { Duplicata{ Surface{4}; } }
Translate {2*L*(1-w/2), 0, 0} { Duplicata{ Surface{4}; } }

Translate {0, L*(1-w/2), 0} { Duplicata{ Surface{4}; } }
Translate {2*L*(1-w/2), L*(1-w/2), 0} { Duplicata{ Surface{4}; } }

Translate {0, 2*L*(1-w/2), 0} { Duplicata{ Surface{4:6}; } }

// Elastic part
BooleanUnion(12) = { Surface{3}; Delete; }{ Surface{4:11}; Delete; };

// Duplicate with overlap
Translate {2*L-w*L, 0, 0} { Duplicata{ Surface{12}; } }
Translate {0, 2*L-w*L, 0} { Duplicata{ Surface{12}; } }
Translate {2*L-w*L, 2*L-w*L, 0} { Duplicata{ Surface{12}; } }

// Complete Elastic part
BooleanUnion(16) = { Surface{12}; Delete; }{ Surface{13:15}; Delete; };

// Add magnets
// Surfaces 17-24
Rectangle(17) = {-0.5*w*L, -0.5*w*L, 0, 2*w*L, 2*w*L, 0};
Translate {L*(1-w/2), 0, 0} { Duplicata{ Surface{17}; } }
Translate {2*L*(1-w/2), 0, 0} { Duplicata{ Surface{17}; } }
Translate {3*L*(1-w/2), 0, 0} { Duplicata{ Surface{17}; } }
Translate {4*L*(1-w/2), 0, 0} { Duplicata{ Surface{17}; } }

Translate {0, L*(1-w/2), 0} { Duplicata{ Surface{17}; } }
Translate {2*L*(1-w/2), L*(1-w/2), 0} { Duplicata{ Surface{17}; } }
Translate {4*L*(1-w/2), L*(1-w/2), 0} { Duplicata{ Surface{17}; } }

// surf 25-29
Translate {0, 2*L*(1-w/2), 0} { Duplicata{ Surface{17:21}; } }

// surf 30-32
Translate {0, 2*L*(1-w/2), 0} { Duplicata{ Surface{22:24}; } }

// surf 33-37
Translate {0, 2*L*(1-w/2), 0} { Duplicata{ Surface{25:29}; } }

BooleanDifference(38) = { Surface{16}; Delete; }{ Surface{17:37}; };

lc = 0.25;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh 2;

Recursive Delete {
  Surface{18}; Surface{20}; Surface{22}; Surface{23}; Surface{24}; Surface{26}; Surface{28}; Surface{30}; Surface{31}; Surface{32}; Surface{34}; Surface{36}; 
}

Mesh 2;
