// Gmsh project created on Mon Nov 13 13:27:42 2023
SetFactory("OpenCASCADE");
Geometry.CopyMeshingMethod = 1;

L = 15; // mm
w = 0.1; // ratio of thickness to beam length

// Unit square
Rectangle(1) = {0, 0, 0, L, L, 0};
Rectangle(2) = {w*L, w*L, 0, L-2*w*L, L-2*w*L, 0};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };

// Duplicate with overlap
// Surfaces 4-11
Translate {L-w*L, 0, 0} { Duplicata{ Surface{3}; } }
Translate {2*(L-w*L), 0, 0} { Duplicata{ Surface{3}; } }
Translate {0, L-w*L, 0} { Duplicata{ Surface{3:5}; } }
Translate {0, 2*(L-w*L), 0} { Duplicata{ Surface{3:5}; } }

// Elastic part
BooleanUnion(12) = { Surface{3}; Delete; }{ Surface{4:11}; Delete; };

// Add magnet slots
// Surfaces 13-28
Rectangle(13) = {-w*L, -w*L, 0, 3*w*L, 3*w*L, 0};
Translate {L-w*L, 0, 0} { Duplicata{ Surface{13}; } }
Translate {2*(L-w*L), 0, 0} { Duplicata{ Surface{13}; } }
Translate {3*(L-w*L), 0, 0} { Duplicata{ Surface{13}; } }
Translate {0, L-w*L, 0} { Duplicata{ Surface{13:16}; } }
Translate {0, 2*(L-w*L), 0} { Duplicata{ Surface{13:16}; } }
Translate {0, 3*(L-w*L), 0} { Duplicata{ Surface{13:16}; } }

BooleanUnion(29) = { Surface{12}; Delete; }{ Surface{13:28}; Delete; };

// Add magnets
// Surfaces 30-45
Rectangle(30) = {-0.5*w*L, -0.5*w*L, 0, 2*w*L, 2*w*L, 0};
Translate {L-w*L, 0, 0} { Duplicata{ Surface{30}; } }
Translate {2*(L-w*L), 0, 0} { Duplicata{ Surface{30}; } }
Translate {3*(L-w*L), 0, 0} { Duplicata{ Surface{30}; } }
Translate {0, L-w*L, 0} { Duplicata{ Surface{30:33}; } }
Translate {0, 2*(L-w*L), 0} { Duplicata{ Surface{30:33}; } }
Translate {0, 3*(L-w*L), 0} { Duplicata{ Surface{30:33}; } }

BooleanDifference(46) = { Surface{29}; Delete; }{ Surface{30:45}; };

lc = 0.25;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;
Mesh 2;
