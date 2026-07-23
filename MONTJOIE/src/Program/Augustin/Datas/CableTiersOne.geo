lc = 0.02;

r1 = 0.0175;
r2 = 0.092;
r3 = 0.108;
r4 = 20.0;
Rp = 3.0;

N1 = 10;
N2 = 3;
N3 = 8;

coef = 0.5;

// points

Point(1) = {0, 0, 0, lc};
// frontière extérieure isolant
Point(2) = {r2 * 0.5, - r2 * Sqrt(3) * 0.5, 0, lc};
Point(3) = {r2, 0, 0, lc};
Point(4) = {r2 * 0.5, r2 * Sqrt(3) * 0.5, 0, lc};
// frontière extérieure acier
Point(5) = {r3 * 0.5, - r3 * Sqrt(3) * 0.5, 0, lc};
Point(6) = {r3, 0, 0, lc};
Point(7) = {r3 * 0.5, r3 * Sqrt(3) * 0.5, 0, lc};
// frontière fictive eau
Point(8) = {r4 * 0.5, - r4 * Sqrt(3) * 0.5, 0, lc};
Point(9) = {r4, 0, 0, lc};
Point(10) = {r4 * 0.5, r4 * Sqrt(3) * 0.5, 0, lc};
// points du conducteur
Point(11) = {r1 * 0.5, - r1 * Sqrt(3) * 0.5, 0, lc};
Point(12) = {r1 * 0.5, r1 * Sqrt(3) * 0.5, 0, lc};


// lignes

// fontières circulaires des matériaux
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {5, 1, 6};
Circle(4) = {6, 1, 7};
Circle(5) = {8, 1, 9};
Circle(6) = {9, 1, 10};
// frontière du conducteur
Circle(7) = {11, 1, 12};
// lignes de maillage
Line(11) = {11, 2};
Line(12) = {12, 4};
Line(13) = {2, 5};
Line(14) = {3, 6};
Line(15) = {4, 7};
Line(16) = {5, 8};
Line(17) = {6, 9};
Line(18) = {7, 10};
Line(19) = {1, 11};
Line(20) = {1, 12};

// lignes transfinies sur le maillage


// cercles
Transfinite Line{1,2,3,4,5,6,7} = N1;
// lignes acier
Transfinite Line{13, 14, 15} = N2;
// lignes eau
Transfinite Line{16, 17, 18} = N3 Using Progression Rp;

// plans isolant (1) et conducteur (2)
Line Loop(1) = {2, -12, -7, 11, 1};
Printf("Coucou");
Line Loop(2) = {7, -20, 19};
Printf("Coucou");

Plane Surface(1) = {1};
Plane Surface(2) = {2};

// plans transfinis
Line Loop(3) = {4, -15, -2, 14};
Plane Surface(3) = {3};

Line Loop(4) = {3, -14, -1, 13};
Plane Surface(4) = {4};

Line Loop(5) = {17, 6, -18, -4};
Plane Surface(5) = {5};

Line Loop(6) = {16, 5, -17, -3};
Plane Surface(6) = {6};

Transfinite Surface 3;
Transfinite Surface 4;
Transfinite Surface 5;
Transfinite Surface 6;

Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 0;

Physical Line(1) = {5, 6};
Physical Line(2) = {3, 4};
Physical Line(3) = {1, 2};
Physical Line(4) = {7};
// lignes de coupes haut et bas
Physical Line(5) = {19, 11, 13, 16};
Physical Line(6) = {20, 12, 15, 18};

Physical Surface(3) = {1};
Physical Surface(2) = {2};
Physical Surface(1) = {3, 4};
Physical Surface(4) = {5, 6};
