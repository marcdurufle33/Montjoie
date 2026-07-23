lc = 0.02;

r1 = 0.01;
r2 = 1;

N1 = 3;
N3 = 8;
Rp = 2.0;

coef = 0.5;

Point(1) = {0, 0, 0, lc};
Point(2) = {r1, 0, 0, lc};
Point(3) = {0, r1, 0, lc};
Point(4) = {-r1, 0, 0, lc};
Point(5) = {0, -r1, 0, lc};
Point(6) = {r2, 0, 0, lc};
Point(7) = {0, r2, 0, lc};
Point(8) = {-r2, 0, 0, lc};
Point(9) = {0, -r2, 0, lc};

Circle(1) = {4, 1, 5};
Circle(2) = {5, 1, 2};
Circle(3) = {2, 1, 3};
Circle(4) = {3, 1, 4};

Circle(5) = {7, 1, 8};
Circle(6) = {8, 1, 9};
Circle(7) = {9, 1, 6};
Circle(8) = {6, 1, 7};

Line(9) = {3, 7};
Line(10) = {2, 6};
Line(11) = {5, 9};
Line(12) = {4, 8};

Line(13) = {1, 3};
Line(14) = {1, 2};
Line(15) = {1, 5};
Line(16) = {1, 4};

Transfinite Line{1,2,3,4} = N3;
Transfinite Line{5,6,7,8,13,14,15,16} = N3;
Transfinite Line{9, 10, 11, 12} = N1 Using Progression Rp;

Curve Loop(1) = {14, 3, -13};
Plane Surface(1) = {1};
Curve Loop(2) = {13, 4, -16};
Plane Surface(2) = {2};
Curve Loop(3) = {16, 1, -15};
Plane Surface(3) = {3};
Curve Loop(4) = {15, 2, -14};
Plane Surface(4) = {4};
Curve Loop(5) = {5, -12, -4, 9};
Plane Surface(5) = {5};
Curve Loop(6) = {6, -11, -1, 12};
Plane Surface(6) = {6};
Curve Loop(7) = {7, -10, -2, 11};
Plane Surface(7) = {7};
Curve Loop(8) = {8, -9, -3, 10};
Plane Surface(8) = {8};


Transfinite Surface 1;
Transfinite Surface 2;
Transfinite Surface 3;
Transfinite Surface 4;
Transfinite Surface 5;
Transfinite Surface 6;
Transfinite Surface 7;
Transfinite Surface 8;

Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 0;

Physical Line(1) = {1, 2, 3, 4};
Physical Line(2) = {5, 6, 7, 8};

Physical Surface(1) = {1, 2, 3, 4};
Physical Surface(2) = {5, 6, 7, 8};
