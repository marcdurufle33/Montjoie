lc = 0.02;

r0 = 0.18;
r1 = 0.2;
r2 = 0.6;

N1 = 4;
N2 = 3;
N3 = 8;
Rp = 3.0;

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
Point(10) = {r0, 0, 0, lc};
Point(11) = {0, r0, 0, lc};
Point(12) = {-r0, 0, 0, lc};
Point(13) = {0, -r0, 0, lc};

Circle(1) = {4, 1, 5};
Circle(2) = {5, 1, 2};
Circle(3) = {2, 1, 3};
Circle(4) = {3, 1, 4};

Circle(5) = {7, 1, 8};
Circle(6) = {8, 1, 9};
Circle(7) = {9, 1, 6};
Circle(8) = {6, 1, 7};

Circle(13) = {10, 1, 11};
Circle(14) = {11, 1, 12};
Circle(15) = {12, 1, 13};
Circle(16) = {13, 1, 10};

Line(9) = {3, 7};
Line(10) = {2, 6};
Line(11) = {5, 9};
Line(12) = {4, 8};
Line(17) = {11, 3};
Line(18) = {10, 2};
Line(19) = {13, 5};
Line(20) = {12, 4};

Transfinite Line{1,2,3,4} = N3;
Transfinite Line{5,6,7,8} = N3;
Transfinite Line{13,14,15,16} = N3;
Transfinite Line{17,18,19,20} = N2;
Transfinite Line{9, 10, 11, 12} = N1 Using Progression Rp;

Curve Loop(1) = {14, 15, 16, 13};
Plane Surface(1) = {1};
Curve Loop(2) = {4, -20, -14, 17};
Plane Surface(2) = {2};
Curve Loop(3) = {1, -19, -15, 20};
Plane Surface(3) = {3};
Curve Loop(4) = {2, -18, -16, 19};
Plane Surface(4) = {4};
Curve Loop(5) = {3, -17, -13, 18};
Plane Surface(5) = {5};
Curve Loop(6) = {5, -12, -4, 9};
Plane Surface(6) = {6};
Curve Loop(7) = {6, -11, -1, 12};
Plane Surface(7) = {7};
Curve Loop(8) = {7, -10, -2, 11};
Plane Surface(8) = {8};
Curve Loop(9) = {8, -9, -3, 10};
Plane Surface(9) = {9};


Transfinite Surface 2;
Transfinite Surface 3;
Transfinite Surface 4;
Transfinite Surface 5;
Transfinite Surface 6;
Transfinite Surface 7;
Transfinite Surface 8;
Transfinite Surface 9;

Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 0;

Physical Line(2) = {1, 2, 3, 4};
Physical Line(3) = {5, 6, 7, 8};
Physical Line(1) = {13, 14, 15, 16};

Physical Surface(1) = {1};
Physical Surface(2) = {2, 3, 4, 5};
Physical Surface(3) = {6, 7, 8, 9};
