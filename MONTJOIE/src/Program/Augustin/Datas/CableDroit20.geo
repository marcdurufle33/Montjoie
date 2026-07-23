lc = 0.02;

r0 = 0.092;
r1 = 0.108;
r2 = 20;
R = 0.05;
rCond = 0.0175;

N1 = 4;
N2 = 3;
N3 = 8;
NCond = 4;
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
// points des conducteurs
Point(14) = {R, 0, 0, coef * lc};
Point(15) = {R+rCond, 0, 0, coef * lc};
Point(16) = {R, rCond, 0, coef * lc};
Point(17) = {R-rCond, 0, 0, coef * lc};
Point(18) = {R, -rCond, 0, coef * lc};

cX = R * Cos(2*Pi/3);
cY = R * Sin(2*Pi/3);
Point(19) = {cX, cY, 0, coef * lc};
Point(20) = {cX+rCond, cY, 0, coef * lc};
Point(21) = {cX, cY + rCond, 0, coef * lc};
Point(22) = {cX-rCond, cY, 0, coef * lc};
Point(23) = {cX, cY - rCond, 0, coef * lc};

cX = R * Cos(-2*Pi/3);
cY = R * Sin(-2*Pi/3);
Point(24) = {cX, cY, 0, coef * lc};
Point(25) = {cX+rCond, cY, 0, coef * lc};
Point(26) = {cX, cY + rCond, 0, coef * lc};
Point(27) = {cX-rCond, cY, 0, coef * lc};
Point(28) = {cX, cY - rCond, 0, coef * lc};

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


// frontière des conducteurs
Circle(21) = {15, 14, 16};
Circle(22) = {16, 14, 17};
Circle(23) = {17, 14, 18};
Circle(24) = {18, 14, 15};

Circle(25) = {20, 19, 21};
Circle(26) = {21, 19, 22};
Circle(27) = {22, 19, 23};
Circle(28) = {23, 19, 20};

Circle(29) = {25, 24, 26};
Circle(30) = {26, 24, 27};
Circle(31) = {27, 24, 28};
Circle(32) = {28, 24, 25};

Transfinite Line{1,2,3,4} = N3;
Transfinite Line{5,6,7,8} = N3;
Transfinite Line{13,14,15,16} = N3;
Transfinite Line{17,18,19,20} = N2;
Transfinite Line{21,22,23,24} = NCond;
Transfinite Line{25,26,27,28} = NCond;
Transfinite Line{29,30,31,32} = NCond;
Transfinite Line{9, 10, 11, 12} = N1 Using Progression Rp;
//Transfinite Line{9, 10, 11, 12} = 6;

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

Curve Loop(10) = {21, 22, 23, 24};
Plane Surface(10) = {10};
Curve Loop(11) = {25, 26, 27, 28};
Plane Surface(11) = {11};
Curve Loop(12) = {29, 30, 31, 32};
Plane Surface(12) = {12};

Curve Loop(1) = {14, 15, 16, 13};
Plane Surface(1) = {1, 10, 11, 12};


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
Physical Line(4) = {21, 22, 23, 24};
Physical Line(5) = {25, 26, 27, 28};
Physical Line(6) = {29, 30, 31, 32};

Physical Surface(1) = {1};
Physical Surface(2) = {2, 3, 4, 5};
Physical Surface(3) = {6, 7, 8, 9};
Physical Surface(4) = {10};
Physical Surface(5) = {11};
Physical Surface(6) = {12};
