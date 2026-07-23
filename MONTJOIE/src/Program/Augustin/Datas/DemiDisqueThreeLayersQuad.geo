//a = 1.0;
//b=0.2;
a = 10;
b = 0.108;
c = 0.092;
d = 0.0175;

h = 1;

Point(1) = {0,0,0,0.5};
Point(2) = {a,0,0,h};
Point(3) = {-a,0,0,h};

Point(4) = {b,0,0,h};
Point(5) = {-b,0,0,h};

Point(6) = {c,0,0,h};
Point(7) = {-c,0,0,h};

Point(8) = {d,0,0,h};
Point(9) = {-d,0,0,h};
Point(10) = {0,d,0,h};
Point(11) = {0,-d,0,h};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 2};

Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 4};

Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 6};

Circle(7) = {8, 1, 10};
Circle(8) = {10, 1, 9};
Circle(19) = {9, 1, 11};
Circle(20) = {11, 1, 8};

// lignes pour le tranfinite line

Line(9) = {2,4};
Line(10) = {4,6};
Line(11) = {6,8};
Line(12) = {3,5};
Line(13) = {5,7};
Line(14) = {7,9};

Line(15) = {8,1};
Line(16) = {9,1};
Line(17) = {10,1};
Line(18) = {11,1};


Transfinite Line {11, 14, 10, 13} =  3;
Transfinite Line {15,16, 17,18} =  5;
Transfinite Line {-9, -12} =  4 Using Progression 3;

Transfinite Line {1,2,3,4,5,6,7,8,19,20} = 5;

Curve Loop(1) = {9, -4, -12, 2};
Plane Surface(1) = {1};
Curve Loop(2) = {9, 3, -12, -1};
Plane Surface(2) = {2};
Curve Loop(3) = {10, -6, -13, 4};
Plane Surface(3) = {3};
Curve Loop(4) = {10, 5, -13, -3};
Plane Surface(4) = {4};
Curve Loop(5) = {11, -20, -19, -14, 6};
Plane Surface(5) = {5};
Curve Loop(6) = {11, 7, 8, -14, -5};
Plane Surface(6) = {6};

Curve Loop(7) = {-15, 7, 17};
Plane Surface(7) = {7};
Curve Loop(8) = {-17, 8, 16};
Plane Surface(8) = {8};
Curve Loop(9) = {-16, 19, 18};
Plane Surface(9) = {9};
Curve Loop(10) = {-18, 20, 15};
Plane Surface(10) = {10};

//Point{1} In Surface{10};

Transfinite Surface 1;
Transfinite Surface 2;
Transfinite Surface 3;
Transfinite Surface 4;
Transfinite Surface 5;
Transfinite Surface 6;
Transfinite Surface 7;
Transfinite Surface 8;
Transfinite Surface 9;
Transfinite Surface 10;

Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 0;

Physical Line(1) = {7, 8};
Physical Line(2) = {5, 6};
Physical Line(3) = {3, 4};
Physical Line(5) = {1, 2};

Physical Surface(2) = {7, 8, 9, 10};
Physical Surface(3) = {5,6};
Physical Surface(1) = {3, 4};
Physical Surface(4) = {1, 2};
