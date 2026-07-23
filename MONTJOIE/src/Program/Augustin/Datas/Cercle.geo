lc = 0.5;

a = 0.2;

N1 = 10;

Point(1) = {a, 0, 0.0, lc};
Point(2) = {0, a, 0, lc};
Point(3) = {-a, 0, 0.0, lc};
Point(4) = {0, -a, 0.0, lc} ;
Point(5) = {0, 0, 0.0, lc} ;

Line(1) = {5,1};
Line(2) = {5,2};
Line(3) = {5,3};
Line(4) = {5,4};
Circle(5) = {1,5,2};
Circle(6) = {2,5,3};
Circle(7) = {3,5,4};
Circle(8) = {4,5,1};

Transfinite Line{1,2,3,4,5,6,7,8} = N1;

Curve Loop(1) = {1, 5, -2};
Curve Loop(2) = {2, 6, -3};
Curve Loop(3) = {3, 7, -4};
Curve Loop(4) = {4, 8, -1};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

//Transfinite Surface{1} = {1};
//Transfinite Surface{2} = {2};
//Transfinite Surface{3} = {3};
//Transfinite Surface{4} = {4};

Mesh.RecombineAll = 1;

Physical Line(1) = {5,6,7,8} ;

Physical Surface(1) = {1, 2, 3, 4};
