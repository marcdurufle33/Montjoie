lc = 0.5;

//a = 0.0175;
//b = 0.092;
//c = 20.0;

a = 0.02;
b = 0.2;

N1 = Ceil(0.25*Pi*b/lc) + 1;
N2 = Ceil(a/lc) + 1;
N3 = Ceil((b-a)/lc) + 1;
//N4 = Ceil((c-b)/lc) + 1;

R1 = 1.0;
N1 = 4; N2 = 3; N3 = 4;

Point(1) = {a, 0, 0.0, lc};
Point(3) = {a*Sqrt(2)/2, a*Sqrt(2)/2, 0, lc};
Point(4) = {b*Sqrt(2.0)/2, b*Sqrt(2.0)/2, 0.0, lc};
Point(5) = {0, a, 0.0, lc} ;
Point(6) = {0, b, 0.0, lc} ;
Point(7) = {0, 0, 0.0, lc} ;

Point(2) = {0.7*a, 0, 0.0, lc} ;
Point(8) = {0.5*a, 0.5*a, 0.0, lc} ;
Point(9) = {0, 0.7*a, 0.0, lc} ;
Point(10) = {b, 0, 0.0, lc} ;

Line(1) = {1,2};
Line(2) = {3,4};
Line(3) = {5,6};
Line(4) = {2,7};
Circle(5) = {4,7,6};
Circle(6) = {1,7,3};
Circle(7) = {3,7,5};
Line(8) = {8,9};
Line(10) = {8,2};
Line(11) = {7, 9};
Line(12) = {9, 5};
Line(13) = {8, 3};
Line(14) = {1, 10};
Circle(15) = {10, 7, 4};

Transfinite Line{4, 5, 6, 7, 8, 10, 11, 15} = N1;
Transfinite Line{12, 13} = N2;
Transfinite Line{1} = N2;
Transfinite Line{2, 3, 14} = N3 Using Progression R1;


Line Loop(1) = {6, 2, -15, -14};
Line Loop(3) = {1, -10, 13, -6};
Line Loop(4) = {7, -12, -8, 13};
Line Loop(5) = {11, -8, 10, 4};
Line Loop(2) = {2, 5, -3, -7};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

Transfinite Surface{1} = {1, 10, 4, 3};
Transfinite Surface{2} = {3, 4, 6, 5};
Transfinite Surface{3} = {1, 2, 3, 8};
Transfinite Surface{4} = {8, 9, 5, 3};
Transfinite Surface{5} = {7, 2, 8, 9};

Mesh.RecombineAll = 1;

Physical Line(2) = {3, 11, 12} ;
Physical Line(4) = {6, 7} ;
Physical Line(3) = {5, 15};

Physical Surface(1) = {1, 2};
Physical Surface(2) = {3, 4, 5};

SymmetryY = 0.0;
SymmetryX = 0.0;
