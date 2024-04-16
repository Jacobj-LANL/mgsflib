// Gmsh project created on Thu Oct 19 12:18:40 2023
cl__1 = 0.1;
Point(1) = {0.5, 0.5, 0.5, cl__1};
Point(2) = {0.55, 0.5, 0.5, cl__1};
Point(3) = {0.45, 0.5, 0.5, cl__1};
Point(4) = {0.5, 0.45, 0.5, cl__1};
Point(5) = {0.5, 0.55, 0.5, cl__1};
Point(6) = {0.5, 0.5, 0.45, cl__1};
Point(7) = {0.5, 0.5, 0.55, cl__1};
//+
Circle(1) = {7, 1, 2};
//+
Circle(2) = {2, 1, 6};
//+
Circle(3) = {6, 1, 3};
//+
Circle(4) = {3, 1, 7};
//+
Circle(5) = {7, 1, 4};
//+
Circle(6) = {4, 1, 6};
//+
Circle(7) = {6, 1, 5};
//+
Circle(8) = {5, 1, 7};
//+
Circle(9) = {2, 1, 4};
//+
Circle(10) = {4, 1, 3};
//+
Circle(11) = {3, 1, 5};
//+
Circle(12) = {5, 1, 2};
//+
Curve Loop(1) = {1, -12, 8};
//+
Surface(1) = {1};
//+
Curve Loop(2) = {2, 7, 12};
//+
Surface(2) = {2};
//+
Curve Loop(3) = {7, -11, -3};
//+
Surface(3) = {3};
//+
Curve Loop(4) = {11, 8, -4};
//+
Surface(4) = {4};
//+
Curve Loop(5) = {1, 9, -5};
//+
Surface(5) = {5};
//+
Curve Loop(6) = {9, 6, -2};
//+
Surface(6) = {6};
//+
Curve Loop(7) = {6, 3, -10};
//+
Surface(7) = {7};
//+
Curve Loop(8) = {4, 5, 10};
//+
Surface(8) = {8};
//+
Surface Loop(1) = {3, 2, 6, 5, 1, 4, 8, 7};
//+
Volume(1) = {1};
