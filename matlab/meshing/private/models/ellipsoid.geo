// Gmsh project created on Thu Oct 19 12:18:40 2023
cl__1 = 0.1;
Point(1) = {0.5, 0.5, 0.5, cl__1};
Point(2) = {0.7, 0.5, 0.5, cl__1};
Point(3) = {0.3, 0.5, 0.5, cl__1};
Point(4) = {0.5, 0.35, 0.5, cl__1};
Point(5) = {0.5, 0.65, 0.5, cl__1};
Point(6) = {0.5, 0.5, 0.4, cl__1};
Point(7) = {0.5, 0.5, 0.6, cl__1};
Ellipse(1) = {4, 1, 2, 2};
Ellipse(2) = {2, 1, 1, 5};
Ellipse(3) = {5, 1, 3, 3};
Ellipse(4) = {3, 1, 3, 4};
Ellipse(5) = {2, 1, 2, 7};
Ellipse(6) = {7, 1, 3, 3};
Ellipse(7) = {3, 1, 3, 6};
Ellipse(8) = {6, 1, 2, 2};
Ellipse(9) = {4, 1, 4, 6};
Ellipse(10) = {6, 1, 5, 5};
Ellipse(11) = {5, 1, 5, 7};
Ellipse(12) = {7, 1, 4, 4};
Curve Loop(1) = {5, 12, 1};
Surface(1) = {1};
Curve Loop(2) = {5, -11, -2};
Surface(2) = {2};
Curve Loop(3) = {3, -6, -11};
Surface(3) = {3};
Curve Loop(4) = {6, 4, -12};
Surface(4) = {4};
Curve Loop(5) = {1, -8, -9};
Surface(5) = {5};
Curve Loop(6) = {8, 2, -10};
Surface(6) = {6};
Curve Loop(7) = {3, 7, 10};
Surface(7) = {7};
Curve Loop(8) = {4, 9, -7};
Surface(8) = {8};
Surface Loop(1) = {6, 5, 1, 2, 3, 7, 8, 4};
Volume(1) = {1};