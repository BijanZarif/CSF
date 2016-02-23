a = 0.001*1000;
x0 = 0.005*1000;
x1 = 0.009*1000;
h = 0.06*1000;
s = 0.003*1000;

Point(1) = {0, 0, 0, a};

Point(2) = {x0, 0, 0,a};
Point(3) = {x1, 0, 0,a};

Point(4) = {0, h, 0, a};
Point(5) = {x0, h, 0,a};
Point(6) = {x1, h, 0,a};

Point(7) = {s,h/6,0,a};
Point(8) = {0,h/6,0,a};
Point(9) = {s,5*h/6,0,a};
Point(10) = {0,5*h/6,0,a};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,6};
Line(4) = {6,5};
Line(5) = {5,4};
Line(6) = {4,10};
Line(7) = {10,8};
Line(8) = {8,1};

Line(9) = {8,7};
Line(10) = {7,9};
Line(11) = {9,10};

Line(12) = {2,5};
Line Loop(13) = {4, -12, 2, 3};
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, -11, -10, -9, 8, 1, 12};
Plane Surface(16) = {15};
Line Loop(17) = {7, 9, 10, 11};
Plane Surface(18) = {17};
