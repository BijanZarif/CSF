a = 0.75;

L = 30;
h = 6;
r = 4;
R = 4;

s = 10;

Point(1) = {0, 0, 0, a};
Point(2) = {s, 0, 0, 0.5*a};
Point(3) = {s+r, R, 0, 0.5*a};
Point(4) = {s+r, 0, 0, 0.5*a};
Point(5) = {s+2*r, 0, 0, 0.5*a};
Point(6) = {L, 0, 0, a};
Point(7) = {L, h, 0, a};
Point(8) = {0, h, 0, a};

Point(9) = {s, h, 0, 0.5*a};
Point(10) = {s + 2*r, h, 0, 0.5*a};

Point(11) = {0,2*h,0,a};
Point(12) = {L, 2*h,0,a};

Line(1) = {1, 2};

Ellipse(2) = {2, 4, 5, 3};
Ellipse(3) = {5, 4, 2, 3};

Line(4) = {5,6};
Line(5) = {6,7};
Line(6) = {7,10};

Line(7) = {10,9};
Line(8) = {9,8};

Line(9) = {8,1};

Line(10) = {8,11};
Line(11) = {11,12};
Line(12) = {12,7};

Line Loop(13) = {11, 12, 6, 7, 8, 10};
Plane Surface(14) = {13};
Line Loop(15) = {9, 1, 2, -3, 4, 5, 6, 7, 8};
Plane Surface(16) = {15};
