a = 1;
h0 = 2;
h1 = 10;
L = 30;

Point(1) = {0,-h0,0,a};
Point(2) = {0,h0,0,a};
Point(3) = {L,-h1,0,a};
Point(4) = {L,h1,0,a};

Point(5) = {-L/10,-h0,0,a};
Point(6) = {-L/10,h0,0,a};
Point(7) = {-L,-h1,0,a};
Point(8) = {-L,h1,0,a};


Line(1) = {1, 3};
Line(2) = {3,4};
Line(3) = {4,2};
Line(4) = {2,1};

//Line(5) = {2,6};
//Line(6) = {5,7};
//Line(7) = {7,8};
//Line(8) = {8,6};
//Line(9) = {5,1};

//Line Loop(10) = {8, -5, -3, -2, -1, -9, 6, 7};
Line Loop(10) = {3, 4, 1, 2};
Plane Surface(11) = {10};

