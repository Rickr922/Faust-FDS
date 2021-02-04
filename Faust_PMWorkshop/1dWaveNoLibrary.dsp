import("stdfaust.lib");

k = 1/ma.SR;
c=344;
h=c*k;
lambda=c*k/h;

/*
u_l^n+1 = 2(1-lambda^2)u_l^n - u_l^n-1 + lambda^2(u_l+1^n + u_l-1^n)
*/

updateEq(fIn, u_w, u, u_e) = 2*(1-lambda^2)*u - u' + lambda^2*(u_w+u_e) + fIn;
nPoints= 3;

build1DScheme(points) = par(i,points,updateEq);

routing = route(6,12, (1,3), (1,6), (2,4), (2,7), (2,10), (3,8), (3,11), (4,1), (5,5), (6,9), (0,2), (0,12));

model(points) = (routing : build1DScheme(points))~si.bus(points);

play = button("Play") : ba.impulsify;
inPoint = hslider("Input Point", 1,0,nPoints-1,1);
outPoint = hslider("Output Point", 1,0,nPoints-1,1);

pointSelectorIn(points,point) = ba.selectoutn(points,point);

pointSelectorOut(points,point) = ba.selectn(points,point);


process = play: pointSelectorIn(nPoints,inPoint) : model(nPoints) : pointSelectorOut(nPoints, outPoint);

