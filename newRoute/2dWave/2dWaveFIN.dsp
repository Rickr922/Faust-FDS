import("stdfaust.lib");

k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 30;
nPointsY = 30;

lambda = c*k/h;

nInputs = inputs(schemeMidPoint);
routing(X, Y, x0, y0) = route(X*Y*2, X*Y*nInputs,
                                par(x, X, par(y, Y, connections(x,y))))
    with {
        in = X*Y + 1; // additional input for signal injection
        connections(x,y) =  P(x,y) + X*Y, F(x,y),
                            P(x,y), S(x,y-1),
                            P(x,y), N(x,y+1),
                            P(x,y), C(x,y),
                            P(x,y), E(x-1,y),
                            P(x,y), W(x+1,y);
        P(x,y) = x*Y+y+1;

        F(x,y) = (1 + 0 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        N(x,y) = (1 + 1 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        S(x,y) = (1 + 2 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        C(x,y) = (1 + 3 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        W(x,y) = (1 + 4 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        E(x,y) = (1 + 5 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    };

model(X,Y) =
        par (x, X, par(y,Y, schemeMidPoint));

schemeMidPoint(fIn,u_n,u_s,u,u_w,u_e) = u*2*(1-2*lambda*lambda) - u' + lambda*lambda*(u_e+u_w+u_n+u_s) + fIn;

inPoint=hslider("input point",1,0,nPointsX*nPointsY,1);
outPoint=hslider("output point",1,0,nPointsX*nPointsY,1);
outPointY=hslider("output point y",0,0,Y-1,1);
forceModel = button("play") : ba.impulsify;

example(X,Y) = (routing(X,Y, 2,2) : par (x, X, par(y,Y, schemeMidPoint))) ~ si.bus(X*Y) : si.bus(X*Y):ba.selectn(X*Y,outPoint);

process = forceModel : ba.selectoutn(nPointsX*nPointsY,inPoint) : example(nPointsX,nPointsY);
//process=model(3,3);
