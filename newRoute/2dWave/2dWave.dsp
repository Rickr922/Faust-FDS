import("stdfaust.lib");

k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);

lambda = c*k/h;

nInputs = inputs(schemeMidPoint);
routing(X, Y, x0, y0) = route(X*Y*4+1, X*Y*4, 
                                par(x, X, par(y, Y, connections(x,y))), 
                                in, C(x0,y0))   
    with {
        in = X*Y*4 + 1; // additional input for signal injection
        connections(x,y) =  N(x,y), S(x,y-1),
                            S(x,y), N(x,y+1),
                            C(x,y), C(x,y),
                            W(x,y), E(x-1,y),
                            E(x,y), W(x+1,y);

        N(x,y) = (1 + 0 + (x+y*X)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        S(x,y) = (1 + 1 + (x+y*X)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        C(x,y) = (1 + 2 + (x+y*X)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        W(x,y) = (1 + 3 + (x+y*X)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        E(x,y) = (1 + 4 + (x+y*X)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    };

model(X,Y, node) = 
        par (x, X, par(y,Y, schemeMidPoint));

schemeMidPoint(u_n,u_s,u,u_w,u_e) = u*2*(1-2*lambda*lambda) - u' + lambda*lambda*(u_e+u_w+u_n+u_s)<:si.bus(nInputs);

outPointX=hslider("output point x",0,0,X*Y-1,1);
outPointY=hslider("output point y",0,0,Y-1,1);

example(X,Y) = (routing(X,Y, 1,1) : model(X, Y, schemeMidPoint)) ~ si.bus(X*Y) : si.bus(X*Y*nInputs):ba.selectn(X*Y*nInputs,outPointX);

process = button("play") : example(6,10);
