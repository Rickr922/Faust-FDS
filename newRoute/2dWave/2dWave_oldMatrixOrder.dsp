import("stdfaust.lib");

k = 1/ma.SR;
c = 344;
h = c * k*sqrt(2);
nPointsX = 3;
nPointsY = 3;

lambda = c*k/h;

nInputs = inputs(schemeMidPoint)-2;
routing(X, Y, x0, y0) = route(X*Y+1, X*Y*nInputs, 
                                par(x, X, par(y, Y, connections(x,y))),
                            in, C(x0,y0))
    with {
        in = X*Y + 1; // additional input for signal injection
        connections(x,y) =  P(x,y), S(x,y+1),
                            P(x,y), N(x,y-1),
                            P(x,y), C(x,y),
                            P(x,y), E(x-1,y),
                            P(x,y), W(x+1,y);
        P(x,y) = x*Y+y+1;
    
        W(x,y) = (1 + 0 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        E(x,y) = (1 + 1 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        C(x,y) = (1 + 2 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        S(x,y) = (1 + 3 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        N(x,y) = (1 + 4 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
        
    };

model(X,Y) = 
        par (x, X, par(y,Y, schemeMidPoint));

schemeMidPoint(i,j, uSide_w, uSide_e,u, uSide_s, uSide_n) = u*2*(1-2*lambda*lambda) - u' + lambda*lambda*(uSide_e+uSide_w+uSide_n+uSide_s);

outPoint=hslider("output point x",1,0,nPointsX*nPointsY,1);
outPointY=hslider("output point y",0,0,Y-1,1);

example(X,Y) = (routing(X,Y, 0,0) : par (x, X, par(y,Y, schemeMidPoint(x,y)))) ~ si.bus(X*Y) : ba.selectn(X*Y,outPoint);

process = button("play"):ba.impulsify : example(nPointsX,nPointsY);
//process=model(3,3);
