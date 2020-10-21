import("stdfaust.lib");

//--------------------------------Model Settings-----------------------------//
k = 1/ma.SR;
K = 20;
s0 = 2;
s1 = 0.05;
c=344;

coeff = c^2*k^2+4*s1*k;
h = sqrt((coeff+sqrt(coeff*coeff+16*K*K*k*k)));

nPointsX = 10;
nPointsY = 10;

lambda = c*k/h;

//----------------------------------Equations--------------------------------//
mu=K*K*k*k/(h^4);
den = 1+s0*k;
A = 2*(1-10*mu-2*lambda*lambda-4*s1*k*k)/den;
B = (s0*k+4*k*k-1)/den;
C = (8*mu + lambda*lambda + 2*s1*k*k)/den;
D = -2*mu/den;
E = -mu/den;
F = -2*s1*k*k/den;

schemeMidPoint(fIn,u_nn,u_n,u_ss,u_s,u,u_w,u_ww,u_e,u_ee,u_nw,u_ne,u_sw,u_se) =
    A*u + B*u' + C*(u_n+u_s+u_w+u_e) + D*(u_nw+u_ne+u_sw+u_se) + E*(u_nn+u_ss+u_ww+u_ee) + F*(u_n'+u_s'+u_w'+u_e') + fIn;

buildScheme(X,Y) = par (x, X,
                par(y,Y, schemeMidPoint));

//----------------------------------Controls---------------------------------//
inPointX=hslider("input point x", floor(nPointsX/2),0,nPointsX-1,1);
inPointY=hslider("input point y", floor(nPointsY/2),0,nPointsY-1,1);
outPointX=hslider("output point x",floor(nPointsX/2),0,nPointsX-1,1);
outPointY=hslider("output point y",floor(nPointsY/2),0,nPointsY-1,1);
forceModel = button("play") : ba.impulsify/3;
stop = button("Stop");

//----------------------------------Interpolation---------------------------------//
linInterpolation2D(nPointsX,nPointsY,pointX,pointY) =
par(i,nPointsX,
    par(j,nPointsY,_*
        select2((i==intX) & (j==intY),
            select2((i==(intX+1)) & (j==intY),
                select2((i==intX) & (j==(intY+1)),
                    select2((i==(intX+1)) & (j==(intY+1)),
                        0,
                        fractionX*fractionY),
                    (1-fractionX)*fractionY),
                fractionX*(1-fractionY)),
            (1-fractionX)*(1-fractionY))))
with
{
    fractionX = ma.frac(pointX);
    fractionY = ma.frac(pointY);
    intX = int(pointX);
    intY = int(pointY);
};

//----------------------------------Force---------------------------------//
stairsForce(X,Y,pointX,pointY) = ba.selectoutn(X*Y,pointY+pointX*Y);

//----------------------------------Output-------------------------------//
stairsOutput(X,Y,pointX,pointY) = ba.selectn(X*Y,pointY+pointX*Y);
linInterpolation2DOut(nPointsX,nPointsY,pointX,pointY) = linInterpolation2D(nPointsX,nPointsY,pointX,pointY):>_;
//----------------------------------Build Model-------------------------------//
nInputs = inputs(schemeMidPoint);
routing(X, Y) = route(X*Y*2, X*Y*nInputs,
                                par(x, X, par(y, Y, connections(x,y))))
with
{
    connections(x,y) =  P(x,y) + X*Y, F(x,y),
                        P(x,y), NN(x,y+2),
                        P(x,y), N(x,y+1),
                        P(x,y), SS(x,y-2),
                        P(x,y), S(x,y-1),
                        P(x,y), C(x,y),
                        P(x,y), WW(x+2,y),
                        P(x,y), W(x+1,y),
                        P(x,y), EE(x-2,y),
                        P(x,y), E(x-1,y),
                        P(x,y), NW(x+1,y+1),
                        P(x,y), NE(x-1,y+1),
                        P(x,y), SW(x+1,y-1),
                        P(x,y), SE(x-1,y-1);
    P(x,y) = x*Y+y+1;

    F(x,y)  = (1 + 0 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    NN(x,y) = (1 + 1 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    N(x,y)  = (1 + 2 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    SS(x,y) = (1 + 3 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    S(x,y)  = (1 + 4 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    C(x,y)  = (1 + 5 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    W(x,y)  = (1 + 6 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    WW(x,y) = (1 + 7 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    E(x,y)  = (1 + 8 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    EE(x,y) = (1 + 9 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    NW(x,y) = (1 + 10 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    NE(x,y) = (1 + 11 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    SW(x,y) = (1 + 12 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
    SE(x,y) = (1 + 13 + (x*Y+y)*nInputs) * (x>=0) * (x<X) * (y>=0) * (y<Y);
};

model(X,Y) = (routing(X,Y) : buildScheme(X,Y)) ~ si.bus(X*Y); //par(i,X*Y,_*(stop==0));*/


process = forceModel <: stairsForce(nPointsX,nPointsY,inPointX,inPointY) : model(nPointsX,nPointsY) : stairsOutput(nPointsX,nPointsY,outPointX,outPointY)<:_,_;

//process=B;
