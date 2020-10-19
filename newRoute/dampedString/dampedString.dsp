import("stdfaust.lib");

nInputs = 3;

surface4i(X, x0) =   route(X*nInputs+1, X*nInputs, 
                                par(x, X, connections(x)), 
                                in, C(x0))   
    with {
        in = X*nInputs + 1; // additional input for signal injection
        connections(x) =  W(x), E(x-1),
                            C(x), C(x),
                            E(x), W(x+1);

        W(x) = (1 + 0 + (x*nInputs)) * (x>=0) * (x<X);
        C(x) = (1 + 1 + (x*nInputs)) * (x>=0) * (x<X);
        E(x) = (1 + 2 + (x*nInputs)) * (x>=0) * (x<X);
    };

model(X, node) = par (x, X, node)
        //west,
        //par (x, X-2, node),
        //east
    with {
        east(e,c,w) =  c*2*(1-lambda*lambda)/C1+c'*C2/C1+ lambda*lambda*(e+w)/C1,
                       c*2*(1-lambda*lambda)/C1+c'*C2/C1+ lambda*lambda*(e+w)/C1,
                       0;
        west(e,c,w) =  0,
                        c*2*(1-lambda*lambda)/C1+c'*C2/C1+ lambda*lambda*(e+w)/C1,
                       c*2*(1-lambda*lambda)/C1+c'*C2/C1+ lambda*lambda*(e+w)/C1;
    };

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 500;

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

lambda = c*k/h;

node(e,c,w) = c*2*(1-lambda*lambda)/C1+c'*C2/C1+ lambda*lambda*(e+w)/C1, 
c*2*(1-lambda*lambda)/C1+c'*C2/C1+ lambda*lambda*(e+w)/C1, 
c*2*(1-lambda*lambda)/C1+c'*C2/C1+ lambda*lambda*(e+w)/C1;

example(X) = (surface4i(X,45) : model(X, node)) ~ si.bus(X*nInputs) : ba.selectn(X*nInputs,outPoint*nInputs+1);

nPoints = 90;
outPoint = hslider("outPoint",5,1,nPoints,1);
process = button("play") :ba.impulsify : example(nPoints);

//model(4,4, _, node); //listen4(3,3, 1,1); //surface4(2,2); //surface4i(4,4,1,1);


/*

    0   1   2
q
    3   4   5

    6   7   8


    0   1

    2   3

*/
