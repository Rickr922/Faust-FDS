

import("stdfaust.lib");



// simple loopback

straight2(X) = par(i, X*3, _);



// surface4i(X,Y,x,y): creates the connections for a surface of X*Y nodes with 3 inputs and 3 outputs 

// one left, one right and one injection into itself

// with an injection point of coord x,y.

surface2i(X, x0) =   route(X*3+1, X*3, 

                                par(x, X, connections(x)), 

                                in, C(x0),

                                in, C(x0)

                            )   

    with {

        in = X*3 + 1; // additional input for signal injection

        connections(x) =    L(x), R(x-1),

                            C(x), C(x),

                            R(x), L(x+1);

        L(x) = (1 + 0 + (x)*3) * (x>=0) * (x<X);

        C(x) = (1 + 1 + (x)*3) * (x>=0) * (x<X);

        R(x) = (1 + 2 + (x)*3) * (x>=0) * (x<X);



    };



listen4(X, x) =   route(X*3,  1,  

                        L(x), 0,

                        C(x), 1,

                        R(x), 0

                        )   

    with {

        L(x) = (1 + 0 + (x)*3) * (x>=0) * (x<X);

        C(x) = (1 + 1 + (x)*3) * (x>=0) * (x<X);

        R(x) = (1 + 2 + (x)*3) * (x>=0) * (x<X);

    };





// Physical parameters for the system

K = 0.4;

Z = 0.001;

invM = 1;



// "precomputed" parameters for the update scheme

A = 2 - 2 * (K+Z)*invM;

B = 2 * Z *invM - 1;

C = (K+Z)*invM;

D = -Z*invM;

E = invM;





model(X, node) = 

    par (x, 1, left),

    par (y, X-2, node),

    par (x, 1, right)

with {

    right(l, c, r) = 0, 0, 0 ;

    left(l, c, r) = 0, 0, 0 ;

};





node(l, c, r) = A*c + B*c' + C*(l+r) + D*(l'+r') <:_,_,_;



example(X) = (surface2i(X, 60) : model(X, node)) ~ straight2(X) : listen4(X, 4);



process = button("play")* 0.1:ba.impulsify: example(100);


