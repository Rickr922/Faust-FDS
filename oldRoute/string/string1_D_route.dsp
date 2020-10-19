import("stdfaust.lib");

k = 1/ma.SR;
c = 344;
h = c*k;

lambda = c*k/h;

string1dMidPoint(fIn, uSide_l, uSide_r) = u
  letrec
  {
      'u = 2*u-u'+ lambda*lambda*(uSide_l-2*u+uSide_r) + fIn;
  };

updatePoint(fIn, uSide_l, uSide_r) = string1dMidPoint(fIn, uSide_l, uSide_r)<:_,_; //u_(l-1)=uSide_l, u_(l+1)=uSide_r

updatePointOut(fIn, uSide_l, uSide_r) = string1dMidPoint(fIn, uSide_l, uSide_r)<:_,_,_; //u_(l-1)=uSide_l, u_(l+1)=uSide_r


fixedPoint(uSide) = 0;

freePoint(uSide) = u
letrec
{
  'u = 2*u-u'+ lambda*lambda*2*(uSide - u);
};

hit = button("hit"):ba.impulsify;

model =
    (route(17,17,
        (1,2),
        (2,1),
        (3,4),
        (4,3),
        (5,6),
        (6,5),
        (7,9),//out 8 is left free for the force in->everything else jumps one place
        (8,7),
        (9,11),
        (10,10),
        (11,13),
        (12,12),
        (13,15),
        (14,14),
        (15,17),
        (16,16),
        (17,8))://force in
            fixedPoint,
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint,
            updatePoint(0),
            updatePointOut(0),
            updatePoint(0),
            fixedPoint :
                route(17,17, //routeOut
                    (1,1),
                    (2,2),
                    (3,3),
                    (4,4),
                    (5,5),
                    (6,6),
                    (7,7),
                    (8,8),
                    (9,9),
                    (10,10),
                    (11,11),
                    (12,12),
                    (13,13),
                    (14,17),//out signal->routed to last slot
                    (15,14),
                    (16,15),
                    (17,16))
            )~par(i, 16, _): par(i, 16,!), par(i, 1, _);


process = hit:model<:_,_;
