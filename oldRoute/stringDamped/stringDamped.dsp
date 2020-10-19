import("stdfaust.lib");

k = 1/ma.SR;
c = 344;
h = c*k;
s0 = 1000;

C1 = 1+(2*s0*k*k/h);
C2 = (2*s0*k*k/h)-1;

lambda = c*k/h;

stringDampMidPoint(fIn, uSide_l, uSide_r) = u //u_(l-1)=uSide_l, u_(l+1)=uSide_
  letrec
  {
      'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*(uSide_l+uSide_r)/C1 + fIn;
  };

//------------------------------------------Boundaries--------------------------------------
stringDampFreePoint(uSide) = u
letrec
{
  'u = u*2*(1-lambda*lambda)/C1+u'*C2/C1+ lambda*lambda*2*uSide/C1;
};

fixedPoint(uSide) = 0;

hit = button("hit"):ba.impulsify;

model =
    (route(37,37,
        (1,2),
        (2,1),(3,4),
        (4,3),(5,6),
        (6,5),(7,8),
        (8,7),(9,10),
        (10,9),(11,12),
        (12,11),(13,14),
        (14,13),(15,16),
        (16,15),(17,19),//out 18 left free for the force in->everything else jumps one place
        (18,17),(19,21),
        (20,20),(21,23),
        (22,22),(23,25),
        (24,24),(25,27),
        (26,26),(27,29),
        (28,28),(29,31),
        (30,30),(31,33),
        (32,32),(33,35),
        (34,34),(35,37),
        (36,36),(37,18))://force in
            stringDampFreePoint,
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint   <:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            (stringDampMidPoint(0)<:_,_),
            fixedPoint :
                route(37,37, //routeOut
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
                    (14,14),
                    (15,15),
                    (16,16),
                    (17,17),
                    (18,18),
                    (19,19),
                    (20,20),
                    (21,21),
                    (22,22),
                    (23,23),
                    (24,24),
                    (25,25),
                    (26,26),
                    (27,27),
                        (28,37),//out signal->routed to last slot
                    (29,28),
                    (30,29),
                    (31,30),
                    (32,31),
                    (33,32),
                    (34,33),
                    (35,34),
                    (36,35),
                    (37,36))
            )~par(i, 36, _): par(i, 36,!), par(i, 1, _);

process = hit:model<:_,_;
