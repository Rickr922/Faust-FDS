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

/*
NOTATION: u_n_m
  n indicates the number of element (update or fixed point)
  m indicates the connection: 1 is left connection, 2 is right connections
  Ex: the extreme left point 0 gives out only x02 because it only has a right
    connection, a generic element gives xn1, xn2, the extreme right point N
    gives out only xN1 because it has only a left connection
*/

//routeOut takes the output and puts it as the last connection

/*routeOut(
    u02,
    u11,u12,
    u21,u22,
    u31,u32,
    u41,u42,
    u51,u52,
    u61,u62,uOut,
    u71,u72,
    u81) = u02,u11,u12,u21,u22,u31,u32,u41,u42,u51,u52,u61,u62,u71,u72,u81,uOut;*/

//routeIn needs to switch all the connections and give them to the proper elements
/*routeIn(u02,u11,u12,u21,u22,u31,u32,u41,u42,u51,u52,u61,u62,u71,u72,u81,fIn) =
    u11,
    u02,u21,
    u12,u31,
    u22,u41,
    fIn,u32,u51,
    u42,u61,
    u52,u71,
    u62,u81,
        u72;*/

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
            fixedPoint,
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint,
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
            updatePointOut(0),
            updatePoint(0),
            updatePoint(0),
            updatePoint(0),
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
