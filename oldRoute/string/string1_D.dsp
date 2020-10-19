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
routeOut(
    u02,
    u11,u12,
    u21,u22,
    u31,u32,
    u41,u42,
    u51,u52,
    u61,u62,uOut,
    u71,u72,
    u81) = u02,u11,u12,u21,u22,u31,u32,u41,u42,u51,u52,u61,u62,u71,u72,u81,uOut;
//routeIn needs to switch all the connections and give them to the proper elements
routeIn(u02,u11,u12,u21,u22,u31,u32,u41,u42,u51,u52,u61,u62,u71,u72,u81,fIn) =
    u11,
    u02,u21,
    u12,u31,
    u22,u41,
    fIn,u32,u51,
    u42,u61,
    u52,u71,
    u62,u81,
        u72;

hit = button("hit"):ba.impulsify;

model =
    (routeIn:
        fixedPoint,
        updatePoint(0),
        updatePoint(0),
        updatePoint(0),
        updatePoint,
        updatePoint(0),
        updatePointOut(0),
        updatePoint(0),
        fixedPoint : routeOut)~par(i, 16, _): par(i, 16,!), par(i, 1, _);


process = hit:model<:_,_;
