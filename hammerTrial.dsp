import("stdfaust.lib");

k = 1/ma.SR;

excit = button("click me"):ba.impulsify*0.006;

KH = 1000;
mH = 0.9;
omega0SqrH = KH/mH;
sigma0H = 14;
alpha = 2.5;

nlHammer(omega0Sqr,sigma0,kH,alpha,K,offset,fIn) =
    (hammerForce<:hammerModel(fIn,K,offset,_),_)~_:_,!
with
{
    hammerModel(in,K,offset) =
        (_,_,_*forceCoeff,in :> _) ~ (_ <: A*_,B*_') :_-offset;
    hammerForce(uh,u)=select2((uh-u)>0,0,((uh-u)^alpha)*(-kH));
    A = (2-omega0Sqr^2*K^2)/(1+sigma0*K);
    B = (-1)*(1-sigma0*K)/(1+sigma0*K);
    forceCoeff = K^2/(1+sigma0*K);
};

process = 0:nlHammer(omega0SqrH,sigma0H,10000,alpha,k,0.23,excit);
