import("stdfaust.lib");
inPoint=hslider("input point", floor(nPoints/2),0,nPoints-1,1);
outPoint=hslider("output point",floor(nPoints/2),0,nPoints-1,1);

nPoints=3;
force = os.osc(440):ba.selectoutn(nPoints,inPoint);//1 <: par(i, nPoints, *(i==inPoint));
process = force:ba.selectn(nPoints,outPoint);
