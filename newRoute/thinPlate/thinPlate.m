clear all; 
close all;

fs = 44100;         % sample rate
k = 1 / fs;         % time step
dur = fs*2;           % duration

K = 20;
s0 = 2;
s1 = 0.05;
c=344;

coeff = c^2*k^2+4*s1*k;
h = sqrt((coeff+sqrt(coeff*coeff+16*K*K*k*k)));
%h=2*sqrt(k*K);
Nx = 70;
Ny = 70;

lambda = c*k/h;

outPos = 5;

uNext = zeros(Nx,Ny);
u = zeros(Nx,Ny);
width = 10;
% u((1:width) + floor(Nx/2)-(width/2),(1:width) + floor(Ny/2)-(width/2)) = window2(width,width,'hann');
%u(floor(Nx/2),floor(Ny/2)) = 1;
%u(5,5) = 1;
uPrev = u;
out = zeros(dur,1);

mu=K*K*k*k/(h^4);
den = 1+s0*k;
A = 2*(1-10*mu-2*lambda*lambda-4*s1*k*k)/den;
B = (s0*k+4*k*k-1)/den;
C = (8*mu + lambda*lambda + 2*s1*k*k)/den;
D = -2*mu/den;
E = -mu/den;
F = -2*s1*k*k/den;

% A=2*(1-10*mu);
% B=-1;
% C=8*mu;
% D=-2*mu;
% E=-mu;

for n=1:dur    
%     u(floor(Nx/2),floor(Ny/2)) = sin(2*pi*440*n/fs)/3;
    fin=0;
    if n==2
        %u((1:width) + floor(Nx/2)-(width/2),(1:width) + floor(Ny/2)-(width/2)) = window2(width,width,'hann');
        fin=1;
    end
    
    %osservo solo UN punto della "corda"
    %out(n)= uNext(outPos,outPos);
    out(n)= u(outPos,outPos);  
    
    for l=3:Nx-2
        for m = 3:Ny-2
            uNext(l,m) = A*u(l,m) + B*uPrev(l,m) + C*(u(l+1,m)+u(l-1,m)+u(l,m+1)+u(l,m-1)) + D*(u(l+1,m+1)+u(l-1,m+1)+u(l+1,m-1)+u(l-1,m-1)) + E*(u(l+2,m)+u(l-2,m)+u(l,m+2)+u(l,m-2)) + F*(uPrev(l+1,m)+uPrev(l-1,m)+uPrev(l,m+1)+uPrev(l,m-1))+fin;
        end
    end  

%  for l=3:Nx-2
%         for m = 3:Ny-2
%             uNext(l,m) = A*u(l,m) + B*uPrev(l,m) + C*(u(l+1,m)+u(l-1,m)+u(l,m+1)+u(l,m-1)) + D*(u(l+1,m+1)+u(l-1,m+1)+u(l+1,m-1)+u(l-1,m-1)) + E*(u(l+2,m)+u(l-2,m)+u(l,m+2)+u(l,m-2));
%         end
%     end 
    
%     surf(uNext)
%     %plot(u);
%     zlim([-1,1]);
%     view([45 45]);
%     %xlim([50,150]);
%     drawnow;
%     %pause(50/1000)
    
    uPrev = u;
    u = uNext;
end

soundsc(out,fs)


figure(1)
plot(out)