clear all; 
close all;

fs = 44100;         % sample rate
k = 1 / fs;         % time step
dur = fs*2;           % duration

% Define variables of the system
c = 344;

% Calculate grid spacing from variables
h = c * k*sqrt(2);
Nx = 5;%floor(1/h);
Ny = 5;

r=50;
%h = 1 / N;

lambda = c*k/h;

outPos = 3;

uNext = zeros(Nx,Ny);
u = zeros(Nx,Ny);
width = 10;
%u((1:width) + floor(Nx/2)-(width/2),(1:width) + floor(Ny/2)-(width/2)) = window2(width,width,'hann');
%u(floor(Nx/2),floor(Ny/2)) = 1;
uPrev = u;
out = zeros(dur,1);

for n=1:dur    
%     u(floor(Nx/2),floor(Ny/2)) = sin(2*pi*440*n/fs)/3;

    if n==2
        %u((1:width) + floor(Nx/2)-(width/2),(1:width) + floor(Ny/2)-(width/2)) = window2(width,width,'hann');
        u(outPos,outPos)=1;
    end
    
    %osservo solo UN punto della "corda"
    %out(n)= uNext(outPos,outPos);
    out(n)= u(outPos,outPos);
    
    %Full Dirichlet
    for l=2:Nx-1
        for m = 2:Ny-1
            uNext(l,m) = 2*(1-2*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1));
        end
    end
    
%     %Full Neumann
%     for l=1:Nx
%         if l==1
%             uNext(l,1) = 2*(1-lambda^2)*u(l,1) - uPrev(l,1) + lambda^2*(u(l+1,1) + u(l,2));
%             for m = 2:(Ny-1) %taking corners into account
%                 uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l,m+1) + u(l,m-1));
%             end
%             uNext(l,Ny) = 2*(1-lambda^2)*u(l,Ny) - uPrev(l,Ny) + lambda^2*(u(l+1,Ny) + u(l,Ny-1));
%         elseif l==Nx
%             uNext(l,1) = 2*(1-lambda^2)*u(l,1) - uPrev(l,1) + lambda^2*(u(l-1,1) + u(l,2));
%             for m = 2:(Ny-1) %taking corners into account
%                 uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l-1,m) + u(l,m+1) + u(l,m-1));
%             end
%             uNext(l,Ny) = 2*(1-lambda^2)*u(l,Ny) - uPrev(l,Ny) + lambda^2*(u(l-1,Ny) + u(l,Ny-1));
%         else
%             for m = 1:Ny
%                 if m==1
%                     uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m+1));
%                 elseif m==Ny
%                     uNext(l,m) = (2-3*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m-1));
%                 else
%                     uNext(l,m) = 2*(1-2*lambda^2)*u(l,m) - uPrev(l,m) + lambda^2*(u(l+1,m) + u(l-1,m) + u(l,m+1) + u(l,m-1));
%                 end
%             end
%         end
%     end    
   
    
%     surf(uNext)
%     %plot(u);
%     zlim([-1,1]);
%     view([45 45]);
%     %xlim([50,150]);
%     drawnow;
%     pause(50/1000)
    
    uPrev = u;
    u = uNext;
end

soundsc(out,fs)


figure(1)
plot(out)