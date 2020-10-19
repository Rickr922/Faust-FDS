clear all; 
close all;

fs = 44100;         % sample rate
fs = 48000;
k = 1 / fs;         % time step
dur = fs*2;           % duration

% Define variables of the system
c = 344;

% Calculate grid spacing from variables
h = c * k;
%N = floor(1/h);
N=100;
%h = 1 / N;

lambda = c*k/h;

%outPos = floor(N/2);
outPos = 7;

uNext = zeros(N,1);
u = zeros(N,1);
width = floor(N/5);
%u((1:width) + floor(N/2)) = hann(width);
%u(floor(N/2)) = 1;
%u(5)=1;
uPrev = u;
out = zeros(dur,1);

%Nota, faccio l'update da 2 a N-1, perchè in questo modo tengo gli estremi
%fissi
range = 2:N-1;

for n=1:dur
    %Not nested loops
    %uNext(range) = 2*u(range) - uPrev(range) + (lambda*lambda)*(u(range+1) - 2*u(range) + u(range-1));
%     
%     %Dirichlet
%     for l=2:N-1
%         uNext(l) = 2*u(l) - uPrev(l) + (lambda*lambda)*(u(l+1) - 2*u(l) + u(l-1));
%     end
    
    %Neumann
    for l=1:N
        if l==1
            uNext(l) = 2*u(l) - uPrev(l) + (lambda*lambda)*(u(l+1) - u(l));
        elseif l==N
            uNext(l) = 2*u(l) - uPrev(l) + (lambda*lambda)*(u(l-1) - u(l));
            %uNext(l)=0;
        else
            uNext(l) = 2*u(l) - uPrev(l) + (lambda*lambda)*(u(l+1) - 2*u(l) + u(l-1));
        end
    end
      
    if n==1
        u(5) = 1;
    end
    
    %osservo solo UN punto della "corda"
    %out(n)= uNext(1);
     out(n)=u(1);
     
    plot(uNext)
    %plot(u);
    %ylim([-1,1]);
    %xlim([50,150]);
    drawnow;
    pause
    
    uPrev = u;
    u = uNext;
end

%sound(out,fs)


figure(1)
plot(out)

%%
% [out,fs] = audioread("faustString.wav");
% out=out(:,1);

windowLength = 1024;
window1 = hann(windowLength);

noverlap = 200;
[signal1,freq1,time1] = spectrogram(out,window1,noverlap,windowLength, fs);
figure(2);
imagesc(time1, freq1, log10(abs(signal1)));