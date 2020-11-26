%{
    Clavinet
%}

clear all; 
close all;
clc;

fs = 44100; % sample rate
k = 1 / fs;  % time step
dur = 5*fs; % duration

%Trial String
T = 150;                 % Tension [N] 
radius = 3.5560e-04;    % Radius (0.016 gauge) [m] 

%
%   SETUP STRING PHYSICAL VALUES
%

rho = 8.05*10^3;        % Density [kg/m^3];
A = pi*radius^2;        % Area of string section
I = (pi*radius^4)/ 4;   % Moment of Inertia
L = 1;                  % String length [m]
E = 174e6;              % Young modulus [Pa]
K = sqrt(E*I/rho/A);    % Stiffness parameter
c = sqrt(T/rho/A);      % Wave speed
sigma1 = 0.01;         % Frequency dependent damping

coeff = c^2*k^2 + 4*sigma1*k;
h =sqrt((coeff + sqrt(coeff^2+((16*E*I*k^2)/(rho*A))))/2);
%N = floor(L/h);
N=50;
sigma0 = 0.0005;

%
%   GRID FUNCTION INITIALIZATION
%

uNext = zeros(N, 1);
u = zeros(N, 1);
uPrev = u;

%hammer
offset = -0.05;
uHNext = 0;
uH=0;
uHPrev=0;
uHNextOut = 0;

KH = 500;
mH = 0.1;
omega0Sqr = KH/mH;
sigma0H = 30;
FIN = 0.02;

% % raised cosine
% tSilence1 = 0;
% tAttack = 1*fs/1000;
% tRelease = 1*fs/1000;
% tSilence2 = 0;
% tSustain = dur-tSilence1-tSilence2-tAttack-tRelease;
% t1 = tSilence1+tAttack;
% t2 = tSilence1+tAttack+tSustain;
% t3 = tSilence1+tAttack+tSustain+tRelease;
% force = 1;
% fInTimeProfile = zeros(dur, 1);
% for t = 1:dur
%     if t>tSilence1 && t<t1
%         fInTimeProfile(t) = 0.5*(1-cos(pi*(t-tSilence1)/tAttack))*force;
%     elseif t1<=t && t<=t2
%         fInTimeProfile(t) = force;
%     elseif t>t2 && t<=t3
%         fInTimeProfile(t) = 0.5*(1-cos(pi*(t-t3)/tRelease))*force;
%     end
% end

% %dirac delta interpolator
% excitPos = 0.5;% (L-excitLength)/h;
% % Dirac Delta approx & J coefficient
% excitPosFloor = floor(excitPos);
% beta = excitPos - excitPosFloor;
% delta = zeros(N,1);
% for i = 1:N
%     if i == excitPosFloor
%         delta(i) = (1-beta)/h;
%     elseif i == excitPosFloor+1
%         delta(i) = beta/h;
%     end
% end

JCoeff = zeros(N,1);
% for i = 1:N
%     JCoeff(i) = k^2*delta(i)/(1-sigma0*k)/(rho*A);
% end

inPos = floor(N/2);
JCoeff(inPos) = k^2/h/(1-sigma0*k)/(rho*A);

% stud force values
epsilon = 0.05;        %distance between stud and string [m]
kStud = 10000;
alpha =3;

%
%   OUTPUT POSITIONS INITALIZATION
%

out1 = zeros(dur, 1);
outPos1 = floor(N/2);

%   OUTPUT VECTOR INITIALIZATION

out = zeros(dur, 1);

%   STRING COEFFICIENT INITIALIZATION

for i = 1:N
      coeff0(i) = 2/(1+k*sigma0);
      coeff1(i) = (k*sigma0 - 1)/(k*sigma0 + 1);
      coeff2(i) = (c^2*k^2 + 2*sigma1*k)/((h^2)*(1+k*sigma0));
      coeff3(i) = 2*sigma1*k/((h^2)*(1+k*sigma0));
      coeff4(i) = (K^2*k^2)/((h^4)*(1+k*sigma0));
end

%
%   SIMULATION
%
extF=0;
for n = 1:dur
    
    %studParam = -(1-beta)*u(excitPosFloor) - beta*u(excitPosFloor+1) - epsilon;
    studParam = uHNextOut-u(inPos);
    if studParam > 0
        fStud = kStud*(studParam)^alpha;
    else
        fStud = 0;
    end
    
    if n==5 extF = FIN;
    else extF = 0; end
    
    uHNext = uH*((2-omega0Sqr^2*k^2)/(1+sigma0H*k))-uHPrev*((1-sigma0H*k)/(1+sigma0H*k)) +extF -  fStud*k^2/(1+sigma0H*k);%+ extF*k^2/(1+sigma0H*k);
    uHNextOut = uHNext+offset;
    
    for l = 3:N-2 
        d4 = u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2);
        uNext(l) = coeff0(l)*u(l) + coeff1(l)*uPrev(l) + coeff2(l)*(u(l+1) - 2 * u(l) + u(l-1)) - coeff3(l)*(uPrev(l+1) - 2*uPrev(l) + uPrev(l-1)) - coeff4(l)*d4 + JCoeff(l)*fStud;%JCoeff(l)*(fStud - fInTimeProfile(n));
    end
    
    % retrieve outputs
    %out(n) = uNext(outPos1);
    out(n) = uHNextOut;
    %out(n) = fStud;
    
%     % % Debug: draw string
%     plot(uNext);
%     %ylim([-1, 1]);
%     drawnow;
%     
    % update grid function
    uPrev = u;
    u = uNext;
    
    uHPrev = uH;
    uH = uHNext;
end
plot(out)
%ylim([-1.5,1]);

%soundsc(out,fs)