%{
    Clavinet
%}

clear all; 
close all;
clc;

fs = 44100; % sample rate
k = 1 / fs;  % time step
dur = 5*fs; % duration

%Set Output Frequency [Hz]
requestedFrequency = 150;

%   STRING SET 
%   T = Tension [N]; gauge [inch]->radius [m]
%   Tension calculator:
%   www.evertune.com/resources/string_tension_gauge_calculator.php
%   String gauges:
%   www.wurlitzer-epservice.nl/clavinet-single-strings

if requestedFrequency <= 75 %String 1 F1-D2
    %T = 25.31;
    T = 100;
    radius = 8.128e-04; % 0.032 gauge
elseif requestedFrequency <= 100 %String 2 D#2 - G2
    T = 60.85;
    radius = 7.112e-04; % 0.028 gauge
elseif requestedFrequency <= 170 %String 3 G#2 - E3 
    T = 67.84;
    radius = 5.588e-04; % 0.022 gauge
elseif requestedFrequency > 170 %String 4 F3 - E6
    T = 37.59;
    radius = 2.286e-04; % gauge 0.009
end

% %Trial String
% T = 75;                 % Tension [N] 
% radius = 3.5560e-04;    % Radius (0.016 gauge) [m] 

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
sigma1 = 0.005;         % Frequency dependent damping

%
%   STRING STABILITY CONDITION
%   (slightly (1.2) greater than the limit)

coeff = c^2*k^2 + 4*sigma1*k;
h =sqrt((coeff + sqrt(coeff^2+((16*E*I*k^2)/(rho*A))))/2)*1.2;
N = floor(L/h);

%
%   FREQUENCY INDEPENDENT DAMPING
%

yarnPos = 10; % percentage of the string covered by the yarn
yarnLimit = yarnPos/100; %lower limit for excitation position
sigma0 = zeros(N, 1);
sigma0(1:floor(N*yarnPos/100)) = 500;
sigma0(floor(N*yarnPos/100):N) = 0.005;

%
%   GRID FUNCTION INITIALIZATION
%

uNext = zeros(N, 1);
u = zeros(N, 1);
uPrev = u;

%
% EXCITATION
%

% force input
tSilence1 = 0;
tAttack = 1*fs/1000;
tRelease = 1*fs/1000;
tSilence2 = 0;
tSustain = dur-tSilence1-tSilence2-tAttack-tRelease;
t1 = tSilence1+tAttack;
t2 = tSilence1+tAttack+tSustain;
t3 = tSilence1+tAttack+tSustain+tRelease;
force = 1;
fInTimeProfile = zeros(dur, 1);
for t = 1:dur
    if t>tSilence1 && t<t1
        fInTimeProfile(t) = 0.5*(1-cos(pi*(t-tSilence1)/tAttack))*force;
    elseif t1<=t && t<=t2
        fInTimeProfile(t) = force;
    elseif t>t2 && t<=t3
        fInTimeProfile(t) = 0.5*(1-cos(pi*(t-t3)/tRelease))*force;
    end
end
% figure(1)
% plot(fInTimeProfile);
%  ylim([0, 2]);
% xline(tAttack,'--r');

excitLength = c/2/requestedFrequency;
%excitLength = 0.8;
if excitLength < yarnLimit
    error('string length incorrect');
end
excitPos = (L-excitLength)/h;

% Dirac Delta approx & J coefficient
excitPosFloor = floor(excitPos);
beta = excitPos - excitPosFloor;
delta = zeros(N,1);
for i = 1:N
    if i == excitPosFloor
        delta(i) = (1-beta)/h;
    elseif i == excitPosFloor+1
        delta(i) = beta/h;
    end
end
JCoeff = zeros(N,1);
for i = 1:N
    JCoeff(i) = k^2*delta(i)/(1-sigma0(i)*k)/(rho*A);
end

% stud force values
epsilon = 0.001;        %distance between stud and string [m]
kStud = 10000;
alpha =1.1;

%
%   OUTPUT POSITIONS INITALIZATION
%

out1 = zeros(dur, 1);
outPos1 = floor(85*N/100);
out2 = zeros(dur, 1); %on the yarn
outPos2 = floor((yarnPos/2)*N/100);

%
%   OUTPUT VECTOR INITIALIZATION
%

out = zeros(dur, 1);

%
%   STRING COEFFICIENT INITIALIZATION
%

for i = 1:N
      coeff0(i) = 2/(1+k*sigma0(i));
      coeff1(i) = (k*sigma0(i) - 1)/(k*sigma0(i) + 1);
      coeff2(i) = (c^2*k^2 + 2*sigma1*k)/((h^2)*(1+k*sigma0(i)));
      coeff3(i) = 2*sigma1*k/((h^2)*(1+k*sigma0(i)));
      coeff4(i) = (K^2*k^2)/((h^4)*(1+k*sigma0(i)));
end

%
%   SIMULATION
%
timeValue=0;
for n = 1:dur
    %studParam = -(1-beta)*u(excitPosFloor) - beta*u(excitPosFloor+1) - epsilon;
     if n == 100 timeValue=1;
    end
    studParam = timeValue*epsilon-(1-beta)*u(excitPosFloor) - beta*u(excitPosFloor+1);
    if studParam >= 0
        fStud = kStud*(studParam)^alpha;
    else
        fStud = 0;
    end
    for l = 3:N-2 
        d4 = u(l+2) - 4*u(l+1) + 6*u(l) - 4*u(l-1) + u(l-2);
        uNext(l) = coeff0(l)*u(l) + coeff1(l)*uPrev(l) + coeff2(l)*(u(l+1) - 2 * u(l) + u(l-1)) - coeff3(l)*(uPrev(l+1) - 2*uPrev(l) + uPrev(l-1)) - coeff4(l)*d4 + JCoeff(l)*fStud;%JCoeff(l)*(fStud - fInTimeProfile(n));
    end
    
    % retrieve outputs
    out1(n) = uNext(outPos1);
    out2(n) = uNext(outPos2);
    out(n) = (out1(n) + out2(n))/2;
    
%     % % Debug: draw string
%     plot(uNext);
%     %ylim([-1, 1]);
%     drawnow;
%     
    % update grid function
    uPrev = u;
    u = uNext;
end

% Normalizing output
% maxOut = max(out);    % find max value of output
% minOut = abs(min(out));
% if maxOut>minOut
%     out = out/maxOut;
% else
%     out = out/minOut;
% end

%
%   PICKUP
%

% out = out - mean(out); % center output on 0
% 
% limiterValue = 0.1; % this value must be between 0 and 1
% for n = (1 : dur)
%     if abs(out(n)) > limiterValue
%         out(n) = limiterValue * sign(out(n));
%     end
% end

%
%   RESULT
%

plot(out)
soundsc(out,fs)
%audiowrite("120HzClavinet.wav",out,fs);