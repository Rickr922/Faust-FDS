%{
    stiff String damp
%}

clear all; 
close all;

fs = 44100;         % sample rate
k = 1 / fs;         % time step
dur = 5*fs;           % duration

% Define variables of the system
T = 150;                 % Tension [N]
radius = 3.5560e-04;    % Radius (0.016 gauge) [m]
rho = 8.05*10^3;        % Density [kg/m^3];
Area = pi*radius^2;        % Area of string section
I = (pi*radius^4)/ 4;   % Moment of Inertia
L = 1;                     % String length [m]
Emod = 174e4;              % Young modulus [Pa]
K = sqrt(Emod*I/rho/Area);    % Stiffness parameter
c = sqrt(T/rho/Area);      % Wave speed
sigma1 = 0.01;         % Frequency dependent damping
sigma0 = 0.0005;

% Calculate grid spacing from variables
A0 = c^2*k^2 + 4*sigma1*k;
h =sqrt((A0 + sqrt((A0)^2 + 16*k^2*K^2))/2);
%N = floor(L/h);
N = 92;
%h = L / N;

FB = 1000000; % bow force/mass (m/s^2)
vB = 0; % bow velocity (m/s)
a=100;

outPos = floor(N/2);
%outPos = 103;

% Initialise u^{n+1} and u^n
uNext = zeros(N, 1);
u = zeros(N, 1);

% Excitation using hann window
% width = floor(N / 5);
% excitationRange = 1:width;
% u(excitationRange + floor(N / 2)) = hann(width);

%u(outPos)=1;

% initialise u^{n-1}
uPrev = u;

% Initialise output and output position
out = zeros(dur, 1);

% A = 2/(1+k*sigma0);
% B = (k*sigma0 - 1)/(k*sigma0 + 1);
% C = (c^2*k^2 + 2*sigma1*k)/((h^2)*(1+k*sigma0));
% D = 2*sigma1*k/((h^2)*(1+k*sigma0));
% E = (K^2*k^2)/((h^4)*(1+k*sigma0));

den = 1+sigma0*k;
A = (2*h^4-2*c^2*k^2*h^2-4*sigma1*k*h^2+6*K^2*k^2)/den/h^4;
B = (sigma0*k*h^2-h^2+4*sigma1*k)/den/h^2;
C = (c^2*k^2*h^2+2*sigma1*k*h^2-4*K^2*k^2)/den/h^4;
D = -2*sigma1*k/den/h^2;
E = K^2*k^2/den/h^4;

%% Loop
for n = 1:dur
    
     if n==10
         vB=0.2;
     end
    dVel = (u(outPos)-uPrev(outPos))/k - vB;
    jCoeff = k^2*FB/den/h;
    f = -jCoeff*1.41*a*dVel*exp(-a*dVel^2+0.5);
        
    %nested for-loop
    for l = 3:N-2
%         d4 = u(l+2) - 4*u(l+2) + 6*u(l) - 4*u(l-1) + u(l-2);
%         uNext(l) = A*u(l) + B*uPrev(l) + C*(u(l+1) - 2 * u(l) + u(l-1)) - D*(uPrev(l+1) - 2*uPrev(l) + uPrev(l-1)) - E*d4;
        uNext(l) = u(l)*A + uPrev(l)*B + (u(l+1)+u(l-1))*C + (uPrev(l+1)+uPrev(l-1))*D + (u(l+2)+u(l-2))*E + f;
    end
    
    % retrieve output
    out(n) = uNext(outPos);
    out2(n) = dVel;
    out3(n)=f;
%     %draw string
%     %plot(uNext);
%     plot(u);
%     ylim([-1, 1]);
%     drawnow;
%     pause;
    
    uPrev = u;
    u = uNext;
end

soundsc(out,fs)
% plot(out)
plot(out3*1000000)
hold on
plot(out2)