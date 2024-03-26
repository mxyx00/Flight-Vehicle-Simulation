








%% External Functions

function [temperature, density, viscosity, a] =  alt(height)

altitude = [0:500:25000]; % creates altitude points from 0-25 km (in m)

% Finding the index of where the troposphere/stratosphere begins/ends
troposphere = find(altitude==11000);
stratosphere = find(altitude==25000);

% Sea Level coefficients
T1 =288.16; % K
p1= 1.01325*10^5; % N/m^2
rho1 = 1.225; %kg/m^2
g0 = 9.80; % m/s^2
R = 287; % J/kgK


% Troposphere Gradient
a = ((216.66-288)/11000);

% For loop that runs through all the altitudes, and determines region to
% calculate through applicable formulas.

if height <= 11000
    temperature = height*a+T1;% altitude times the temperature gradient plus temp at sea level (288K)
    density = rho1*((temperature/T1)^-((g0/(a*R))+1));
elseif height > 11000
    temperature = temperature(troposphere);
    density = density(troposphere)*exp(-(g0/(R*temperature))*(altitude-altitude(troposphere)));
end
viscosity = 1.54*(1+0.0039*(temperature-250))*10^-5;
a = sqrt(1.4*R*temperature);



