[temperature, density, viscosity, a] = alt(2000)
T1 =288.16; % K
p1= 1.01325*10^5; % N/m^2
rho1 = 1.225; %kg/m^2
g0 = 9.80; % m/s^2
R = 287;

% Variables
S = 554
b = 68.45
AR=b^2/S;
WTO = 4390000
e = 0.82
CD0 = 0.015
k=1/(pi*e*AR);



TAsl = 100000

TA= TAsl*(density/rho1)



%% Minimum and Maximum Flight Speeds
qmax=((TA)/(2*S*CD0))*( 1  +  sqrt( 1 - (4*k*CD0*WTO^2)/(TA^2)  ));
qmin=((TA)/(2*S*CD0))*( 1  -  sqrt( 1 - (4*k*CD0*WTO^2)/(TA^2)  ));
vmax=sqrt(2*qmax/rho1)
vmin=sqrt(2*qmin/rho1)



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

end

%% External Functions for HW

% Modified function from previous homework
function [rate]= rateofClimb(k, v, rho, Wo, S, CD0, e, AR, Tmax, condition, J, tc, a)

    % Determining whether compressible effects included or not. 
    for i = 1:length(v)
        dynamic(i)=0.5*rho*v(i)^2;
        CL(i)=Wo/(S*dynamic(i));
        CD(i)=CD0+(k.*(CL(i)).^2);

        if condition == 0
            D(i)=CD(i)*dynamic(i)*S;
        else
            M(i)=v(i)/a;
            Mcc1(i)=0.87-(0.175*CL(i))-(0.83*tc);
            m(i)=0.83-0.583*CL(i)+0.111*CL(i)^2;
            Mcc2(i)=Mcc1(i)/((cosd(J))^m(i));

            x(i)=M(i)/Mcc2(i);

            if CL(i) >1.4
                dCDc(i)= 0;
            else
                dCDc(i)=(cosd(J))^3*(3.97*(10^-9)*exp(12.7*x(i))+(10^-40)*exp(81*x(i)));
            end

            D(i)=(CD(i)+dCDc(i))*dynamic(i)*S;
        end
    end

    rho1 = 1.225;

    for i = 1:length(v)
        rate(i)=v(i)*(Tmax*(rho/rho1)-D(i))/Wo;
    end

end
