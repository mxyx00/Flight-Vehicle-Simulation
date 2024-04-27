clear; clear vars; close all

T1 =288.16; % K
p1= 1.01325*10^5; % N/m^2
rho1 = 1.225; %kg/m^2
g0 = 9.80; % m/s^2
R = 287;

% Variables
S = 554;
b = 78;
AR=b^2/S;
e = 0.87;
CD0 = 0.0145;
k=1/(pi*e*AR);
TAsl = 4*413*1000;
J = 20;
tc = 0.12;
C_L_max = 2.2;
cj = 0.605/3.6;


W_empty = 1368270;
W_fuel = 1035774;
minpayload = (295000/2.205)*9.81;
maxpayload = (430000/2.205)*9.81;
% WTO = W_empty+ +W_fuel + minpayload
% WTO = 3329181;
WTO = 4392000;


%% Minimum and Maximum Flight Speeds
% For a range of heights, calculating the density, dynamic pressure,
% maximum/minimum velocity, and calculating the vstall. 
altitude = 10500;
altitude_msg = [" at Altitude = " + altitude " m"];
height_range = 0:100:altitude+10000;
vstall = 0:length(height_range)-1;

for i = 1:length(height_range)

    [temperature, rho, viscosity, a] = alt(height_range(i));
    TA= TAsl*(rho/rho1);
    Tmax = TA;

    qmax=((TA)/(2*S*CD0))*( 1  +  sqrt( 1 - (4*k*CD0)/(TA/WTO)^2));
    qmin=((TA)/(2*S*CD0))*( 1  -  sqrt( 1 - (4*k*CD0)/((TA/WTO)^2)));
    vmax(i)=sqrt(2*qmax/rho);
    vmin(i)=sqrt(2*qmin/rho);

%     if vstall(i) < vmin(i)
        vstall(i) = sqrt((2*WTO/(rho*S*C_L_max)));
%     end

    vsc(i) = sqrt((2*WTO)/(rho*S*sqrt(pi*e*AR*CD0)));
    em = 1/(sqrt(4*k*CD0));
    qfc = ((TA)/(6*S*CD0))*(1+sqrt(1+(3/(em^2*(TA/WTO)^2))));
    vfc(i) = sqrt(2*qfc/rho);
    vec(i) = (vfc(i)+vsc(i))/2;


end

vstall_height = vstall(1)

ceiling = (WTO*(4*k*CD0)^0.5)/(TAsl)*rho1; % Using appendix A this is 18,800 m

plot(vmax, height_range)
hold on
plot(vmin, height_range)
plot(vstall, height_range)
plot(vsc, height_range,'--')
plot(vfc, height_range,'--')
plot(vec, height_range,':')
title("Flight envelope for HLA")
xlabel("Velocity (m/s)")
ylabel("Altitude (m)")
yyaxis right
ylabel("Altitude (ft)")
ftlimit = max(height_range)*3.28084;
ylim([ 0,ftlimit])
yticks([0:2500:ftlimit]);
msg = ["Maximum Velocity (m/s)", "Minimum Velocity(m/s)", "Stall Velocity (m/s)", "Velocity of Steepest Climb (m/s)", "Velocity of Fastest Climb (m/s)", "Velocity of Most Economical Climb (m/s)"];
legend(msg)
legend('Location','southoutside')


%% Thrust vs Velocity with Ranges
% Different values of area, looped through different values of v. 
% Altitude constant
S_ranges = S;
AR_ranges=  (b^2)./S_ranges;
k_ranges= 1./(pi*e.*AR_ranges);
[temperature, rho, viscosity, a] =  alt(altitude);
Tmax = TAsl*(rho/rho1);
L = WTO;


figure()
for j = 1:length(S_ranges)

    vstall = sqrt((2*WTO/(rho*S_ranges(j)*C_L_max)));
    vrange = vstall:2:vmax(1);
    
    for i = 1:length(vrange)    
        q = 0.5*rho*vrange(i)^2;
        CL = L/(S_ranges(j)*q);
        CD = CD0 + k_ranges(j)*(CL)^2;
        TR(i) = q*S_ranges(j)*CD0 + (k*WTO^2)/(q*S_ranges(j));
    end

    vdmin = sqrt((2*WTO)/(rho*S_ranges(j)*sqrt(pi*e*AR*CD0)))
    tratmin = min(TR);
    plot(vrange,TR)
    hold on
    plot(vdmin,tratmin,'o')
  
end

yline(Tmax)
title("Thrust required for different velocities" + altitude_msg)
xlabel("Velocity (m/s)")
ylabel("Thrust Required (N)")
msg = ["Thrust (N) with S = " + S_ranges(1) + " m^2", "Vmin = " + vdmin + "m/s at TR = " + tratmin + " N", "TA = " + Tmax];
% msg = ["Thrust (N) with S = " + S_ranges(1) + " m^2", "Thrust (N) with S = " + S_ranges(2) + " m^2","Thrust (N) with S = " + S_ranges(3) + " m^2", "Thrust (N) with S = " + S_ranges(4) + " m^2", "Thrust (N) with S = " + S_ranges(5) + " m^2","Tmax = 826kN"];
legend(msg)


%% Velocity at minimum drag
figure()
alt_ranges = [altitude, altitude+5000];
vrange = vstall_height:2:vmax(1);
colors = ["blue","red"];

for j = 1:length(alt_ranges)

    [temperature, rho, viscosity, a] =  alt(alt_ranges(j));

    for i = 1:length(vrange)
        q = 0.5*rho*vrange(i)^2;
        dminp(i) = 2*CD0*q*S;
        dmini(i) = (WTO^2)/(q*b^2*pi*e);
    end

    plot(vrange,dminp,colors(j))

    hold on
    plot(vrange,dmini,colors(j))

end

title("Altitude Effects on Minimum Drag")
xlabel("Velocity (m/s)")
ylabel("Drag (N)")
msg = ["Min. Drag" + altitude_msg, "Min. Drag at Altitude = " + alt_ranges(2) + " m"];
legend(msg)


%% Rate of Climb + Velocity at Climb

[temperature, rho, viscosity, a] =  alt(altitude);

for i = 1:length(vrange)
    % With compressible effects
    rate1_comp = rateofClimb(k, vrange, rho, WTO, S, CD0, e, AR, TAsl, 1, J, tc, a);

    % Without
    rate1_incomp = rateofClimb(k, vrange, rho, WTO, S, CD0, e, AR, TAsl, 0);
end

vsc = sqrt((2*WTO)/(rho*S*sqrt(pi*e*AR*CD0)));

figure()
plot(vrange, rate1_incomp)
hold on
plot(vrange,rate1_comp, '--')
ylim([0 100])
xlabel("Velocity (m/s)")
ylabel("Rate of Climb (m/s)")
title("Rate of Climb vs. Velocity of Aircraft" + altitude_msg)
legend(" Incompressible", " Compressible")


%% Range
figure()
W0 = WTO;
W1 = W0-W_fuel;

Em = 1/sqrt(4*k*CD0);
Sfuel = W_fuel/W0

% For constant height and CL
% Variable velocity 
for i = 1:length(vrange)
    q = 0.5*rho*vrange(i)^2;
    CL = (WTO)/(q*S);
    CD = CD0 + k*CL^2;
    E = CL/CD;
    sj = vrange(i)*E/(cj*WTO);
    xhv(i) = sj*W_fuel;
    WD = q*S*sqrt(CD0/k);
    W0star = W0/WD;
end

% For constant velocity and CL
% Variable height
hrange = [0:100:15000];
CL = sqrt(pi*e*AR*CD0/3);
CD = (4/3)*CD0;
E = CL/CD;

for i = 1:length(hrange)
    vbestrange = sqrt((2*WTO)/(rho*S*(pi*e*AR*CD0/3)^0.5));
    xvcl(i) = ((vbestrange*E)/cj)*log(W0/W1);
end

subplot(2,1,1) 
plot(vrange,xhv)
title("Range of HLA at different velocities")
xlabel("Velocity (m/s)")
ylabel("Range (km)")

subplot(2,1,2) 
plot(hrange,xvcl)
title("Range of HLA at different heights")
xlabel("Height (m)")
ylabel("Range (km)")


% end

%% Takeoff Performance




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
temperature = height*a+T1;
density = rho1*((temperature/T1)^-((g0/(a*R))+1));

if height <= 11000
    % altitude times the temperature gradient plus temp at sea level (288K)
    density = rho1*((temperature/T1)^-((g0/(a*R))+1));
elseif height > 11000
    temperature = 11000*a+T1;
    density = density*exp(-(g0/(R*temperature))*(height-11000));
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


