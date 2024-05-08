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
C_D_0 = 0.0145;
k=1/(pi*e*AR);
throttle = 1;
TAsl = throttle*(4*419.25)*1000; %ge90b
J = 20;
tc = 0.12;
C_L_max = 2.2;
CL = 1.1;
cj = 0.520;
wingtipheight = 5.4;
g=9.81;

W_empty = 888328.20;
W_fuel = 2673660.00;
max_payload = 1912735.29; % calculated to N already
min_payload = 1312225;        

WTO = W_empty+W_fuel; 

%% Weight/Thrust Trade Study
% 
% weights = [WTO+4*81225,WTO+4*81225,WTO+4*35586,WTO+4*31582];
% thrust =  [4*419.25*1000,4*366.75*1000,4*193*1000,4*194.54*1000];
% heights = 0:500:13000;
% colors = ['b','r','g','m']
% 
% for j = 1:length(weights)
%     for i = 1:length(heights)
%         [temp, rho, visc, a] = alt(heights(i));
%         TA(i)= thrust(j)*(rho/rho1);
% 
% %         cldmin = sqrt(pi*e*AR*CD0);
% %         vdmin(i) = sqrt(2*weights(j)/rho*S*C_L_max);
% 
%         qmax=((TA(i))/(2*S*CD0))*( 1  +  sqrt( 1 - (4*k*CD0)/(TA(i)/weights(j))^2));
%         qmin=((TA(i))/(2*S*CD0))*( 1  -  sqrt( 1 - (4*k*CD0)/((TA(i)/weights(j))^2)));
%         vmax(i)=sqrt(2*qmax/rho);
%         vmin(i)=sqrt(2*qmin/rho);
% 
% 
%     end
%     plot(vmax,heights,colors(j))
%     hold on
% %     plot(vmin,heights,colors(j))
% 
% end
% 
% xlabel("Velocity (m/s)")
% ylabel("Altitude (m)")
% yyaxis right
% ylabel("Altitude (ft)")
% ftlimit = max(heights)*3.28084;
% ylim([ 0,ftlimit])
% yticks([0:2500:ftlimit]);
% legend("GE90-90b","GE90-77B","GE-TF39","PW-2000")
% set(gca, 'FontName', 'Times')

%% Wingspan Trade Study

% V = 150;
% [temp, rho, visc, a] = alt(0);
% q = (0.5)*rho*V^2;
% b_ranges = linspace(40,80,30);
% S_ranges = linspace(10,700,30);
% 
% for i = 1:30
%     AR = b^2/S_ranges(i);
%     Dp(i) = C_D_0*q*S_ranges(i);
%     Di(i) = (WTO)^2/(pi*e*AR*q*S_ranges(i));
%     Dt(i) = Dp(i)+Di(i);
% end
% 
% plot(S_ranges, Dp)
% hold on
% plot(S_ranges, Di)
% plot(S_ranges, Dt)
% ylabel("Drag (N)")
% xlabel("Wing area(m^2)")
% legend("Parasitic Drag", "Induced Drag", "Total Drag")
% set(gca, 'FontName', 'Times')


%% Minimum and Maximum Flight Speeds
% For a range of heights, calculating the density, dynamic pressure,
% maximum/minimum velocity, and calculating the vstall. 
figure()
altitude = 0;
altitude_msg = [" at Altitude = " + altitude " m"];
height_range = 0:100:altitude+10000;
vstall = 0:length(height_range)-1;

for i = 1:length(height_range)

    [temperature, rho, viscosity, a] = alt(height_range(i));
    TA= TAsl*(rho/rho1);
    Tmax = TA;

    qmax=((TA)/(2*S*C_D_0))*( 1  +  sqrt( 1 - (4*k*C_D_0)/(TA/WTO)^2));
    qmin=((TA)/(2*S*C_D_0))*( 1  -  sqrt( 1 - (4*k*C_D_0)/((TA/WTO)^2)));
    vmax(i)=sqrt(2*qmax/rho);
    vmin(i)=sqrt(2*qmin/rho);

%     if vstall(i) < vmin(i)
        vstall(i) = sqrt((2*WTO/(rho*S*C_L_max)));
%     end

    vsc(i) = sqrt((2*WTO)/(rho*S*sqrt(pi*e*AR*C_D_0)));
    em = 1/(sqrt(4*k*C_D_0));
    qfc = ((TA)/(6*S*C_D_0))*(1+sqrt(1+(3/(em^2*(TA/WTO)^2))));
    vfc(i) = sqrt(2*qfc/rho);
    vec(i) = (vfc(i)+vsc(i))/2;
end

vstall_height = vstall(1);

ceiling = (WTO*(4*k*C_D_0)^0.5)/(TAsl)*rho1; % 


plot(vmax, height_range)
hold on
plot(vmin, height_range)
plot(vstall, height_range)
plot(vsc, height_range,'--')
plot(vfc, height_range,'--')
plot(vec, height_range,':')
xlabel("Velocity (m/s)")
ylabel("Altitude (m)")
yyaxis right
ylabel("Altitude (ft)")
ftlimit = max(height_range)*3.28084;
ylim([ 0,ftlimit])
yticks(0:2500:ftlimit);
msg = ["Maximum Velocity (m/s)", "Minimum Velocity(m/s)", "Stall Velocity (m/s)", "Velocity of Steepest Climb (m/s)", "Velocity of Fastest Climb (m/s)", "Velocity of Most Economical Climb (m/s)"];
legend(msg)
legend('Location','southoutside')
set(gca, 'FontName', 'Times')



%% Thrust vs Velocity with Ranges
% Different values of area, looped through different values of v. 
% Altitude constant
S_ranges = S;
b_ranges = linspace(b-20,b+2,5);
AR_ranges=  (b_ranges.^2)./S_ranges;
k_ranges= 1./(pi*e.*AR_ranges);
[temperature, rho, viscosity, a] =  alt(altitude);
Tmax = TAsl*(rho/rho1);
L = WTO;

figure()
for j = 1:length(b_ranges)

    vstall = sqrt((2*WTO/(rho*S_ranges*C_L_max)));
    vrange = vstall:2:vmax(1);
    
    for i = 1:length(vrange)    
        q = 0.5*rho*vrange(i)^2;
        TR(i) = q*S_ranges*C_D_0 + (k_ranges(j)*WTO^2)/(q*S_ranges);
    end
    plot(vrange,TR)
    hold on  
end

yline(Tmax)
xlabel("Velocity (m/s)")
ylabel("Thrust Required (N)")
% msg = ["Thrust (N) with S = " + S_ranges(1) + " m^2", "Vmin = " + vdmin + "m/s at TR = " + tratmin + " N", "TA = " + Tmax];
msg = ["Thrust (N) with b = " + b_ranges(1) + " m", "Thrust (N) with b = " + b_ranges(2) + " m","Thrust (N) with b = " + b_ranges(3) + " m", "Thrust (N) with b = " + b_ranges(4) + " m", "Thrust (N) with b = " + b_ranges(5) + " m","Tmax = " + Tmax/100 + "kN"];
legend(msg)
set(gca, 'FontName', 'Times')


%% Velocity at minimum drag
figure()
alt_ranges = [altitude, altitude+5000];
vrange = vstall_height:2:vmax(1);
colors = ["blue","red"];

for j = 1:length(alt_ranges)

    [temperature, rho, viscosity, a] =  alt(alt_ranges(j));

    for i = 1:length(vrange)
        q = 0.5*rho*vrange(i)^2;
        dminp(i) = 2*C_D_0*q*S;
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
set(gca, 'FontName', 'Times')

%% Rate of Climb + Velocity at Climb
figure()
alt_range = [0:4500:9000];
rate_range = [100,40,3];

for j = 1:length(alt_range)
    subplot(3,1,j)
    [temperature, rho, viscosity, a] =  alt(alt_range(j));
    TA= TAsl*(rho/rho1);
    for i = 1:length(vrange)
        % With compressible effects
        rate1_comp = rateofClimb(k, vrange, rho, WTO, S, C_D_0, e, AR, TA, 1, J, tc, a);
        % Without
        rate1_incomp = rateofClimb(k, vrange, rho, WTO, S, C_D_0, e, AR, TA, 0);
    end
    plot(vrange, rate1_incomp)
    hold on
    plot(vrange,rate1_comp, '--')
    xlabel("Velocity (m/s)")
    ylabel("Rate of Climb (m/s)")
    ylim([0 rate_range(j)])
    subtitle("Altitude = " + alt_range(j))
end

vsc = sqrt((2*WTO)/(rho*S*sqrt(pi*e*AR*C_D_0)));

legend(" Incompressible", " Compressible")

%% Range
figure()

W0 = W_empty+W_fuel;
W1 = W0-W_fuel;
W1_reserves = W0-0.85*W_fuel;

Em = 1/sqrt(4*k*C_D_0);
Sfuel = W_fuel/W0;
height_range = 0:100:15000;
weightrange = [W0,W0+min_payload,W0+max_payload]

% For constant height and CL
% Variable velocity 
for j = 1:length(weightrange)
    for i = 1:length(height_range)
        W1_reserves = weightrange(j)-0.85*W_fuel;
        [temperature, density, viscosity, a] =  alt(10000);
        CL = (weightrange(j))/(q*S);
        v = sqrt((2*weightrange(j))/rho*S*CL);
        q = 0.5*rho*vrange(i)^2;
        CD = C_D_0 + k*CL^2;
        E = CL/CD;
        Sj = (vrange(i)*E/cj*weightrange(j))*3.6;%in KM
        
        xhcl(i) = (2/cj)*sqrt(2/(rho*S))*sqrt(CL)/CD*(sqrt(weightrange(j))-sqrt(W1_reserves));
    end

    subplot(3,1,j)
    plot(height_range,xhcl)
    xlabel("Altitude (m)")
    ylabel("Range (km)")
    set(gca, 'FontName', 'Times')
end



% subplot(2,1,2) 
% plot(hrange,xvcl)
% title("Range of HLA at different heights")
% xlabel("Height (m)")
% ylabel("Range (km)")
% set(gca, 'FontName', 'Times')



%% Ground Effect Comparison
figure()
% wingtipheight_range = 0:0.1:17;
b_range = 5:0.1:78;

for i = 1:length(b_range)
    hbratio(i) = wingtipheight/b_range(i);
    phi(i) = (16*hbratio(i))^2/(1+(16*hbratio(i))^2);
end

plot(hbratio,phi)
title("Ground effect vs h/b ratio")
subtitle("Wingspan held constant, height varied")
xlabel("H/b ratio")
ylabel("Ground effect")
set(gca, 'FontName', 'Times')

%% Takeoff Performance
figure()
phi = 0.9;
vlo = 1.2*vstall;
q = 0.5*rho*vlo^2;
Dp = q*S*C_D_0;
Di = phi*k*CL^2*q*S;
uf = [0.02,0.07,0.12]; % for slightly unsmooth surface
L = CL*q*S;
Dphi = Dp+Di;

throttle = 0.8:0.02:1;

for j = 1:length(uf)
    subplot(3,1,j)
    for i = 1:length(throttle)
        TAsl(i) = throttle(i)*4*413*1000;
        
        dlo(i) = (1.44*WTO^2)/(g*rho*S*C_L_max*(TAsl(i)-(Dphi + uf(j)*(WTO-L))));
        chi = WTO^2/(S*C_L_max*TAsl(i));
        dto(i) = 0.217*chi+183;
    end
    plot(TAsl,dlo)
    hold on
    plot(TAsl,dto)
    subtitle("U_f = " + uf(j))
    xlabel("Tmax (N)")
    ylabel("Liftoff Distance (m)")
    set(gca, 'FontName', 'Times')

end

% title("Throttle varying between 80-100%")

legend("Liftoff Distance", "Takeoff Distance")

%% Turning Flight
figure()
[temperature, rho, viscosity, a] =  alt(10000);

q = [0:100:25000];
sigma = rho/rho1;
TA_1 = sigma*TAsl;
n_1 = 1/cos(45);

r_n = (2.*q)./(rho*g*sqrt(n_1^2-1));
r_clmax = ((2.*q).*(WTO/S))./((rho*g)*((C_L_max.*q).^2-(WTO/S)^2).^0.5);
  
plot(q,r_n)
hold on
plot(q,r_clmax)
ylim([0,8000])
xlim([0,25000])
xlabel("Dynamic Pressure (N/m^2)")
ylabel("Turn Radius (m)")
legend("Turning radius (n_struct)","Turning radius (C_L_max)")
set(gca, 'FontName', 'Times')

%% Landing Performance
[temperature, rho, viscosity, a] =  alt(0);

vt = 1.3*vstall;
v_landing = 0.7*vt;
q = 0.5*rho*v_landing^2;
L = C_L_max*q*S;
D = (C_D_0+k*C_L_max^2)*q*S;

dg = 1.69*(WTO)^2/(rho*g*S*C_L_max*((D+(WTO-L))));


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


