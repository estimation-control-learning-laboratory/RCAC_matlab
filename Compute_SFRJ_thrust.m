function [FNET, GAMMA6] = Compute_SFRJ_thrust(PHI)
% Solid Fuel Ramjet Design
% External CEA software is used to calculate the combustion properties.
% O. Tumuklu 03/05/2022
% Additional materials would be useful: AGARD 323 report.


% Input parameters.
% FUEL TYPE set a number that is not equal to zero for the solid fuel. For this particular case, Teflon (C2F4)n is selected. 

Tfuel = 800;

%PHI = 0.99;    % Definition  of equivalance ratio  
               % PHI = (FA)/(FA)s  where (FA_s = MMfuel/(4.76MM0(X+Y/4))
               % F/A = mdot_fuel/ mdot_0
FA_s = 0.05;
A0 = 0.020325; % in m^2 == area of intake
A4 = 0.022698; % Combusiton chamber area in m^2
A5 = 0.032;    % Cross-section of C-D nozzle (i.e. choked area)
A6 = 0.048;    % Nozzle exit area in m^2
M0 = 1.0;      % Mach number of the freestream
T0 = 150;      % Free-stream temperature in K
GAMMA0 = 1.4;  % Freestream specific heat ratio.
P0 = 6700;     % Freestream pressure in Pa corresponding to a 30 km altitude where the scramjet can operate;  
ETA_C = 1.0 ; 
MM0 = 29;      % Mass of air
CD_5  = 0.86; % Nozzle design efficiency 
% GAS Constant
R_u = 8314;
% FINDING mdot_0 
mdot_0 = M0*P0*A0*sqrt(MM0*GAMMA0/(R_u*T0));
mdot_fuel = PHI*FA_s*mdot_0;
mdot_5= mdot_0 + mdot_fuel;

% Calculate inlet total temperature and pressure to estimate temperature
% and pressure before combustion. The total pressure can be smaller due to nonideal shock effects. 
Tt2 = (1+(GAMMA0-1)/2*M0^2)*T0;
Pt2 = (1+(GAMMA0-1)/2*M0^2)^(GAMMA0/(GAMMA0-1))*P0;

% GENERALLY  PT4 is measured quantity and given as an input condition; 
% BUT PT5 can be assumed to be equal to PT4. For simplicity, it is assumed
% to be equal to the total pressure.

Pt4 = Pt2/101300;% atm 

% Creating CEA input file text
fileID = fopen('cea.inp','w');
fprintf(fileID,'%s','prob case=_______________1512 hp ');
fprintf(fileID,'%s\r\n','');
fprintf(fileID,'%s','p,atm=');
fprintf(fileID,'%6.3f',Pt4);
fprintf(fileID,'%s\r\n','');
fprintf(fileID,'%s','phi=');
fprintf(fileID,'%4.2f',PHI);
fprintf(fileID,'%s\r\n','');
fprintf(fileID,'%s\r\n','');
fprintf(fileID,'%s\r\n',' reac');


fprintf(fileID,'%s',' fuel C2F4              wt%=');
fprintf(fileID,'%7.4f',100);
fprintf(fileID,'%s','  t,k= ');
fprintf(fileID,'%6.3f',Tfuel);
fprintf(fileID,'%s\r\n','');
fprintf(fileID,'%s',' oxid Air               wt%=');
fprintf(fileID,'%7.4f',100);
fprintf(fileID,'%s','  t,k= ');
fprintf(fileID,'%6.3f',Tt2);
fprintf(fileID,'%s\r\n','');


fprintf(fileID,'%s\r\n','');
fprintf(fileID,'%s\r\n',' output  short ');
fprintf(fileID,'%s\r\n',' output massf');
fprintf(fileID,'%s\r\n',' output siunits');
fprintf(fileID,'%s\r\n',' end');
fclose(fileID);

%RUN  CEA2
% cmd = fullfile(pwd,'CEA2.exe');
cmd = './CEA2.exe';
% cmd = 'cea.exe'
system(cmd);

%%

fileID = 'cea.out';
CEA_Output = fileread(fileID);

% Strings to be searched for
A = CEA_Output; 
B = strfind(A,'PHI,EQ.RATIO= ');
C = strfind(A,'T, K');
D = strfind(A,' M, (1/n)');
E = strfind(A,'GAMMAs');
F = strfind(A,'O/F');
G = strfind(A,'*O2');
H = strfind(A,'MACH NUMBER');
I = strfind(A,'P, BAR');

for i = 1:numel(B)
     phi = str2num(A(B(i)+13:B(i)+20)); 
%      disp(phi_ratio(i));
end

% Oxidizer to Fuel Ratio as reported by CEA
for i = 1:numel(F)
    O_F = str2num(A(F(i)+6:F(i)+14)); 
end

% Extracting Temps
for ii = 1:numel(C)
    TT5 = str2num(A(C(i)+16:C(i)+23));   
end

% Extracting Molecular Weights
for i = 1:numel(D)
    MM5 = str2num(A(D(i)+19:D(i)+24)); 
end
% Extracting gamma values
for i = 1:numel(E)
    GAMMA5 = str2num(A(E(i)+18:E(i)+23)); 
end



PT5_N = mdot_5*sqrt(R_u*TT5/MM5);
PT5_D =CD_5*A5*GAMMA5^0.5*(2/(GAMMA5+1))^(GAMMA5+1/(2*(GAMMA5-1)));

PT5= PT5_N/PT5_D;


%%%%%%%%%%%% THE FOLLOWING ASSUMPTION SHOULD BE CHECKED%%%%%%%%%%%%%%%%%%%%
GAMMA6 =GAMMA5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms x
M6 = vpasolve((1/x)*(((2+(GAMMA6-1)*x^2)/(GAMMA6+1))^((GAMMA6+1)/(2*(GAMMA6-1)))) == A6/A5, x);
% M6 = solve((1/x)*(((2+(GAMMA6-1)*x^2)/(GAMMA6+1))^((GAMMA6+1)/(2*(GAMMA6-1)))) == A6/A5, x);

P6 = PT5/(1+(GAMMA6-1)/2*M6^2)^((GAMMA6)/(GAMMA6-1));

FGROSS = P6*A6*GAMMA6*M6^2+A6*(P6-P0);
FDRAG = P0*A0*GAMMA0*M0^2;
FNET = FGROSS-FDRAG;
FNET = double(FNET);

end