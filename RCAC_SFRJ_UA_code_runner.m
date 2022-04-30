clc
clear all
close all
format long
warning('off','all')
tic
addpath './RCAC_functions/'
LoadFigurePrintingProperties
randn('state',2)
steps           = 5000;

FLAG.Nc         = 2;
FLAG.window     = steps; 

% RCAC weights
FLAG.Rz         = 1;
FLAG.Ru         = 0e+1;
FLAG.R0         = 1e+4;
FLAG.lambda     = 1;


% Controller structure
FLAG.RegZ       = 1;            % 0->y goes in the controller.
FLAG.FF         = 0;
FLAG.Integrator = 1;
FLAG.FIR        = 0;
FLAG.ContType   = 'dense';
% FLAG.ContType   = 'PID';

%% System parameters
nn = abs(randn);
static_map = @(x)(x^nn);

%% Filter settings (FLAG)
lu = 1;
lz = 1;
ly = 1;

FILT.TYPE   = 'TF'; 
FILT.Nu     = -1;
FILT.Nf     = size(FILT.Nu,2)/lu;

%% Memory
ltheta      = CalculateRegSize( FLAG.Nc, lu, lz, ly, FLAG);
u           = zeros(lu, steps);
y           = zeros(ly,steps);
z           = zeros(lz,steps);
theta       = zeros(ltheta,steps);
r           = zeros(ly,steps)+1.2;


%% Simulation
PHI0 = 0.5;
r_jump = randn(3,1);
r_jump = [0.1 0.1];
for ii = 1:steps
    
    if ii>steps/3
        r(ii) = r(ii) + r_jump(1);
    end
    if ii>2*steps/3
        r(ii) = r(ii) + r_jump(2);
    end
%     r(:,ii) = 1.2;
    if ii == 1
        %% Initialization
        
        u(:,1) = zeros(lu,1);
        [FNET, GAMMA6] = Compute_SFRJ_thrust(PHI0+u(:,1),0);
        y(:,ii) = GAMMA6;
        z(:,ii) = y(:,ii) - r(:,ii);
        
        [u(:,ii), theta(:,ii)] = ...
            RCAC_V6(ii, zeros(lu,1), 0*z(:,ii), 0*z(:,ii), r(:,ii), FILT, FLAG);
        
    else
        y(:,ii) = static_map(u(:,ii-1));
        [FNET, GAMMA6] = Compute_SFRJ_thrust(PHI0+u(:,ii-1),0);
        y(:,ii) = GAMMA6;
        z(:,ii) = y(:,ii) - r(:,ii);
        
        [u(:,ii), theta(:,ii)] = ...
            RCAC_V6(ii, u(:,ii-1), z(:,ii-1), 0*z(:,ii-1), r(:,ii), FILT, FLAG);
        
    end
    %[ii GAMMA6 u(:,ii)]
    u(:,ii) = 1*Saturate(u(:,ii), -1, 2);
end

%%
close all
tic
Ts = 1;
time = (1:ii)*Ts;
%dTheta = A(up);
irow    = 2;
icol    = 2;
   
    
subplot(irow,icol,2)
    stairs(time, u(:,1:ii)', 'linewidth',2)
    grid on; axis tight; hold on
    ylabel('$u$')
    hold off
subplot(irow,icol,3)
    stairs(time, log10(abs(z(:,1:ii))'), 'linewidth',2)
    grid on; axis tight; hold on
    ylabel('$\log_{10}|z|$')
    hold off

%     xlabel('(b)')
subplot(irow,icol,1)
    stairs(time, r(:,1:ii)', 'k--','linewidth',2)    
    hold on; grid on; axis tight
    stairs(time, y(:,1:ii)', 'linewidth',2)
    hold off
    ylabel('$y$')
    aa = legend('Command', 'Output');
    set(aa, 'box','off','location','best')
    
    tt = ['$y=u^{' num2str(nn) '}$'];
    title(tt)
%     fix_legend(aa)
%     xlabel('(a)')
    
subplot(irow,icol,4)
    stairs(time, theta(:,1:ii)', 'linewidth',2)
    grid on; axis tight; 
    ylabel('$\theta$')
%     xlabel('(c)')
    
    


