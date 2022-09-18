clc
clear all
close all
format longe
% format
warning('off','all')
tic

% randn('state',9)
% rand('state',2)
addpath '../'
addpath './RCAC_functions/'
LoadFigurePrintingProperties


steps           = 50000;
FLAG.steps      = steps;

FLAG.Nc         = 0;
FLAG.window     = steps; 300;

FLAG.optimType  = 'Batch';
FLAG.optimType  = 'RLS';
FLAG.R0         = 1e+1;
FLAG.Rz         = 1;
FLAG.Ru         = 0e+1;

% Controller structure
FLAG.RegZ       = 1;            % 0->y goes in the controller.
FLAG.FF         = 0;
FLAG.Integrator = 1;
FLAG.FIR        = 0;
FLAG.ContType   = 'dense';
% FLAG.ContType   = 'sparse';
% FLAG.ContType   = 'supersparse';
% FLAG.ContType   = 'ManySISO';
% FLAG.ContType   = 'StepInjectionFIR';


%%
A   = [0.3 0.2; 0.1 0.6];
D   = [0; 1];
E   = [1 0];[1 2];
lx  = size(A,1);
lz  = size(E,1);
ly  = size(E,1);



%%

up          = [1 2];    % Uncertain Parameter
lu          = numel(up);
ltheta      = CalculateRegSize( FLAG.Nc, lu, lz, ly, FLAG);
A_hat       = A;
A_hat(up)   = 0;
D_hat       = D;
E_hat       = E;

ustar       = vec(A(up));
uhat        = 0*ustar;

FLAG.Ru     = 0e-5;
FLAG.R0     = 1e+7;
FLAG.lambda = 0.9995;

FILT.TYPE   = 'TF';
FILT.Nf     = 2;
FILT.Nu     = [1 0 0 1];


x           = zeros(lx, steps);
u           = zeros(lu, steps);
y0          = zeros(ly,steps);
x_hat       = x;
y0_hat      = y0;
z           = zeros(lz,steps);
theta       = zeros(ltheta,steps);

%% Simulation
for ii = 1:steps
    Aup(:,ii)   = ustar;

    w(ii) = 10+InputRCPE(ii,pi/50,1,15);
    if ii == 1
        x(:,ii)     = 1*ones(lx,1);
        y0(:,ii)    = E*x(:,ii);

        x_hat(:,ii) = 0*ones(lx,1);
        y0_hat(:,ii)= E*x_hat(:,ii);

        z(:,ii)     = y0(:,ii) - y0_hat(:,ii);

        [u(:,ii), theta(:,ii)] = ...
            RCAC_V6(ii, zeros(lu,1), 0*z(:,ii), y0(:,ii), 0, FILT, FLAG);
    else
        x(:,ii)     = A*x(:,ii-1) + D*w(:,ii-1);
        y0(:,ii)    = E*x(:,ii);

        x_hat(:,ii) = A_hat*x_hat(:,ii-1) + D_hat*w(:,ii-1);
        y0_hat(:,ii)= E_hat*x_hat(:,ii);

        z(:,ii)     = y0(:,ii) - y0_hat(:,ii);

        [u(:,ii), theta(:,ii),OptData] = ...
            RCAC_V6(ii, u(:,ii-1), z(:,ii-1), y0(:,ii-1), 0, FILT, FLAG);

        u(:,ii)     = Saturate(u(:,ii), -1,+1);
        para_est    = abs(u(:,ii));
        A_hat(up)   = para_est;

    end
    if abs(sum(u(:,ii))>1e2) || abs(z(:,ii))>1e3
        break
    end
end

disp([u(:,ii)' log10(abs(z(:,ii)')) toc]')




%%
close all
umin = 0;
umax = .7;
fname = ['JAIS_RCMR_lin_lu3_Ex1'];


set(0, 'DefaultFigureColor', 'default');
set(0, 'DefaultTextFontSize', 16);
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 12);
figure(43)
time = 1:ii;
irow    = 2;
icol    = 2;
co = [0 0 1;
    1 0 0;
    0 1 0];

% u = uTraj(:,:,25);
subplot(irow,icol,1)
plot(time, (abs(z(:,1:ii)))', 'b')
set(gca, 'yscale', 'log')
grid on; axis tight; hold on
ylabel('$| z(k) |$')
set(gca, 'ytick', [1e-12 1e-9 1e-6 1e-3 1])

subplot(irow,icol,2)

plot(time, Aup(:,1:ii)','k-', 'linewidth',2)
grid on; axis tight; hold on
plot(time, u(:,1:ii)', '--')
for pp = 1:lu
    H(pp) = plot(time, abs(u(pp,1:ii))', 'color', co(pp,:));
end

aa = legend([H],{'$\hat \mu_1$', '$\hat \mu_2$', '$\hat \mu_3$'});
set(aa, 'interpreter', 'latex', 'box', 'off', 'location', 'best', 'fontsize',8, 'orientation', 'horizontal')
%ylabel('$\hat \mu(k)$')
hold off



subplot(irow,icol,3)
for jj = 1:ii
    xerr(jj) = norm(x(:,jj)-x_hat(:,jj));
end
plot(time, xerr, 'b')
grid on; axis tight; hold on
ylabel('$\| x(k)-\hat{x}(k) \|_2$')
set(gca, 'yscale', 'log')
set(gca, 'ytick', [1e-12 1e-9 1e-6 1e-3 1])

subplot(irow,icol,4)
plot(time, theta(:,1:ii)')
grid on; axis tight;
ylabel('$\theta(k)$')

ABCD = 'abcd';
for pp = 1:2
    subplot(irow,icol,pp)
    set(gca,'XTickLabel', [])
end
for pp = 1:4
    subplot(irow,icol,pp)
    %set(gca, 'xtick', [10000 20000])
    %xlabel(['(' ABCD(pp) ') ' 'Step, $k$'])
    xlabel(['(' ABCD(pp) ')'])

end













%% functions
function [w] = InputRCPE(ii,omega_d, type, j_f)

% Input for RCPE

if type == 1
    w  = -2;
    for jj = 1:1:j_f
        w = w + sin(1*omega_d*jj*ii);
    end
    w   = w;
elseif type==2
    w  = -2;
    for jj = 1:1:j_f
        w = w + sin(1*omega_d*jj*ii+jj^2)/jj;   % shroeder signal
    end
    w   = w;
end

end

