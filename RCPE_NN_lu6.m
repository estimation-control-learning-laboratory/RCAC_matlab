clc
clear all
close all
format longe
format
warning('off','all')
tic

% randn('state',9)
% rand('state',2)
addpath '../'
addpath './RCAC_functions/'
LoadFigurePrintingProperties


steps           = 10000;
FLAG.steps      = steps;

FLAG.Nc         = 0;
FLAG.window     = steps; 300;

FLAG.optimType  = 'Batch';
FLAG.optimType  = 'RLS';
FLAG.R0         = 1e-8;
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


X = [0 0 1 1
    0 1 0 1];
Y = [0 1 1 0];

lx  = size(X,2);

% The XOR function has 4 possible values
% Take one of these by random and use it for error calculations
lz  = size(1,1);
ly  = size(1,1);


% theta_nn0 = randn(1,6);

%manually set 6 of them
theta_nn0(1) = 4.510866097865148e+01;

theta_nn0(2)= -4.481706269421458e+01;
theta_nn0(3)= -2.311066575047726e+01;
theta_nn0(4)= -4.584524260712102e+01;
theta_nn0(5)=   4.522449892595682e+01;
theta_nn0(6)= -2.293930732337328e+01;

Theta = [theta_nn0(1:3)' theta_nn0(4:6)'];
for kk = 1:4
    Yhat_f(1,kk) = [1 1]*neural_layer(X(:,kk),Theta);
end



%randomly initialize for rcac
% theta_nn0(1) = .800;
%%

up = [1 2 3 4 5 6];
lu          = numel(up);
ltheta      = CalculateRegSize( FLAG.Nc, lu, lz, ly, FLAG);


FLAG.Ru     = 0e-5;
FLAG.R0     = 1e6;
FLAG.lambda = 1;

FILT.TYPE   = 'TF';

FILT.Nu     = -1;
FILT.Nu     = -[1 0 0 0 1 0 0 0 1];
FILT.Nu     = -[1 0 0 0 0 1 0 1 0];
FILT.Nu     = 1*vec(eye(6)*randn(6,1))';
FILT.Nf     = size(FILT.Nu,2)/lu;


u           = zeros(lu, steps);
y0          = zeros(ly,steps);


y0_hat      = y0;
z           = zeros(lz,steps);
theta       = zeros(ltheta,steps);
theta_nn = theta_nn0;
theta_nn(up) = 0;
%% Simulation
for ii = 1:steps

    if ii == 1


        Theta = [theta_nn(1:3)' theta_nn(4:6)'];
        z(:,ii) = compute_NN_RCPE_error(Theta);
        
        [u(:,ii), theta(:,ii)] = ...
            RCAC_V6(ii, zeros(lu,1), 0*z(:,ii), 0*z(:,ii), 0, FILT, FLAG);
    else
        theta_nn(up) =  u(:,ii-1)/100;
        Theta = [theta_nn(1:3)' theta_nn(4:6)'];
        z(:,ii) = compute_NN_RCPE_error(Theta);

        [u(:,ii), theta(:,ii)] = ...
            RCAC_V6(ii, u(:,ii-1), z(:,ii-1), 0*z(:,ii-1), 0, FILT, FLAG);

    end
%     if abs(sum(u(:,ii))>1e2) || abs(z(:,ii))>1e3
%         break
%     end
end



for kk = 1:4
    Yhat_f(1,kk) = [1 1]*neural_layer(X(:,kk),Theta);
end
Theta
[Yhat_f' Y']
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
xlabel('\rm{iterations}');
% set(gca, 'ytick', [1e-12 1e-9 1e-6 1e-3 1])
%
subplot(irow,icol,2)
%
% plot(time, Aup(:,1:ii)','k-', 'linewidth',2)
grid on; axis tight; hold on
plot(time, u(:,1:ii)', '--')
% for pp = 1:lu
%     H(pp) = plot(time, abs(u(pp,1:ii))', 'color', co(pp,:));
% end
%
% aa = legend([H],{'$\hat \mu_1$', '$\hat \mu_2$', '$\hat \mu_3$'});
aa = legend("$u(1)$");
set(aa, 'interpreter', 'latex', 'box', 'off', 'location', 'best', 'fontsize',8, 'orientation', 'horizontal')
%ylabel('$\hat \mu(k)$')
hold off


subplot(irow,icol,3)
plot(time,theta')
aa = legend("$\theta$");
set(aa, 'interpreter', 'latex', 'box', 'off', 'location', 'best', 'fontsize',8, 'orientation', 'horizontal')

grid on; axis tight; hold on;

% plot(time,theta(2,1:ii)')
% subplot(irow,icol,3)
% for jj = 1:ii
%     xerr(jj) = norm(x(:,jj)-x_hat(:,jj));
% end
% plot(time, xerr, 'b')
% grid on; axis tight; hold on
% ylabel('$\| x(k)-\hat{x}(k) \|_2$')
% set(gca, 'yscale', 'log')
% set(gca, 'ytick', [1e-12 1e-9 1e-6 1e-3 1])
%
% subplot(irow,icol,4)
% plot(time, theta(:,1:ii)')
% grid on; axis tight;
% ylabel('$\theta(k)$')
%
% ABCD = 'abcd';
% for pp = 1:2
%     subplot(irow,icol,pp)
%     set(gca,'XTickLabel', [])
% end
% for pp = 1:4
%     subplot(irow,icol,pp)
%     %set(gca, 'xtick', [10000 20000])
%     %xlabel(['(' ABCD(pp) ') ' 'Step, $k$'])
%     xlabel(['(' ABCD(pp) ')'])
%
% end
%
%
%
%
%
%
%
%
%
%
%
%


