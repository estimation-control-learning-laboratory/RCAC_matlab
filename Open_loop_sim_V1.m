clc
clear all
close all


PHI = 0.20:0.01:1.2;
% PHI = [.8 1.1];
for ii = 1:length(PHI)

    [FNET(ii), GAMMA6(ii)] = Compute_SFRJ_thrust(PHI(ii));
end

%%
FNET
subplot(2,1,1)
stairs(PHI, FNET,'b*-','linewidth',2)
xlabel('PHI') 
ylabel('F_{NET}')
axis tight; grid on

subplot(2,1,2)
stairs(PHI, GAMMA6,'b*-','linewidth',2)
xlabel('PHI')
ylabel('GAMMA6')
axis tight; grid on