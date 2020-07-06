clear

%% Time interval
Tmax = .8;   % [s]

%% Calcium transient
c0 = .1;     % [micro M]
cmax = 1.1;  % [micro M]
tau1 = .02;  % [s]
tau2 = .05;  % [s]
t0 = 0.01;   % [s]
beta = (tau1/tau2)^(-1/(tau1/tau2 - 1)) - (tau1/tau2)^(-1/(1 - tau2/tau1));
Ca_base = @(t) c0 + (t>=t0) .* ((cmax - c0) / beta * (exp(-(t-t0)/tau1) - exp(-(t-t0)/tau2)));

%% SL transient
SL0 = 2.2;     % [micro m]
SL1 = SL0*.95; % [micro m]
SLt0 = .15 ;   % [s]
SLt1 = .55;    % [s]
SLtau0 = .05;  % [s]
SLtau1 = .02;  % [s]
SL_base = @(t) SL0+ (SL1-SL0) * (max(0,1-exp((SLt0-t)/SLtau0)) - max(0,1-exp((SLt1-t)/SLtau1)));
    
%% Input definition
input.times = 0:1e-4:Tmax;
input.Ca = Ca_base(input.times);
input.SL = SL_base(input.times);

%% Simulation
output_RDQ18 = model_RDQ18(input);
output_ANN   = model_ANN(input);

%% Postprocessing
figure('units','pixel','outerposition', [100, 100, 250, 600]);

subplot(3,1,1);
plot(output_RDQ18.times, output_RDQ18.Ca, 'k', 'linewidth', 2)
ylim([0 1.5])
xlabel('time [s]')
ylabel('[Ca^{2+}]_i [\muM]')

subplot(3,1,2);
plot(output_RDQ18.times, output_RDQ18.SL, 'k', 'linewidth', 2)
ylim([2 2.3])
xlabel('time [s]')
ylabel('SL [\mum]')

subplot(3,1,3);
plot(output_RDQ18.times, output_RDQ18.P, 'linewidth', 2)
hold on
plot(output_ANN.times, output_ANN.P, '--', 'linewidth', 2)
legend('RDQ18', 'ANN model')
ylim([0 .4])
xlabel('time [s]')
ylabel('permissivity [-]')

set(gcf, 'PaperPosition', [0 0 3 5]); 
set(gcf, 'PaperSize', [3 5]);
saveas(gcf, 'output.pdf')