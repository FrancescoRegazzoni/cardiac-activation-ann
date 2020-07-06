function output = model_ANN(input, ANN_path)
%__________________________________________________________________________
%
% This code implements the ANN-based model for sarcomere activation
% proposed in [1]. The ANN model is trained from a database of simulations
% obtained with the RDQ18 model [2]. The ANN model has been trained by
% means of the Machine Learning method proposed in [3], using the library
% model-learning (https://github.com/FrancescoRegazzoni/model-learning).
%
% Inputs:
%  - input: struct containing the input data
%     - input.times: time instants [s]
%     - input.Ca:    intracellular calcium ions concentration [micro M]
%     - input.SL:    sarcomere length [micro m]
%  - ANN_path: path of the folder containing the ANN weights, ANN biases
%              and the initial state of the model.
%
% Outputs:
%  - output: struct containing the output data
%     - output.times: time instants [s]
%     - output.Ca:    intracellular calcium ions concentration [micro M]
%     - output.SL:    sarcomere length [micro m]
%     - output.P:     permissivity [-]
%
%__________________________________________________________________________
%
% Author: Francesco Regazzoni - MOX, Politecnico di Milano
% Email:  francesco.regazzoni@polimi.it
% Date:   2020
%__________________________________________________________________________
%
% References:
%
% [1] F. Regazzoni, L. Dedè, A. Quarteroni "Machine learning of multiscale
%     active force generation models for the efficient simulation of
%     cardiac electromechanics", Computer Methods in Applied Mechanics and
%     Engineering (2020)
%
% [2] F. Regazzoni, L. Dedè, A. Quarteroni "Active contraction of cardiac
%     cells: a reduced model for sarcomere dynamics with cooperative
%     interactions", Biomechanics and Modeling in Mechanobiology (2018)
%     https://doi.org/10.1007/s10237-018-1049-0
%
% [3] F. Regazzoni, L. Dedè, A. Quarteroni "Machine learning for fast and
%     reliable solution of time-dependent differential equations", Journal
%     of Computational Physics (2019)
%     https://doi.org/10.1016/j.jcp.2019.07.050
%__________________________________________________________________________


    %% ANN definition
    if nargin == 1
        ANN_path = '../ANN';
    end

    x0 = readmatrix([ANN_path '/initial_state.csv'])';
    W0 = readmatrix([ANN_path '/weights_0.csv']);
    W1 = readmatrix([ANN_path '/weights_1.csv']);
    W2 = readmatrix([ANN_path '/weights_2.csv']);
    T0 = readmatrix([ANN_path '/biases_0.csv'])';
    T1 = readmatrix([ANN_path '/biases_1.csv'])';
    T2 = readmatrix([ANN_path '/biases_2.csv'])';

    f = @(x,Ca,SL) W2*tanh(W1*tanh(W0*[Ca;SL;x]-T0)-T1)-T2;

    dt = 1e-3; % [s]

    %% Time loop
    times = min(input.times):dt:max(input.times);
    Ca    = interp1(input.times, input.Ca, times);
    SL    = interp1(input.times, input.SL, times);
    nT    = length(times);
    P     = zeros(1,nT);

    x = x0;
    P(1) = x(1);

    time_init = tic();
    fprintf('ANN model.   Computing... ')
    for iT = 2:nT
        x = x + dt*f(x,Ca(iT),SL(iT));
        P(iT) = x(1);
    end
    fprintf('done. Time elapsed: %1.3f s\n', toc(time_init))

    output.times = times;
    output.Ca    = Ca;
    output.SL    = SL;
    output.P     = P;

end