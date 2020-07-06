#!/usr/bin/env python3

# Author: Francesco Regazzoni - MOX, Politecnico di Milano
# Email:  francesco.regazzoni@polimi.it
# Date:   2020

import numpy as np
import time

class model_ANN:
    """
    Class implementing the ANN-based model for sarcomere activation
    proposed in [1]. The ANN model is trained from a database of simulations
    obtained with the RDQ18 model [2]. The ANN model has been trained by
    means of the Machine Learning method proposed in [3], using the library
    model-learning (https://github.com/FrancescoRegazzoni/model-learning).

    References
    ----------

    [1] F. Regazzoni, L. Dedè, A. Quarteroni "Machine learning of multiscale
        active force generation models for the efficient simulation of
        cardiac electromechanics", Computer Methods in Applied Mechanics and
        Engineering (2020)

    [2] F. Regazzoni, L. Dedè, A. Quarteroni "Active contraction of cardiac
        cells: a reduced model for sarcomere dynamics with cooperative
        interactions", Biomechanics and Modeling in Mechanobiology (2018)
        https://doi.org/10.1007/s10237-018-1049-0

    [3] F. Regazzoni, L. Dedè, A. Quarteroni "Machine learning for fast and
        reliable solution of time-dependent differential equations", Journal
        of Computational Physics (2019)
        https://doi.org/10.1016/j.jcp.2019.07.050
    """

    def __init__(self, ANN_path = '../ANN'):
        """
        Constructor.

        Parameters
        ----------
        ANN_path : path of the folder containing the ANN weights, ANN biases
                   and the initial state of the model.
        """

        self.x0 = np.genfromtxt(ANN_path + '/initial_state.csv', delimiter = ',')
        self.W0 = np.genfromtxt(ANN_path + '/weights_0.csv', delimiter = ',')
        self.W1 = np.genfromtxt(ANN_path + '/weights_1.csv', delimiter = ',')
        self.W2 = np.genfromtxt(ANN_path + '/weights_2.csv', delimiter = ',')
        self.T0 = np.genfromtxt(ANN_path + '/biases_0.csv', delimiter = ',')
        self.T1 = np.genfromtxt(ANN_path + '/biases_1.csv', delimiter = ',')
        self.T2 = np.genfromtxt(ANN_path + '/biases_2.csv', delimiter = ',')

        self.dt = 1e-3 # [s]

    def f(self, x, Ca, SL):
        inp = np.concatenate((np.array(Ca)[None],np.array(SL)[None],x))
        return np.matmul( self.W2, np.tanh( np.matmul(self.W1, np.tanh( np.matmul(self.W0, inp) - self.T0 ) ) - self.T1 ) ) - self.T2

    def solve(self, inputs):
        """
        Perform a simulation with the model.

        Parameters
        ----------
        inputs: dictionary containing the input data
          - inputs['times']: time instants [s]
          - inputs['Ca']:    intracellular calcium ions concentration [micro M]
          - inputs['SL']:    sarcomere length [micro m]

        Returns
        -------
        output: dictionary containing the output data
          - output['times']: time instants [s]
          - output['Ca']:    intracellular calcium ions concentration [micro M]
          - output['SL']:    sarcomere length [micro m]
          - output['P']:     permissivity [-]
        """

        times = np.arange(np.min(inputs['times']), np.max(inputs['times']), self.dt)
        Ca    = np.interp(times, inputs['times'], inputs['Ca'])
        SL    = np.interp(times, inputs['times'], inputs['SL'])
        nT    = len(times)
        P     = np.zeros(nT)

        x = self.x0
        P[0] = x[0]

        time_init = time.time();
        print('ANN model.   Computing... ', end = '')
        for iT in range(1, nT):
            x = x + self.dt * self.f(x, Ca[iT], SL[iT])
            P[iT] = x[0]
        print('done. Time elapsed: %1.3f s' % (time.time() - time_init))

        output = dict()
        output['times'] = times
        output['Ca']    = Ca
        output['SL']    = SL
        output['P']     = P

        return output