#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from model_RDQ18 import model_RDQ18
from model_ANN   import model_ANN

#%% Time interval
Tmax = .8   # [s]

#%% Calcium transient
c0 = .1     # [micro M]
cmax = 1.1  # [micro M]
tau1 = .02  # [s]
tau2 = .05  # [s]
t0 = 0.01   # [s]
beta = (tau1/tau2)**(-1/(tau1/tau2 - 1)) - (tau1/tau2)**(-1/(1 - tau2/tau1))
Ca_base = lambda t: c0 + (t>=t0) * ((cmax - c0) / beta * (np.exp(-(t-t0)/tau1) - np.exp(-(t-t0)/tau2)))

#%% SL transient
SL0 = 2.2;     # [micro m]
SL1 = SL0*.95; # [micro m]
SLt0 = .15;    # [s]
SLt1 = .55;    # [s]
SLtau0 = .05;  # [s]
SLtau1 = .02;  # [s]

SL_base = lambda t: SL0 + (SL1-SL0) * (np.clip(1-np.exp((SLt0-t)/SLtau0),0,None) - np.clip(1-np.exp((SLt1-t)/SLtau1),0,None))

#%% Input definition
inputs = dict()
inputs['times'] = np.arange(0,Tmax,1e-4)
inputs['Ca'] = Ca_base(inputs['times']);
inputs['SL'] = SL_base(inputs['times']);

#%% Simulation
output_RDQ18 = model_RDQ18().solve(inputs)
output_ANN   = model_ANN().solve(inputs)

#%% Postprocessing

fig, axes = plt.subplots(3, 1, figsize=(3, 5))

axes[0].plot(output_RDQ18['times'], output_RDQ18['Ca'], 'k')
axes[0].set_xlabel('time [s]')
axes[0].set_ylabel('$[\mathrm{Ca}^{2+}]_i \, [\mu M]$')
axes[0].set_ylim([0, 1.5])

axes[1].plot(output_RDQ18['times'], output_RDQ18['SL'], 'k')
axes[1].set_xlabel('time [s]')
axes[1].set_ylabel('$SL \, [\mu m]$')
axes[1].set_ylim([2, 2.3])

axes[2].plot(output_RDQ18['times'], output_RDQ18['P'], label = 'RDQ18')
axes[2].plot(output_ANN['times'], output_ANN['P'], '--', label = 'ANN model')
axes[2].set_xlabel('time [s]')
axes[2].set_ylabel('permissivity [-]')
axes[2].set_ylim([0, 0.4])
axes[2].legend()

fig.tight_layout()
fig.savefig('output.pdf')