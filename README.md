# Neural Networks based model of cardiac tissue activation #

This repository contains the trained data of the ANN (Artificial Neural Network) model of cardiac tissue activation presented in [1].
The ANN model is trained from a database of simulations obtained with the RDQ18 model [2].
The ANN model has been trained by means of the Machine Learning method proposed in [3], using the library [model-learning](https://github.com/FrancescoRegazzoni/model-learning).


### Content
- `ANN`: weights and biases of the ANN model (see Eq. (17) of [1]).
- `params`: parameters of the RDQ18 model [2].
- `matlab`: Matlab code with an example of simulation.
- `python`: Python code with an example of simulation.

### References

- [1] F. Regazzoni, L. Dedè, A. Quarteroni ["Machine learning of multiscale active force generation models for the efficient simulation of cardiac electromechanics"](https://doi.org/10.1016/j.cma.2020.113268), *Computer Methods in Applied Mechanics and Engineering* (2020).
- [2] F. Regazzoni, L. Dedè, A. Quarteroni ["Active contraction of cardiac cells: a reduced model for sarcomere dynamics with cooperative interactions"](https://doi.org/10.1007/s10237-018-1049-0), *Biomechanics and Modeling in Mechanobiology* (2018).
- [3] F. Regazzoni, L. Dedè, A. Quarteroni ["Machine learning for fast and reliable solution of time-dependent differential equations"](https://doi.org/10.1016/j.jcp.2019.07.050), *Journal of Computational Physics* (2019).

### Author

Francesco Regazzoni, MOX - Politecnico di Milano (<francesco.regazzoni@polimi.it>)
