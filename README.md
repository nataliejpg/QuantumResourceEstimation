# Quantum Resource Estimation

Q# code implementing dynamics of the Fermi-Hubbard and long range Ising models on digital quantum hardware in order to facilitate gate count estimates. The evolution can be run from the jupyter notebooks and compared with the results from the basic exact diagonalisation code written in python which can be found in the exact_diagonalisation folder. The comparison is intended for proof of principle and is only possible for small systems (<10 qubits). Similarly the Trotterisation of the models for implementation on quantum hardware is only 1st order which results in substantial Trotter errors but crucuially still allows for the resource estimation functions of Q# to be useful in determining the exact gate counts for the systems per Trotter time step sweep which can then be used to calculate the full evolution estimates.

### Requires:
- python 3.6.7
- numpy 1.18.1
- scipy 1.4.1
- qsharp >= 0.11.2004.2825
- notebook 6.0.1
optional for jupyter notebook example:
- jupyter
- matplotlib
- progressbar2

### Improvements
- Implement 2nd and 4th order Trotter decompositions so that the full evolution estimates can be directly calculated (and that the dynamics have lower errors so can be run with larger timesteps more quickly)
- Implement Q# depth counter
- Implement resource estimation for hardware with only nearest neighbour connectivity
- Extend to include other models (incl nearest neighbour ising with variable coupling coefficients)
