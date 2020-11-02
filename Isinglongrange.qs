namespace Quantum.Isinglongrange {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Quantum.Migrating;

    // Implement evolution under Trotterised Ising model as defined by
    //     H ≔ - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ - Σᵢ hᵢ Zᵢ - Σᵢ gᵢ Xᵢ


    /// # Summary
    /// Applies coupling term - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ for all qubits one site at a time (ie sweep j for each i)
    /// # Input
    /// ## nSites
    /// Number of sites in the Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## nnonly
    /// bool whether or not there is only nearest neighbour qubit interactions
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    operation EvolveCouplings(nSites: Int, dt: Double, J: Double[][], nnonly: Bool, qubits: Qubit[]): Unit {
        for (i in 0 .. nSites - 2) {
            for (j in i + 1 .. nSites - 1) {
                if  (nnonly) {
                    Juxtapose(i, j, qubits);
                    ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[i + 1]);
                    Separate(i, j, qubits);
                } else {
                    ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
                }
            }
        }
    } 

    /// # Summary
    /// Applies coupling term - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ for all qubits minimising circuit depth
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## nnonly
    /// bool whether or not there is only nearest neighbour qubit interactions
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    operation EvolveCouplingsNested(nSites: Int, dt: Double, J: Double[][], nnonly: Bool, qubits: Qubit[]): Unit {
        let num = nSites/2 + nSites%2;
        for (step in 0 .. num) {
            for (ind in 0 .. nSites/2 - 1) {
                let i = ind - step;
                let j = nSites - (ind + step + 1);
                if  (nnonly) {
                    Juxtapose(i, j, qubits);
                    ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[i + 1]);
                    Separate(i, j, qubits);
                } else {
                    ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
                }
            }
            if ((nSites%2 == 0) or (step < nSites/2)){
                for (ind in 0 .. nSites/2 - 1) {
                    let i = ind - step;
                    let j = nSites - (step + ind + 2);
                    if (ind == (num - 1)) {
                        if  (nnonly) {
                            Juxtapose(i, j + num, qubits);
                            ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j + num] * dt, _), qubits[i + 1]);
                            Separate(i, j + num, qubits);
                        } else {
                            ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j + num] * dt, _), qubits[j + num]);
                        }
                    } else {
                        if  (nnonly) {
                            Juxtapose(i, j, qubits);
                            ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[i + 1]);
                            Separate(i, j, qubits);
                        } else {
                            ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
                        }
                    }
                    
                }
            }
        }
    }

    /// # Summary
    /// Applies evolution for a single timestep.
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## g
    /// 1d array of transverse field coefficients, gᵢ
    /// ## h
    /// 1d array of longitudinal field coefficients hᵢ
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## nested
    /// bool value of whether or not to reorder coupling terms to minimise circuit depth
    /// ## nnonly
    /// bool whether or not there is only nearest neighbour qubit interactions
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    operation EvolveSingleTimestep(nSites: Int, dt: Double, g: Double[], h: Double[], J: Double[][], nested: Bool, nnonly: Bool, qubits : Qubit[]): Unit {
        for (idxSite in 0 .. nSites - 1) {
            Rx((-2.0 * g[idxSite]) * dt, qubits[idxSite]);
            Rz((-2.0 * h[idxSite]) * dt, qubits[idxSite]);
            }
        if (nested) {
            EvolveCouplingsNested(nSites, dt, J, nnonly, qubits);
        } else {
            EvolveCouplings(nSites, dt, J, nnonly, qubits);
        }
    }

    /// # Summary
    /// Applies evolution for a single timestep starting with freshly initialised qubits
    /// (ie useful for gate count but not evolution)
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## g
    /// 1d array of transverse field coefficients, gᵢ
    /// ## h
    /// 1d array of longitudinal field coefficients hᵢ
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## nested
    /// bool value of whether or not to reorder coupling terms to minimise circuit depth
    /// ## nnonly
    /// bool whether or not there is only nearest neighbour qubit interactions
    operation EvolveSingleTimestepDummy(nSites: Int, dt: Double, g: Double[], h: Double[], J: Double[][], nested: Bool, nnonly: Bool): Unit {
        using (qubits = Qubit[nSites]) {
            EvolveSingleTimestep(nSites, dt, g, h, J, nested, nnonly, qubits);
        }
    }
    /// # Summary
    /// Applies full evolution in steps of dt.
    ///
    /// # Input
    /// ## initialState
    /// 1d array of initial states of each qubit in z basis (0 or 1)
    /// ## time
    /// Total time for evolution
    /// ## dt
    /// Trotter time step size
    /// ## g
    /// 1d array of transverse field coefficients, gᵢ
    /// ## h
    /// 1d array of longitudinal field coefficients hᵢ
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## nested
    /// bool value of whether or not to reorder coupling terms to minimise circuit depth
    /// ## nnonly
    /// bool whether or not there is only nearest neighbour qubit interactions
    /// ## xinit
    /// bool whether to apply an initial Hadamard so as to initialise in the x basis
    /// ## xmeas
    /// bool whether to measure in the x basis
    ///
    /// # Output
    /// ## finalState
    /// 1d array of final states of each qubit in z basis (0 or 1)
    operation Evolve(initialState: Int[], time: Double, dt: Double, g: Double[], h: Double[], J: Double[][], nested: Bool, nnonly: Bool, xinit: Bool, xmeas:Bool): Result[] {
        let nSites = Length(initialState);
        using (qubits = Qubit[nSites]) {
            for (idxSite in 0 .. nSites - 1) {
                if (initialState[idxSite] == 1) {
                    X(qubits[idxSite]);
                }
                if (xinit) {
                    H(qubits[idxSite]);
                }
            }
            let nSteps = Floor(time / dt);
            for (idxIter in 0 .. nSteps - 1) {
                EvolveSingleTimestep(nSites, dt, g, h, J, nested, nnonly, qubits);
            }
            if (xmeas) {
                for (q in qubits) {
                    H(q);
                }
            }
            return ForEach(MResetZ, qubits);
        }
    }
}
