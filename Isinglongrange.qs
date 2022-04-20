namespace Quantum.Isinglongrange {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    // Implement evolution under Trotterised Ising model as defined by
    //     H ≔ - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ - Σᵢ hᵢ Zᵢ - Σᵢ gᵢ Xᵢ


    /// # Summary
    /// Applies coupling term - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ for all qubits one site at a time (ie sweep j for each i)
    //
    /// # Input
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    operation EvolveCouplingsNaive(dt: Double, J: Double[][], qubits: Qubit[]): Unit {
        let nSites = Length(qubits);
        for (i in 0 .. nSites - 2) {
            for (j in i + 1 .. nSites - 1) {
                ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
            }
        }
    } 

    /// # Summary
    /// Applies coupling term - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ for all qubits minimising circuit depth and assuming all to all coupling
    ///
    /// # Input
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## mapping
    /// register of the site indices respresented by the qubits (-1 signifies an ancilla)
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    ///
    /// # Output
    /// ## newMapping
    /// the updated mapping register
    operation EvolveCouplingsAllToAll(dt: Double, J: Double[][], mapping: Int[], qubits: Qubit[]): Int[] {
        let nSites = Length(qubits);
        mutable idxTarget = 0;
        mutable newMapping = new Int[Length(mapping)];
        set newMapping = mapping;
        for (step in 0..(nSites - 1)) {
            if (step % 2 == 0) {
                for (idxSource in 0..(nSites/2 - 1)) {
                    set idxTarget = nSites - idxSource - 1;
                    ApplyWithCA(CNOT(qubits[newMapping[idxSource]], _),
                                Rz(-2.0 * J[newMapping[idxSource]][newMapping[idxTarget]] * dt, _),
                                qubits[newMapping[idxTarget]]);
                }
                for (idx in 0..((nSites + 1)/2 - 1)) {
                    set newMapping = Swapped(nSites - idx, nSites - idx - 1, newMapping);
                }
            }
            else {
                for (idxSource in 0..((nSites - 1)/2 - 1)) {
                    set idxTarget = nSites - idxSource - 1;
                    ApplyWithCA(CNOT(qubits[newMapping[idxSource]], _),
                                Rz(-2.0 * J[newMapping[idxSource]][newMapping[idxTarget]] * dt, _),
                                qubits[newMapping[idxTarget]]);
                }
                for (ind in 0..(nSites/2 - 1)) {
                    set newMapping = Swapped(nSites/2 - ind, nSites/2 - ind - 1, newMapping);
                }
                set newMapping = Swapped(0, nSites, newMapping);
            }
        }
        return newMapping;
    }



    /// # Summary
    /// Applies evolution for a single timestep.
    ///
    /// # Input
    /// ## dt
    /// Trotter time step size
    /// ## g
    /// 1d array of transverse field coefficients, gᵢ
    /// ## h
    /// 1d array of longitudinal field coefficients hᵢ
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## parallel
    /// bool whether or not to execute SWAP operationsin parallel
    /// ## mapping
    /// register of the site indices respresented by the qubits (-1 signifies an ancilla)
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    ///
    /// # Output
    /// ## newMapping
    /// the updated mapping register
    operation EvolveSingleTimestep(dt: Double, g: Double[], h: Double[], J: Double[][], parallel: Bool,
                                   mapping: Int[], qubits: Qubit[]): Int[] {
        let nMap = Length(mapping);
        for (idx in 0 .. nMap - 1) {
            if (mapping[idx] >= 0) {
                Rx((-2.0 * g[mapping[idx]]) * dt, qubits[mapping[idx]]);
                Rz((-2.0 * h[mapping[idx]]) * dt, qubits[mapping[idx]]);
            }
        }
        mutable newMapping = new Int[nMap];
        set newMapping = mapping;
        if (parallel) {
            set newMapping = EvolveCouplingsAllToAll(dt, J, mapping, qubits);
        } else {
            EvolveCouplingsNaive(dt, J, qubits);
        }
        return newMapping;
    }

    /// # Summary
    /// Applies evolution for a single timestep starting with freshly initialised qubits
    /// (ie useful for gate count but not evolution)
    ///
    /// # Input
    /// ## dt
    /// Trotter time step size
    /// ## g
    /// 1d array of transverse field coefficients, gᵢ
    /// ## h
    /// 1d array of longitudinal field coefficients hᵢ
    /// ## J
    /// 2d array of coupling coefficients Jᵢⱼ
    /// ## parallel
    /// bool whether or not to execute SWAP operationsin parallel
    operation EvolveSingleTimestepDummy(dt: Double, g: Double[], h: Double[], J: Double[][], parallel: Bool): Unit {
        let nSites = Length(g); 
        mutable nQubits = nSites;
        mutable nMap = nSites;
        if (parallel) {
            set nMap += 1;
        }
        mutable mapping = new Int[nMap];
        for (idx in 0 .. nSites - 1) {
            if (idx < nSites) {
                set mapping w/= idx <- idx;
            } else {
                set mapping w/= idx <- -1;
            }
        }
        using (qubits = Qubit[nQubits]) {
            set mapping = EvolveSingleTimestep(dt, g, h, J, parallel, mapping, qubits);
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
    /// ## parallel
    /// bool whether or not to execute SWAP operationsin parallel
    /// ## xinit
    /// bool whether to apply an initial Hadamard so as to initialise in the x basis
    /// ## xmeas
    /// bool whether to measure in the x basis
    ///
    /// # Output
    /// ## finalState
    /// 1d array of final states of each qubit in z basis (false or true)
    operation Evolve(initialState: Int[], time: Double, dt: Double, g: Double[],
        h: Double[], J: Double[][], parallel: Bool, xinit: Bool, xmeas:Bool): Bool[] {
        let nSites = Length(initialState); 
        mutable nQubits = nSites;
        mutable nMap = nSites;
        if (parallel) {
            set nMap += 1;
        }
        mutable mapping = new Int[nMap];
        using (qubits = Qubit[nQubits]) {
            for (idx in 0 .. nMap - 1) {
                if (idx < nSites) {
                    set mapping w/= idx <- idx;
                    if (initialState[idx] == 1) {
                        X(qubits[idx]);
                    }
                    if (xinit) {
                        H(qubits[idx]);
                    }
                } else {
                    set mapping w/= idx <- -1;
                }
            }
            let nSteps = Floor(time / dt);
            for (timestep in 0 .. nSteps - 1) {
                set mapping = EvolveSingleTimestep(dt, g, h, J, parallel, mapping, qubits);
            }
            mutable boolArray = new Bool[nSites];
            for (idx in 0 .. nMap - 1) {
                if (mapping[idx] >= 0) {
                    if (xmeas) {
                        H(qubits[idx]);
                    }
                    if (MResetZ(qubits[idx]) == One) {
                        set boolArray w/= mapping[idx] <- true;
                    }
                }
            }
            return boolArray;
        }
    }
}
