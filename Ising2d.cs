namespace Quantum.Ising2d {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    // open Quantum.Migrating;


    /// # Summary
    /// Applies coupling term - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ for all qubits minimising circuit depth
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## Jh
    /// 2d array of horizontal coupling coefficients
    /// ## Jv
    /// 2d array of vertical coupling coefficients
    /// ## nnonly
    /// bool whether or not there is only nearest neighbour qubit interactions
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    operation EvolveCouplings(nSites: Int, dt: Double, Jh: Double[][], Jv: Double[][], qubits: Qubit[]): Unit {
        let rows = Length(Jh);
        let cols = Length(Jv);
        for (r in 0 .. (rows - 1)) {
            for (c in 0..2..(cols - 2)) {
                let ind = r  * c;
                ApplyWithCA(CNOT(qubits[ind], _), Rz(-2.0 * Jh[r][c] * dt, _), qubits[ind + 1]);
            }
        }
        for (r in 0 .. rows - 1) {
            for (c in 1..2..(cols - 2)) {
                let ind = r  * c;
                ApplyWithCA(CNOT(qubits[ind], _), Rz(-2.0 * Jh[r][c] * dt, _), qubits[ind + 1]);
            }
        }
        for (c in 0 .. (cols - 1)) {
            for (r in 0..2..(rows - 2)) {
                let ind = r  * c;
                ApplyWithCA(CNOT(qubits[ind], _), Rz(-2.0 * Jv[r][c] * dt, _), qubits[ind + cols]);
            }
        }
        for (c in 0 .. (cols - 1)) {
            for (r in 0..1..(rows - 2)) {
                let ind = r  * c;
                ApplyWithCA(CNOT(qubits[ind], _), Rz(-2.0 * Jv[r][c] * dt, _), qubits[ind + cols]);
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
    /// ## Jh
    /// 2d array of horizontal coupling coefficients
    /// ## Jv
    /// 2d array of vertical coupling coefficients
    /// ## qubits
    /// Qubits that the encoded Ising Hamiltonian acts on.
    operation EvolveSingleTimestep(nSites: Int, dt: Double, g: Double[], h: Double[], Jh: Double[][], Jv: Double[][], qubits : Qubit[]): Unit {
        for (idxSite in 0 .. nSites - 1) {
            Rx((-2.0 * g[idxSite]) * dt, qubits[idxSite]);
            Rz((-2.0 * h[idxSite]) * dt, qubits[idxSite]);
            }
        EvolveCouplings(nSites, dt, Jh, Jv, qubits);
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
    /// ## Jh
    /// 2d array of horizontal coupling coefficients
    /// ## Jv
    /// 2d array of vertical coupling coefficients
    operation EvolveSingleTimestepDummy(nSites: Int, dt: Double, g: Double[], h: Double[], Jh: Double[][], Jv: Double[][]): Unit {
        using (qubits = Qubit[nSites]) {
            EvolveSingleTimestep(nSites, dt, g, h, Jh, Jv, qubits);
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
    /// ## Jh
    /// 2d array of horizontal coupling coefficients
    /// ## Jv
    /// 2d array of vertical coupling coefficients
    /// ## xinit
    /// bool whether to apply an initial Hadamard so as to initialise in the x basis
    /// ## xmeas
    /// bool whether to measure in the x basis
    ///
    /// # Output
    /// ## finalState
    /// 1d array of final states of each qubit in z basis (0 or 1)
    operation Evolve(initialState: Int[], time: Double, dt: Double, g: Double[], h: Double[], Jh: Double[][], Jv: Double[][], xinit: Bool, xmeas:Bool): Bool[] {
        nSites = Length(initialState);
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
                EvolveSingleTimestep(nSites, dt, g, h, J, qubits);
            }
            if (xmeas) {
                for (q in qubits) {
                    H(q);
                }
            }
            let resultArray = ForEach(MResetZ, qubits);
            mutable boolArray = Bool[nSites];
            for (idxSite in 0 .. nSites - 1) {
                if (resultArray[idxSite] == One) {
                   set boolArray w/= idxSite <- true; 
                }
            }
            return boolArray;
        }
    }
}
