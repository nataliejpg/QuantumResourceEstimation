namespace Quantum.Isingext {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;


    operation Juxtapose(index1: Int, index2: Int, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in index2..-1..index1+2) {
            SWAP(qubits[idxSite], qubits[idxSite - 1]);
        }
    }

    operation Separate(index1: Int, index2: Int, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in index1+1..index2-1) {
            SWAP(qubits[idxSite], qubits[idxSite + 1]);
        }
    }


    /// # Summary
    /// Applies coupling term - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ for all qubits minimising circuit depth and assuming nearest neighbour coupling
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
    operation EvolveCouplingsParallel(dt: Double, J: Double[][], mapping: Int[], qubits: Qubit[]): Int[] {
        let nSites = (Length(qubits) - 1) / 2;
        mutable idxTarget = 0;
        mutable newMapping = new Int[Length(mapping)];
        set newMapping = mapping;
        for (step in 0..(nSites - 1)) {
            if (step % 2 == 0) {
                for (idxSource in 0..(nSites/2 - 1)) {
                    set idxTarget = nSites - idxSource - 1;
                    ApplyWithCA(CNOT(qubits[idxSource], _),
                                Rz(-2.0 * J[newMapping[idxSource]][newMapping[idxTarget]] * dt, _),
                                qubits[idxTarget]);
                }
                for (idx in 1..((nSites + 1)/2 - 1)) {
                    SWAP(qubits[nSites - idx], qubits[2 * nSites - idx + 1]);
                    SWAP(qubits[2 * nSites - idx + 1], qubits[nSites - idx - 1]);
                    set newMapping = Swapped(nSites - idx, nSites - idx - 1, newMapping);
                }
                SWAP(qubits[nSites], qubits[nSites - 1]);
                set newMapping = Swapped(nSites, nSites - 1, newMapping);

            } else {
                for (idxSource in 0..((nSites - 1)/2 - 1)) {
                    set idxTarget = nSites - idxSource - 1;
                    ApplyWithCA(CNOT(qubits[idxSource], _),
                                Rz(-2.0 * J[newMapping[idxSource]][newMapping[idxTarget]] * dt, _),
                                qubits[idxTarget]);
                }
                for (ind in 1..((nSites)/2 - 1)) {
                    SWAP(qubits[nSites/2 - ind], qubits[nSites/2 - ind + nSites + 1]);
                }
                SWAP(qubits[0], qubits[nSites + 1]);
                for (ind in 1..((nSites)/2 - 1)) {
                    SWAP(qubits[nSites/2 - ind + nSites + 1], qubits[nSites/2 - ind - 1]);
                    set newMapping = Swapped(nSites/2 - ind, nSites/2 - ind - 1, newMapping);
                }
                SWAP(qubits[nSites + 1], qubits[nSites]);
                set newMapping = Swapped(0, nSites, newMapping);
                SWAP(qubits[nSites/2], qubits[nSites/2 - 1]);
                set newMapping = Swapped(nSites/2, nSites/2 - 1, newMapping);
            }
        }
        return newMapping;
    }



    /// # Summary
    /// Applies coupling term - Σᵢⱼ Jᵢⱼ Zᵢ Zⱼ for all qubits minimising circuit depth in rotation gates
    /// but not in SWAPS and assuming nearest  neighbour coupling
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
    operation EvolveCouplingsSequential(dt: Double, J: Double[][], mapping: Int[], qubits: Qubit[]): Int[] {
        let nSites = Length(qubits) - 1;
        mutable idxTarget = 0;
        mutable newMapping = new Int[Length(mapping)];
        set newMapping = mapping;
        for (step in 0..(nSites - 1)) {
            if (step % 2 == 0) {
                for (idxSource in 0..(nSites/2 - 1)) {
                    set idxTarget = nSites - idxSource - 1;
                    ApplyWithCA(CNOT(qubits[idxSource], _),
                                Rz(-2.0 * J[newMapping[idxSource]][newMapping[idxTarget]] * dt, _),
                                qubits[idxTarget]);
                }
                for (idx in 0..((nSites + 1)/2 - 1)) {
                    SWAP(qubits[nSites - idx], qubits[nSites - idx - 1]);
                    set newMapping = Swapped(nSites - idx, nSites - idx - 1, newMapping);
                }
            } else{
                for (idxSource in 0..((nSites - 1)/2 - 1)) {
                    set idxTarget = nSites - idxSource - 1;
                    ApplyWithCA(CNOT(qubits[idxSource], _),
                                Rz(-2.0 * J[newMapping[idxSource]][newMapping[idxTarget]] * dt, _),
                                qubits[idxTarget]);
                }
                for (ind in 0..((nSites)/2 - 1)) {
                    SWAP(qubits[nSites/2 - ind], qubits[nSites/2 - ind - 1]);
                    set newMapping = Swapped(nSites/2 - ind, nSites/2 - ind - 1, newMapping);
                }
                SWAP(qubits[0], qubits[nSites]);
                set newMapping = Swapped(0, nSites, newMapping);
            }
        }
        return newMapping;
    }
}