namespace Quantum.Hubbard {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    // Implement evolution under Trotterised Hubbard model as defined by
    //     H ≔ - JΣ'ᵢⱼ Σσ (cᵢσ^\dagger cⱼσ + h.c) + U Σᵢnᵢup nᵢdown + μ Σᵢ Σσ nᵢσ
    // where Σ'ᵢⱼ means only nearest neighbour hopping

    /// # Summary
    /// Applies Clifford X/2 rotation gate so 0 -> Y and Y -> 1
    operation ClifY(qubit: Qubit): Unit is Adj + Ctl {
        Rx(PI() / 2., qubit);
    }

    /// # Summary
    /// Applies Clifford -X/2 rotation gate so 1 -> Y and Y -> 0
    operation ClifYInv(qubit: Qubit): Unit is Adj + Ctl {
        Rx(PI() / -2., qubit);
    }

    /// # Summary
    /// Applies Controlled-Z  gate between source q1 and target q2
    operation CZ(q1: Qubit, q2: Qubit): Unit is Adj + Ctl {
        H(q2);
        CNOT(q1, q2);
        H(q2);
    }

    /// # Summary
    /// Applies chemical potential terms (~μ Σᵢ Σσ nᵢσ)
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hubbard Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## mu
    /// Coefficient of the chemical potential, μ
    /// ## U
    /// Coefficient of the Coulomb repulsion, U
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on.
    operation ApplyChemicalPotentialTerms(nSites: Int, dt: Double, mu: Double, U: Double, qubits: Qubit[]): Unit {
        for (idxSite in 0 .. nSites - 1) {
            Rz((mu - U / 2.) * dt, qubits[idxSite]);
            Rz((mu - U / 2.) * dt, qubits[idxSite + nSites]);
            }
    } 

    /// # Summary
    /// Applies Coulomb repulsion terms (~U Σᵢnᵢup nᵢdown)
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hubbard Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## U
    /// Coefficient of the Coulomb repulsion, U
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on.
    operation ApplyCoulumbRepulsionTerms(nSites: Int, dt: Double, U: Double, qubits : Qubit[]): Unit {
        for (idxSite in 0 .. nSites - 1) {
            let q1 = qubits[idxSite];
            let q2 = qubits[idxSite + nSites];
            CNOT(q1, q2);
            Rz((0.5 * U) * dt, q2);
            CNOT(q1, q2);
            }
    }

    /// # Summary
    /// Applies controlled rotation part of a hopping term between two sites 
    /// of one species(ie w/o JW string), this is the full hopping term for sites where the 
    /// qubits representing the two sites are also adjacent.
    ///
    /// # Input
    /// ## q1
    /// Source qubit
    /// ## q1
    /// Target qubit
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// Coefficient of the hopping, J
    operation ApplyRotationSequence(q1: Qubit, q2: Qubit, dt: Double, J: Double): Unit {
        H(q1);
        H(q2);
        CNOT(q1, q2);
        Rz((-1. * J) * dt, q2);
        CNOT(q1, q2);
        H(q1);
        H(q2);
        ClifY(q1);
        ClifY(q2);
        CNOT(q1, q2);
        Rz((-1. * J) * dt, q2);
        CNOT(q1, q2);
        ClifYInv(q1);
        ClifYInv(q2);
    }

    /// # Summary
    /// Applies hopping term to a chain of qubits with nearest neighbour hopping
    /// - JΣ'ᵢⱼ Σσ (cᵢσ^\dagger cⱼσ + h.c)
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hubbard Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on.
    operation ApplyHoppingTermsChain(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit {
        for (idxSite in 0..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
    }

    /// # Summary
    /// Applies the hopping term between two nearest neighbour sites on the 'vertical ladder'
    /// which are represented by non adjacent qubits so a short JW string must be computed
    /// - JΣ'ᵢⱼ Σσ (cᵢσ^\dagger cⱼσ + h.c)
    ///
    /// # Input
    /// ## q1Idx
    /// Index of the source qubit
    /// ## q2Idx
    /// Index of the target qubit
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on.
    operation ApplyVerticalHoppingTerm(q1Idx: Int, q2Idx: Int, dt: Double, J: Double, qubits: Qubit[]): Unit {
        CNOT(qubits[q1Idx + 1], qubits[q1Idx + 2]);
        CZ(qubits[q2Idx - 1], qubits[q2Idx]);
        ApplyRotationSequence(qubits[q1Idx], qubits[q2Idx], dt, J);
        CZ(qubits[q2Idx - 1], qubits[q2Idx]);
        CNOT(qubits[q1Idx + 1], qubits[q1Idx + 2]);
    }

    /// # Summary
    /// Applies the hopping term between two rows of qubits representing two rows of sites of one 
    /// species. The hopping is nearest neighbour and the rotations are executed in parallel.
    /// The cumulative parity is stored on an ancilla qubit.
    /// - JΣ'ᵢⱼ Σσ (cᵢσ^\dagger cⱼσ + h.c)
    ///
    /// # Input
    /// ## q1Idxs
    /// Indices of the source qubits
    /// ## q2Idxs
    /// Indices of the target qubit
    /// ## ancillaIdx
    /// Index of the ancilla qubit
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on plus ancilla(s)
    operation ApplyRowHoppingTerms(q1Idxs: Int[], q2Idxs: Int[], ancillaIdx: Int, dt: Double, J: Double, qubits: Qubit[]): Unit {
        let rungs = Length(q1Idxs);
        for (rIdx in 0..rungs - 1) {
            CNOT(qubits[q1Idxs[rIdx] + 1], qubits[ancillaIdx]);
            CNOT(qubits[q2Idxs[rIdx] - 1], qubits[ancillaIdx]);
            CZ(qubits[ancillaIdx], qubits[q2Idxs[rIdx]]);
        }
        for (rIdx in 0..(rungs - 1)) {
            let q1 = qubits[q1Idxs[rIdx]];
            let q2 = qubits[q2Idxs[rIdx]];
            ApplyRotationSequence(q1, q2, dt, J);
        }
        for (rIdx in (rungs - 1)..-1..0) {
            CZ(qubits[ancillaIdx], qubits[q2Idxs[rIdx]]);
            CNOT(qubits[q2Idxs[rIdx] - 1], qubits[ancillaIdx]);
            CNOT(qubits[q1Idxs[rIdx] + 1], qubits[ancillaIdx]);
        }
    }

    /// # Summary
    /// Applies the hopping terms for a ladder architecture of the Hubbard model assuming that the
    /// qubits are ordered in a snake on the 'vertical' ladder.
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hubbard Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on.
    operation ApplyHoppingTermsLadderVertical(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit {
        ApplyHoppingTermsChain(nSites, dt, J, qubits);
        for (idxSite in 0..4..nSites - 4) {
            ApplyVerticalHoppingTerm(idxSite, idxSite + 3, dt, J, qubits);
            ApplyVerticalHoppingTerm(idxSite + nSites, idxSite + nSites + 3, dt, J, qubits);
        }
        for (idxSite in 2..4..nSites - 4) {
            ApplyVerticalHoppingTerm(idxSite, idxSite + 3, dt, J, qubits);
            ApplyVerticalHoppingTerm(idxSite + nSites, idxSite + nSites + 3, dt, J, qubits);
        }
    }

    /// # Summary
    /// Applies the hopping terms for a ladder architecture of the Hubbard model assuming that the
    /// qubits are ordered in a snake on the 'horizontal' ladder.
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hubbard Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on plus ancilla
    operation ApplyHoppingTermsLadderHorizontal(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit {
        ApplyHoppingTermsChain(nSites, dt, J, qubits);
        let length = nSites / 2;
        let upperRowIndsup = Reversed(SequenceI(0, length - 2));
        let upperRowIndsdown = Reversed(SequenceI(nSites, nSites + length - 2));
        let lowerRowIndsup = SequenceI(length + 1, nSites - 1);
        let lowerRowIndsdown = SequenceI(nSites + length + 1, 2 * nSites - 1);  
        ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, 2 * nSites, dt, J, qubits);
        ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, 2 * nSites + 1, dt, J, qubits);
    }

    /// # Summary
    /// Applies the hopping terms for a 2D lattice with qubits ordered in snake pattern.
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hubbard Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on plus ancillas
    operation ApplyHoppingTermsLattice(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit {
        ApplyHoppingTermsChain(nSites, dt, J, qubits);
        let length = Floor(Sqrt(IntAsDouble(nSites)));
        for (row in 0..2..(length - 2)){
            let upperRowIndsup = Reversed(SequenceI(row * length, row * length + length - 2));
            let upperRowIndsdown = Reversed(SequenceI(nSites + row * length, nSites + row * length + length - 2));
            let lowerRowIndsup = SequenceI(row * length + length + 1, row * length + 2 * length - 1);
            let lowerRowIndsdown = SequenceI(nSites + row * length + length + 1, nSites + row * length + 2 * length - 1);  
            let ancillaIdxup = 2 * nSites + row;
            let ancillaIdxdown = 2 * nSites + row + 1;
            ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, ancillaIdxup, dt, J, qubits);
            ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, ancillaIdxdown, dt, J, qubits);
        }
        for (row in 1..2..(length - 2)){
            let upperRowIndsup = Reversed(SequenceI(row * length, row * length + length - 2));
            let upperRowIndsdown = Reversed(SequenceI(nSites + row * length, nSites + row * length + length - 2));
            let lowerRowIndsup = SequenceI(row * length + length + 1, row * length + 2 * length - 1);
            let lowerRowIndsdown = SequenceI(nSites + row * length + length + 1, nSites + row * length + 2 * length - 1); 
            let ancillaIdxup = 2 * nSites + row - 1;
            let ancillaIdxdown = 2 * nSites + row;
            ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, ancillaIdxup, dt, J, qubits);
            ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, ancillaIdxdown, dt, J, qubits);
        }
    }

    /// # Summary
    /// Applies evolution for a single timestep.
    ///
    /// # Input
    /// ## nSites
    /// Number of sites in the Hubbard Hamiltonian.
    /// ## dt
    /// Trotter time step size
    /// ## mu
    /// Coefficient of the chemical potential, μ
    /// ## U
    /// Coefficient of the Coulomb repulsion, U
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## structure
    /// integer labeling structure, 1: chain, 2: vertical ladder, 3: horizontal ladder, 4: 2D lattice
    /// ## qubits
    /// Qubits that the encoded Hubbard Hamiltonian acts on plus any ancilla qubits
    operation EvolveSingleTimestep(nSites: Int, dt: Double, mu: Double, U: Double, J: Double, structure: Int, qubits: Qubit[]): Unit {
        ApplyChemicalPotentialTerms(nSites, dt, mu, U, qubits);
        ApplyCoulumbRepulsionTerms(nSites, dt, U, qubits);
        if (structure == 1) {
            ApplyHoppingTermsChain(nSites, dt, J, qubits);
        } elif (structure == 2) {
            ApplyHoppingTermsLadderVertical(nSites, dt, J, qubits);
        } elif (structure == 3) {
            ApplyHoppingTermsLadderHorizontal(nSites, dt, J, qubits);
        } else {
            ApplyHoppingTermsLattice(nSites, dt, J, qubits);
        }
    }

    /// # Summary
    /// Applies evolution for a single timestep starting with freshly initialised qubits
    /// (ie useful for gate count but not evolution)
    ///
    /// # Input
    /// ## initialState
    /// List of initial states of each site in the z basis (0: down down, 1: down up, 2: up down, 3: up up)
    /// ## dt
    /// Trotter time step size
    /// ## mu
    /// Coefficient of the chemical potential, μ
    /// ## U
    /// Coefficient of the Coulomb repulsion, U
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## structure
    /// integer labeling structure, 1: chain, 2: vertical ladder, 3: horizontal ladder, 4: 2D lattice
    operation EvolveSingleTimestepDummy(initialState: Int[], dt: Double, mu: Double, U: Double, J: Double, structure: Int): Unit {
        let nSites = Length(initialState);
        mutable qubitNum = 2 * nSites;
        if (structure == 3) {
            set qubitNum += 2;
        } elif (structure == 4) {
            let length = Floor(Sqrt(IntAsDouble(nSites)));
            set qubitNum += 2 * (length / 2);
        }
        using (qubits = Qubit[qubitNum]) {
            EvolveSingleTimestep(nSites, dt, mu, U, J, structure, qubits);
        }
    }

    /// # Summary
    /// Applies full evolution in steps of dt on structure of choice, allocating ancillas where necessary.
    ///
    /// # Input
    /// ## initialState
    /// List of initial states of each site in the z basis (0: down down, 1: down up, 2: up down, 3: up up)
    /// ## time
    /// Total time for evolution
    /// ## dt
    /// Trotter time step size
    /// ## mu
    /// Coefficient of the chemical potential, μ
    /// ## U
    /// Coefficient of the Coulomb repulsion, U
    /// ## J
    /// Coefficient of the hopping term, J
    /// ## structure
    /// integer labeling structure, 1: chain, 2: vertical ladder, 3: horizontal ladder, 4: 2D lattice
    ///
    /// # Output
    /// ## finalState
    /// list comparible to initialState of states of each site measured in z basis
    operation Evolve(initialState: Int[], time: Double, dt: Double, mu: Double, U: Double, J: Double, structure: Int): Int[] {
        let nSites = Length(initialState);
        mutable qubitNum = 2 * nSites;
        if (structure == 3) {
            set qubitNum += 2;
        } elif (structure == 4) {
            let length = Floor(Sqrt(IntAsDouble(nSites)));
            set qubitNum += 2 * (length / 2);
        }
        using (qubits = Qubit[qubitNum]) {
            for (idxSite in 0 .. nSites - 1) {
                if (initialState[idxSite] % 2 == 1) {
                    X(qubits[idxSite + nSites]);
                }
                if (initialState[idxSite] > 1) {
                    X(qubits[idxSite]);
                }
            }
            let nSteps = Floor(time / dt);
            for (idxIter in 0 .. nSteps - 1) {
                EvolveSingleTimestep(nSites, dt, mu, U, J, structure, qubits);
            }
            let result = ForEach(MResetZ, qubits);
            mutable finalState = new Int[nSites];
            for (idxSite in 0 .. nSites - 1) {
                if (result[idxSite] == One) {
                    set finalState w/= idxSite <- 2;
                }
                if (result[idxSite + nSites] == One) {
                    let newState = finalState[idxSite] + 1;
                    set finalState w/= idxSite <- newState;
                }
            }
            return finalState;
        }
    }

}