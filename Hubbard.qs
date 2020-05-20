namespace Quantum.Hubbard {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;


    operation ClifY(qubit: Qubit): Unit is Adj + Ctl {
        Rx(PI() / 2., qubit);
    }

    operation ClifYInv(qubit: Qubit): Unit is Adj + Ctl {
        Rx(PI() / -2., qubit);
    }

    operation CZ(q1: Qubit, q2: Qubit): Unit is Adj + Ctl {
        H(q2);
        CNOT(q1, q2);
        H(q2);
    }

    operation ApplyChemicalPotentialTerms(nSites: Int, dt: Double, mu: Double,
    U: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0 .. nSites - 1) {
            Rz(-(mu + U / 2.) * dt, qubits[idxSite]);
            Rz(-(mu + U / 2.) * dt, qubits[idxSite + nSites]);
            }
    } 

    operation ApplyCoulumbRepulsionTerms(nSites: Int, dt: Double, U: Double,
    qubits : Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0 .. nSites - 1) {
            let q1 = qubits[idxSite];
            let q2 = qubits[idxSite + nSites];
            CNOT(q1, q2);
            Rz((0.5 * U) * dt, q2);
            CNOT(q1, q2);
            }
    }

    operation ApplyRotationSequence(q1: Qubit, q2: Qubit, dt: Double, J: Double): Unit is Adj + Ctl {
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

    operation ApplyHoppingTermsChain(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
    }

    operation ApplySingleHoppingTerm(q1Idx: Int, q2Idx: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        let dist = q2Idx - q1Idx;
        if (dist > 1) {
            for (idx in 1 .. (dist - 3)) {
                CNOT(qubits[q1Idx + idx], qubits[q1Idx + idx + 1]);
            }
            CZ(qubits[q2Idx - 1], qubits[q2Idx]);
        }
        ApplyRotationSequence(qubits[q1Idx], qubits[q2Idx], dt, J);
        if (dist > 1) {
            CZ(qubits[q2Idx - 1], qubits[q2Idx]);
            for (idx in (dist - 3) .. 1) {
                CNOT(qubits[q1Idx + idx], qubits[q1Idx + idx + 1]);
            }
        }
    }

    operation ApplyRowHoppingTerms(q1Idxs: Int[], q2Idxs: Int[], ancillaIdx: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        let rungs = Length(q1Idxs);
        for (rung in 0 .. rungs - 1) {
            CNOT(qubits[q1Idxs[rung] + 1], qubits[ancillaIdx]);
            CNOT(qubits[q2Idxs[rung] - 1], qubits[ancillaIdx]);
            CZ(qubits[ancillaIdx], qubits[q2Idxs[rung]]);
        }
        for (rung in 0 .. (rungs - 1)) {
            let q1 = qubits[q1Idxs[rung]];
            let q2 = qubits[q2Idxs[rung]];
            ApplyRotationSequence(q1, q2, dt, J);
        }
        for (rung in (rungs - 1) .. 0) {
            CZ(qubits[ancillaIdx], qubits[q2Idxs[rung]]);
            CNOT(qubits[q2Idxs[rung] - 1], qubits[ancillaIdx]);
            CNOT(qubits[q1Idxs[rung] + 1], qubits[ancillaIdx]);
        }
    }

    operation ApplyHoppingTermsLadderVertical(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        for (idxSite in 0..4..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 3, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 3, dt, J, qubits);
        }
        for (idxSite in 2..4..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 3, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 3, dt, J, qubits);
        }
    }

    operation ApplyHoppingTermsLadderHorizontal(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        let length = nSites / 2;
        let upperRowIndsup = Reversed(SequenceI(0, length - 2));
        let upperRowIndsdown = Reversed(SequenceI(nSites, nSites + length - 2));
        let lowerRowIndsup = Reversed(SequenceI(length + 1, nSites - 1));
        let lowerRowIndsdown = Reversed(SequenceI(nSites + length + 1, 2 * nSites - 1));  
        ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, nSites, dt, J, qubits);
        ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, nSites + 1, dt, J, qubits);
    }

    operation ApplyHoppingTermsLattice(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        let length = Floor(Sqrt(IntAsDouble(nSites)));
        for (row in 0..2..(length - 2)){
            let upperRowIndsup = Reversed(SequenceI(row * length, length * (row + 1) - 2));
            let upperRowIndsdown = Reversed(SequenceI(nSites + row * length, nSites + length * (row + 1) - 2));
            let lowerRowIndsup = Reversed(SequenceI(length * (row + 1) + 1, length * (row + 2) - 1));
            let lowerRowIndsdown = Reversed(SequenceI(nSites + length * (row + 1) + 1, nSites + length * (row + 2) - 1));  
            let ancillaIdxup = 2 * nSites + row;
            let ancillaIdxdown = 2 * nSites + length - 1 + row;
            ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, ancillaIdxup, dt, J, qubits);
            ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, ancillaIdxdown, dt, J, qubits);
        }
        for (row in 1..2..(length - 2)){
            let upperRowIndsup = Reversed(SequenceI(row * length, length * (row + 1) - 2));
            let upperRowIndsdown = Reversed(SequenceI(nSites + row * length, nSites + length * (row + 1) - 2));
            let lowerRowIndsup = Reversed(SequenceI(length * (row + 1) + 1, length * (row + 2) - 1));
            let lowerRowIndsdown = Reversed(SequenceI(nSites + length * (row + 1) + 1, nSites + length * (row + 2) - 1));  
            let ancillaIdxup = 2 * nSites + row;
            let ancillaIdxdown = 2 * nSites + length - 1 + row;
            ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, ancillaIdxup, dt, J, qubits);
            ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, ancillaIdxdown, dt, J, qubits);
        }
    }

    operation EvolveSingleTimestep(nSites: Int, dt: Double, mu: Double, U: Double, J: Double, structure: Int, qubits: Qubit[]): Unit is Adj + Ctl {
        ApplyChemicalPotentialTerms(nSites, dt, mu, U, qubits);
        ApplyCoulumbRepulsionTerms(nSites, dt, U, qubits);
        if (structure == 1) {
            ApplyHoppingTermsChain(nSites, dt, J, qubits);
        } elif (structure == 2) {
            ApplyHoppingTermsLadderVertical(nSites, dt, J, qubits);
        } elif (structure == 3) {
            ApplyHoppingTermsLadderHorizontal(nSites, dt, J, qubits);
        } else{
            ApplyHoppingTermsLattice(nSites, dt, J, qubits);
        }
    }

    operation EvolveSingleTimestepDummy(nSites: Int, dt: Double, mu: Double, U: Double, J: Double, structure: Int): Unit is Adj + Ctl {
        using (qubits = Qubit[nSites * 2]) {
            EvolveSingleTimestep(nSites, dt, mu, U, J, structure, qubits);
        }
    }

    operation Evolve(initialState: Int[], time: Double, dt: Double, mu: Double, U: Double, J: Double, structure: Int): Int[] {
        let nSites = Length(initialState);
        if (structure == 3) {
            let qubitNum = 2 * nSites + 2;
        } elif (structure == 4) {
            let length = Floor(Sqrt(IntAsDouble(nSites)));
            let qubitNum = 2 * (nSites + length - 1);
        } else {
            let qubitNum = 2 * nSites;
        }
        using (qubits = Qubit[nSites * 2]) {
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
                    let newState = finalState[idxSite] + 2;
                    set finalState w/= idxSite <- newState;
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