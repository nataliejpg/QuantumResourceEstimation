namespace Quantum.Hubbardladder {
    open Quantum.Hubbard;

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

    operation ApplyRowHoppingTerms(q1Idxs: Int[], q2Idxs: Int[], ancillaIdx: Int, dt: Double, J: Double, qubits: Qubits[]): Unit is Adj + Ctl {
        let rungs = Length(q1Idxs);
        for rung in (0 .. rungs - 1) {
            CNOT(qubits[q1Idxs[rung] + 1], qubits[ancillaIdx]);
            CNOT(qubits[q2Idxs[rung] - 1], qubits[ancillaIdx]);
            CZ(qubits[ancillaIdx], qubits[q2Idxs[rung]]);
        }
        for rung in (0 .. (rungs - 1) {
            let q1 = qubits[q1Idxs[rung]];
            let q2 = qubits[q2Idxs[rung]];
            ApplyRotationSequence(q1, q2, dt, J)
        }
        for rung in ((rungs - 1) .. 0) {
            CZ(qubits[ancillaIdx], qubits[q2Idxs[rung]]);
            CNOT(qubits[q2Idxs[rung] - 1], qubits[ancillaIdx]);
            CNOT(qubits[q1Idxs[rung] + 1], qubits[ancillaIdx]);
        }
    }

    operation ApplyHoppingTerms(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplyRotationSequence(qubits[idxSite], qubits[idxSite + 1], dt, J);
            ApplyRotationSequence(qubits[idxSite + nSites], qubits[idxSite + nSites + 1], dt, J);
        }
        let length = Floor(Sqrt(IntasDouble(N)));
        mutable upperRowIndsup = Int[length - 1];
        mutable lowerRowIndsup = Int[length - 1];
        mutable upperRowIndsdown = Int[length - 1];
        mutable lowerRowIndsdown = Int[length - 1];
        for row in (0..2..(length - 2)){
            for (idxSite in 0 .. length - 2) {
                set upperRowIndsup w/= length - 1 - idxSite <- idxSite + (row * length);
                set lowerRowIndsup w/= length - 1 - idxSite <- nSites - idxSite - 1 + (row * length);
                set upperRowIndsdown w/= length - 1 - idxSite <- idxSite + (row * length) + nSites;
                set lowerRowIndsdown w/= length - 1 - idxSite <- 2 * nSites - idxSite - 1 + (row * length);
            }
            let ancillaIdxup = 2 * nSites + row;
            let ancillaIdxdown = 2 * nSites + length - 1 + row;
            ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, ancillaIdxup, dt, J, qubits);
            ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, ancillaIdxdown, dt, J, qubits);
        }
        for row in (1..2..(length - 2)){
            for (idxSite in 0 .. length - 2) {
                set upperRowIndsup w/= idxSite <- idxSite + (row * length);
                set lowerRowIndsup w/= idxSite <- nSites - idxSite - 1 + (row * length);
                set upperRowIndsdown w/= idxSite <- idxSite + (row * length) + nSites;
                set lowerRowIndsdown w/= idxSite <- 2 * nSites - idxSite - 1 + (row * length);
            }
            let ancillaIdxup = 2 * nSites + row;
            let ancillaIdxdown = 2 * nSites + length - 1 + row;
            ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, ancillaIdxup, dt, J, qubits);
            ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, ancillaIdxdown, dt, J, qubits);
        }
    }
}