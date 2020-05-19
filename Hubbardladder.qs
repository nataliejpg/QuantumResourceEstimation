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

    operation ApplySingleHoppingTerm(q1Idx: Int, q2Idx: Int, dt: Double, J: Double, qubits: Qubits[]): Unit is Adj + Ctl {
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

    operation ApplyHoppingTermsLadderVertical(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0..2..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 1, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 1, dt, J, qubits);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 1, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 1, dt, J, qubits);
        }
        for (idxSite in 0..4..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 3, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 3, dt, J, qubits);
        }
        for (idxSite in 2..4..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 3, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 3, dt, J), qubits;
        }
    }

    operation ApplyHoppingTermsLadderHorizontal(nSites: Int, dt: Double, J: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0..2..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 1, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 1, dt, J, qubits);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 1, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 1, dt, J, qubits);
        }
        let length = nSites / 2;
        mutable upperRowIndsup = Int[length - 1];
        mutable lowerRowIndsup = Int[length - 1];
        mutable upperRowIndsdown = Int[length - 1];
        mutable lowerRowIndsdown = Int[length - 1];
        for (idxSite in 0 .. length - 2) {
            set upperRowIndsup w/= length - 1 - idxSite <- idxSite;
            set lowerRowIndsup w/= length - 1 - idxSite <- nSites - idxSite - 1;
            set upperRowIndsdown w/= length - 1 - idxSite <- idxSite + nSites;
            set lowerRowIndsdown w/= length - 1 - idxSite <- 2 * nSites - idxSite - 1;
        }
        ApplyRowHoppingTerms(upperRowIndsup, lowerRowIndsup, nSites, dt, J, qubits);
        ApplyRowHoppingTerms(upperRowIndsdown, lowerRowIndsdown, nSites + 1, dt, J, qubits);
    }

}