namespace Quantum.Hubbardchain {
    open Quantum.Hubbard;

    operation ApplySingleHoppingTerm(q1Idx: Int, q2Idx: Int, dt: Double, J: Double, qubits: Qubits): Unit is Adj + Ctl {
        let q1 = qubits[q1Idx]
        let q2 = qubits[q2Idx]
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
            ApplySingleHoppingTerm(idxSite, idxSite + 1, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 1, dt, J, qubits);
        }
        for (idxSite in 1..2..nSites - 2) {
            ApplySingleHoppingTerm(idxSite, idxSite + 1, dt, J, qubits);
            ApplySingleHoppingTerm(idxSite + nSites, idxSite + nSites + 1, dt, J, qubits);
        }
    }

}