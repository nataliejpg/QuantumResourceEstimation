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

}