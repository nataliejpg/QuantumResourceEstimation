namespace Quantum.Isinglongrange {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    //     H ≔ - J Σ'ᵢⱼ Zᵢ Zⱼ - hZ Σᵢ Zᵢ - hX Σᵢ Xᵢ

    operation ApplyZZ(phi : Double, q1 : Qubit, q2 : Qubit) : Unit is Adj + Ctl {
        ApplyWithCA(CNOT(q1, _), Rz(-2.0 * phi, _), q2);
    }

    operation ApplySingleQubitTerms(nSites: Int, dt: Double, hx: Double, hz: Double, qubits : Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0 .. nSites - 1) {
            Rx((-2.0 * hx) * dt, qubits[idxSite]);
            Rz((-2.0 * hz) * dt, qubits[idxSite]);
            }
    }

    operation ApplyTwoQubitTerms(nSites: Int, dt: Double, J: Double, dist: Int, qubits : Qubit[]): Unit is Adj + Ctl {
        for (step in 0 .. 1) {
            let nBatches = (nSites - (dist * step) / (2 * dist)) + ((nSites - (dist * step) % (2 * dist)) / (dist + 1));
            for (batch in 0 .. nBatches - 1) {
                for (i in 0 .. dist - 1) {
                    let start = (dist * step) + batch * 2 * dist + i;
                    if (start + dist < nSites) {
                        ApplyZZ(J * dt, qubits[start], qubits[start + dist]);
                    }
                }
            }
        }
    }

    operation ApplyTwoQubitTermsUnordered(nSites: Int, dt: Double, J: Double, dist: Int, qubits : Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0 .. nSites - (dist + 1)) {
            ApplyZZ(J * dt, qubits[idxSite], qubits[idxSite + dist]);
            }
    }

    operation EvolveSingleTimestep(nSites: Int, dt: Double, hx: Double, hz: Double, J: Double, scalings: Double[], ordered: Bool , qubits : Qubit[]): Unit is Adj + Ctl {
        let nCouplings = Length(scalings);
        ApplySingleQubitTerms(nSites, dt, hx, hz, qubits);
        for (dist in 1 .. nCouplings) {
            let scaledJ = J * scalings[dist - 1];
            if (ordered) {
                ApplyTwoQubitTerms(nSites, dt, scaledJ, dist, qubits);
            } else {
                ApplyTwoQubitTermsUnordered(nSites, dt, scaledJ, dist, qubits);
            }
        }
    }

    operation EvolveSingleTimestepDummy(nSites: Int, dt: Double, hx: Double, hz: Double, J: Double, scalings: Double[], ordered: Bool): Unit is Adj + Ctl {
        let nCouplings = Length(scalings);
        using (qubits = Qubit[nSites]) {
            ApplySingleQubitTerms(nSites, dt, hx, hz, qubits);
            for (dist in 1 .. nCouplings) {
                let scaledJ = J * scalings[dist - 1];
                if (ordered) {
                    ApplyTwoQubitTerms(nSites, dt, scaledJ, dist, qubits);
                } else {
                    ApplyTwoQubitTermsUnordered(nSites, dt, scaledJ, dist, qubits);
                }
            }
        }
    }

    operation Evolve(initialState: Int[], time: Double, dt: Double, hx: Double, hz: Double, J: Double, scalings: Double[], ordered: Bool): Result[] {
        let nSites = Length(initialState);
        using (qubits = Qubit[nSites]) {
            for (idxSite in 0 .. nSites - 1) {
                if (initialState[idxSite] == 1) {
                    X(qubits[idxSite]);
                }
            }
            let nSteps = Floor(time / dt);
            for (idxIter in 0 .. nSteps - 1) {
                EvolveSingleTimestep(nSites, dt, hx, hz, J, scalings, ordered, qubits);
            }
            return ForEach(MResetZ, qubits);
        }
    }
}
