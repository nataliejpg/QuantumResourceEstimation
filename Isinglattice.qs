namespace Quantum.Isinglattice {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    //     H ≔ - J Σ'ᵢⱼ Zᵢ Zⱼ - hZ Σᵢ Zᵢ - hX Σᵢ Xᵢ

    operation EvolveSingleTimestep(nSites: Int, dt: Double, hx: Double[], hz: Double[], J: Double[][], qubits : Qubit[]): Unit is Adj + Ctl {
        for (idxSite in 0 .. nSites - 1) {
            Rx((-2.0 * hx[idxSite]) * dt, qubits[idxSite]);
            Rz((-2.0 * hz[idxSite]) * dt, qubits[idxSite]);
            }
        for (i in 0 .. nSites - 2) {
            for (j in i + 1 .. nSites - 1) {
                ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
            }
        }
    }

    operation EvolveSingleTimestepDummy(nSites: Int, dt: Double, hx: Double[], hz: Double[], J: Double[][]): Unit is Adj + Ctl {
        using (qubits = Qubit[nSites]) {
            EvolveSingleTimestep(nSites, dt, hx, hz, J, qubits)
        }
    }

    operation Evolve(initialState: Int[], time: Double, dt: Double, hx: Double[], hz: Double[], J: Double[][]): Result[] {
        let nSites = Length(initialState);
        using (qubits = Qubit[nSites]) {
            for (idxSite in 0 .. nSites - 1) {
                if (initialState[idxSite] == 1) {
                    X(qubits[idxSite]);
                }
            }
            let nSteps = Floor(time / dt);
            for (idxIter in 0 .. nSteps - 1) {
                EvolveSingleTimestep(nSites, dt, hx, hz, J, qubits);
            }
            return ForEach(MResetZ, qubits);
        }
    }
}