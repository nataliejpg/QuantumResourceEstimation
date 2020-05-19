namespace Quantum.Isinglongrange {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    //     H ≔ - Σ'ᵢⱼ Jᵢⱼ Zᵢ Zⱼ - Σᵢ hᵢ Zᵢ - Σᵢ gᵢ Xᵢ

    operation EvolveCouplings(nSites: Int, dt: Double, J: Double[][], qubits: Qubit[]): Unit is Adj + Ctl {
        for (i in 0 .. nSites - 2) {
            for (j in i + 1 .. nSites - 1) {
                ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
            }
        }
    } 

    operation EvolveCouplingsNested(nSites: Int, dt: Double, J: Double[][], qubits: Qubit[]): Unit is Adj + Ctl {
        let num = nSites/2 + nSites%2;
        for step in (0 .. num) {
            for ind in (0 .. nSites/2 - 1) {
                let i = ind - step;
                let j = nSites - (ind + step + 1);
                ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
            }
            if ((nSites%2 == 0) or (step < nSites/2)){
                for ind in (0 .. nSites/2 - 1) {
                    let i = ind - step;
                    if (ind == (num - 1)) {
                        let j = nSites - (step + 1);
                    } else {
                        let j = nSites - (step + ind + 2);
                    }
                    ApplyWithCA(CNOT(qubits[i], _), Rz(-2.0 * J[i][j] * dt, _), qubits[j]);
                }
            }
        }
    }

    operation EvolveSingleTimestep(nSites: Int, dt: Double, g: Double[], h: Double[], J: Double[][], qubits : Qubit[], nested: Bool): Unit is Adj + Ctl {
        for (idxSite in 0 .. nSites - 1) {
            Rx((-2.0 * g[idxSite]) * dt, qubits[idxSite]);
            Rz((-2.0 * h[idxSite]) * dt, qubits[idxSite]);
            }
        if nested {
            EvolveCouplingsNested(nSites, dt, J, qubits);
        } else {
            EvolveCouplings(nSites, dt, J, qubits)
        }
    }

    operation EvolveSingleTimestepDummy(nSites: Int, dt: Double, g: Double[], h: Double[], J: Double[][]): Unit is Adj + Ctl {
        using (qubits = Qubit[nSites]) {
            EvolveSingleTimestep(nSites, dt, g, h, J, qubits)
        }
    }

    operation Evolve(initialState: Int[], time: Double, dt: Double, g: Double[], h: Double[], J: Double[][], nested: Bool): Result[] {
        let nSites = Length(initialState);
        using (qubits = Qubit[nSites]) {
            for (idxSite in 0 .. nSites - 1) {
                if (initialState[idxSite] == 1) {
                    X(qubits[idxSite]);
                }
            }
            let nSteps = Floor(time / dt);
            for (idxIter in 0 .. nSteps - 1) {
                EvolveSingleTimestep(nSites, dt, g, h, J, qubits, nested);
            }
            return ForEach(MResetZ, qubits);
        }
    }
}