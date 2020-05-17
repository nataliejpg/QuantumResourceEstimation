namespace Quantum.Hubbardevolve {
    open Quantum.Hubbard;
    open Quantum.Hubbardchain;
    open Quantum.Hubbardladder;

    operation EvolveSingleTimestep(nSites: Int, dt: Double, mu: Double, U: Double, J: Double, structure: Int, qubits: Qubit[]): Unit is Adj + Ctl {
        ApplyChemicalPotentialTerms(nSites, dt, mu, U, qubits);
        ApplyCoulumbRepulsionTerms(nSites, dt, U, qubits);
        if (structure == 1) {
            ApplyHoppingTermsChain(nSites, dt, J, qubits);
        } else-if (structure == 2) {
            ApplyHoppingTermsLadderVertical(nSites, dt, J, qubits);
        } else-if (structure == 3) {
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