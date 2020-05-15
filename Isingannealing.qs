// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
namespace Quantum.Isingannealing {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    operation ApplyZZ(phi : Double, q1 : Qubit, q2 : Qubit) : Unit is Adj + Ctl {
        ApplyWithCA(CNOT(q1, _), Rz(-2.0 * phi, _), q2);
    }

    operation SimulateIsingEvolution(nSites : Int, time : Double, dt : Double, hx: Double, hz: Double, J: Double) : Result[] {
        using (qs = Qubit[nSites]) {
            ApplyToEach(H, qs);
            let nSteps = Floor(time / dt);
            for (idxIter in 0 .. nSteps - 1) {
                let sweepParameter = IntAsDouble(idxIter) / IntAsDouble(nSteps);
                for (idxSite in 0 .. nSites - 1) {

                    // Evolve under the transverse field for φx ≔ (1 - s) hx dt.
                    Rx(((-2.0 * (1.0 - sweepParameter)) * hx) * dt, qs[idxSite]);

                    // Evolve under the longitudinal field for φz ≔ s hz dt.
                    Rz(((-2.0 * sweepParameter) * hz) * dt, qs[idxSite]);

                    // If we aren't the last qubit, evolve under the Ising
                    // coupling for φJ ≔ s J dt.
                    if (idxSite < nSites - 2) {
                        ApplyZZ((sweepParameter * J) * dt, qs[idxSite], qs[idxSite + 1]);
                    }
                }
            }
            return ForEach(MResetZ, qs);
        }
    }

}