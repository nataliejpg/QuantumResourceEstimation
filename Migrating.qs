namespace Quantum.Hubbard {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Measurement;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    // 

    /// # Summary
    /// 
    operation Juxtapose(index1: Int, index2: Int, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in index2..-1..index1+2) {
            SWAP(qubits[idxSite], qubits[idxSite - 1]);
        }
    }

    operation Separate(index1: Int, index2: Int, qubits: Qubit[]): Unit is Adj + Ctl {
        for (idxSite in index1+1..index2-1) {
            SWAP(qubits[idxSite], qubits[idxSite + 1]);
        }
    }

}