import time
from typing import Any

from qiskit import QuantumCircuit
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeBrisbane


def get_backend(use_noise: bool = False) -> AerSimulator:
    """
    Restituisce il backend AerSimulator configurato in modalità ideale
    oppure rumorosa (emulando l'hardware reale IBM FakeBrisbane).
    """
    if use_noise:
        fake_hw = FakeBrisbane()
        return AerSimulator.from_backend(fake_hw)
    return AerSimulator(basis_gates=["rz", "sx", "x", "cx"])


def count_two_qubit_gates(qc: QuantumCircuit) -> int:
    """
    Conta il numero totale di gate a 2 qubit nel circuito transpilato.
    Nota: i chip IBM Eagle (come FakeBrisbane) usano 'ecr' come gate nativo
    a 2 qubit al posto di 'cx', quindi contiamo qualsiasi istruzione su 2 qubit.
    """
    return sum(1 for instruction in qc.data if instruction.operation.num_qubits == 2 and instruction.operation.name not in ("barrier", "measure"))


def decode_results(counts: dict[str, int], threshold: float | None = None, ) -> tuple[list[int], dict[int, float]]:
    """
    Decodifica i conteggi di Qiskit (little-endian) in indici interi del DNA
    e filtra il rumore di fondo tramite una soglia statistica.

    Args:
        counts: Dizionario {bitstring: occorrenze} restituito da Qiskit.
        threshold: Soglia minima di probabilità (0.0 - 1.0) per considerare
                   un indice come match valido. Se None, usa una soglia dinamica
                   pari a max(0.10, 2.0 / numero_stati_possibili).

    Returns:
        Tupla (found_indices, probabilities):
        - found_indices: lista ordinata degli indici sopra soglia.
        - probabilities: dizionario {indice_intero: probabilità}.
    """
    total_shots = sum(counts.values())
    if total_shots == 0:
        return [], {}

    probabilities: dict[int, float] = {}
    n_bits = 0

    for raw_bitstring, count in counts.items():
        # Rimuove eventuali spazi nel caso ci siano più registri classici
        clean_bits = raw_bitstring.replace(" ", "")
        n_bits = max(n_bits, len(clean_bits))

        # Qiskit usa little-endian (q0 è il bit più a destra):
        # '001' -> int('001', 2) == 1
        idx = int(clean_bits, 2)
        probabilities[idx] = probabilities.get(idx, 0.0) + (count / total_shots)

    # Soglia statistica di default se non specificata:
    # con distribuzione uniforme ogni stato avrebbe probabilità 1 / 2^n_bits.
    # Richiediamo almeno il doppio della probabilità uniforme (con minimo 0.10).
    if threshold is None:
        num_states = 2 ** n_bits if n_bits > 0 else 1
        threshold = max(0.10, 2.0 / num_states)

    found_indices = sorted(idx for idx, prob in probabilities.items() if prob >= threshold)

    return found_indices, probabilities


def run_circuit(grover_qc: QuantumCircuit, shots: int = 1024, use_noise: bool = False, threshold: float | None = None, seed: int | None = 42, ) -> dict[str, Any]:
    """
    Transpila (optimization_level=3) ed esegue il circuito di Grover su
    AerSimulator (ideale o rumoroso), restituendo il contratto ExecutionResult.

    Args:
        grover_qc: Circuito quantistico con misure sui qubit di indice.
        shots: Numero di esecuzioni (campionamenti).
        use_noise: Se True, usa il modello di rumore di FakeBrisbane.
        threshold: Soglia di probabilità per filtrare il rumore.
        seed: Seed per la riproducibilità del simulatore e del transpiler.

    Returns:
        Dizionario ExecutionResult:
        {
            "counts": dict[str, int],
            "found_indices": list[int],
            "probabilities": dict[int, float],
            "depth": int,
            "cx_count": int,
            "exec_time_ms": float
        }
    """
    backend = get_backend(use_noise=use_noise)

    # Transpilazione ottimizzata al massimo livello (level 3)
    pass_manager = generate_preset_pass_manager(optimization_level=3, backend=backend, seed_transpiler=seed)
    transpiled_qc = pass_manager.run(grover_qc)

    # Estrazione metriche hardware dal circuito transpilato
    depth = transpiled_qc.depth()
    cx_count = count_two_qubit_gates(transpiled_qc)

    # Esecuzione e misurazione del tempo in millisecondi
    start_time = time.perf_counter()
    job = backend.run(transpiled_qc, shots=shots, seed_simulator=seed)
    result = job.result()
    exec_time_ms = (time.perf_counter() - start_time) * 1000.0

    counts: dict[str, int] = result.get_counts()
    found_indices, probabilities = decode_results(counts, threshold=threshold)

    return {
        "counts": counts, "found_indices": found_indices, "probabilities": probabilities, "depth": depth, "cx_count": cx_count, "exec_time_ms": round(exec_time_ms, 3),
    }
