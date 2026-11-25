import time
from typing import Any

from qiskit import QuantumCircuit
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_aer import AerSimulator
from qiskit_ibm_runtime.fake_provider import FakeBrisbane


def get_backend(use_noise: bool = False) -> AerSimulator:
    if use_noise:
        fake_hw = FakeBrisbane()
        return AerSimulator.from_backend(fake_hw)
    return AerSimulator(basis_gates=["rz", "sx", "x", "cx"])


def count_two_qubit_gates(qc: QuantumCircuit) -> int:
    return sum(1 for instruction in qc.data if instruction.operation.num_qubits == 2 and instruction.operation.name not in ("barrier", "measure"))


def decode_results(counts: dict[str, int], threshold: float | None = None, ) -> tuple[list[int], dict[int, float]]:
    total_shots = sum(counts.values())
    if total_shots == 0:
        return [], {}

    probabilities: dict[int, float] = {}
    n_bits = 0

    for raw_bitstring, count in counts.items():
        clean_bits = raw_bitstring.replace(" ", "")
        n_bits = max(n_bits, len(clean_bits))

        idx = int(clean_bits, 2)
        probabilities[idx] = probabilities.get(idx, 0.0) + (count / total_shots)

    if threshold is None:
        num_states = 2 ** n_bits if n_bits > 0 else 1
        threshold = max(0.10, 2.0 / num_states)

    found_indices = sorted(idx for idx, prob in probabilities.items() if prob >= threshold)

    return found_indices, probabilities


def run_circuit(grover_qc: QuantumCircuit, shots: int = 1024, use_noise: bool = False, threshold: float | None = None, seed: int | None = 42, ) -> dict[str, Any]:
    backend = get_backend(use_noise=use_noise)

    pass_manager = generate_preset_pass_manager(optimization_level=3, backend=backend, seed_transpiler=seed)
    transpiled_qc = pass_manager.run(grover_qc)

    depth = transpiled_qc.depth()
    cx_count = count_two_qubit_gates(transpiled_qc)

    start_time = time.perf_counter()
    job = backend.run(transpiled_qc, shots=shots, seed_simulator=seed)
    result = job.result()
    exec_time_ms = (time.perf_counter() - start_time) * 1000.0

    counts: dict[str, int] = result.get_counts()
    found_indices, probabilities = decode_results(counts, threshold=threshold)

    return {
        "counts": counts, "found_indices": found_indices, "probabilities": probabilities, "depth": depth, "cx_count": cx_count, "exec_time_ms": round(exec_time_ms, 3),
    }
