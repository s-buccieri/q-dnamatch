from __future__ import annotations
from src.dna_core import process_dna
from src.q_grover import build_grover_circuit
from src.q_runtime import run_circuit
import pytest
import math

@pytest.mark.parametrize("sequence, pattern", [
    ("ACGTTAGC", "TT"), ("CGAATC", "AA"), ("ACGTACGT", "CG"),
], )
def test_end_to_end_pipeline_ideal_monitored(sequence: str, pattern: str) -> None:
    print(f"\n{'=' * 68}")
    print(f"[END-TO-END IDEAL] Sequenza: {sequence} | Pattern: {pattern}")

    payload = process_dna(sequence, pattern)
    classical_indices = payload["classical_res"]["indices"]
    classical_queries = payload["classical_res"]["queries"]
    print(f"[1. dna_core] valid_len={payload['valid_len']} | "
          f"n_idx={payload['n_idx']} | n_data={payload['n_data']}")
    print(f"[1. dna_core] Baseline classica -> indici={classical_indices} | "
          f"query O(N)={classical_queries}")

    qc, k = build_grover_circuit(payload, num_solutions=len(classical_indices), measure=True)
    print(f"[2. quantum_core] qubit={qc.num_qubits} | clbit={qc.num_clbits} | "
          f"iterazioni Grover O(sqrt(N))={k}")

    res = run_circuit(qc, shots=1024, use_noise=False, seed=42)
    sorted_probs = dict(sorted(res["probabilities"].items()))
    print(f"[3. q_runtime] found_indices={res['found_indices']} | "
          f"depth={res['depth']} | 2q_gates={res['cx_count']} | "
          f"time={res['exec_time_ms']} ms")
    print(f"[3. q_runtime] probabilità decodificate={sorted_probs}")

    assert res["found_indices"] == classical_indices
    assert math.isclose(sum(res["probabilities"].values()), 1.0, abs_tol=1e-9)
    assert res["depth"] > 0
    assert res["cx_count"] > 0


def test_end_to_end_pipeline_ideal_vs_noisy_monitored() -> None:
    sequence = "ACGT"
    pattern = "G"

    print(f"\n{'=' * 68}")
    print(f"[END-TO-END IDEAL vs NOISY] Sequenza: {sequence} | Pattern: {pattern}")

    payload = process_dna(sequence, pattern)
    classical_indices = payload["classical_res"]["indices"]
    qc, k = build_grover_circuit(payload, num_solutions=len(classical_indices), measure=True)

    res_ideal = run_circuit(qc, shots=1024, use_noise=False, seed=42)
    res_noisy = run_circuit(qc, shots=1024, use_noise=True, seed=42)

    target = classical_indices[0]
    p_ideal = res_ideal["probabilities"].get(target, 0.0)
    p_noisy = res_noisy["probabilities"].get(target, 0.0)

    print(f"[IDEAL] found={res_ideal['found_indices']} | P(idx={target})={p_ideal:.4f} | "
          f"depth={res_ideal['depth']} | 2q_gates={res_ideal['cx_count']} | "
          f"time={res_ideal['exec_time_ms']} ms")
    print(f"[NOISY] found={res_noisy['found_indices']} | P(idx={target})={p_noisy:.4f} | "
          f"depth={res_noisy['depth']} | 2q_gates={res_noisy['cx_count']} | "
          f"time={res_noisy['exec_time_ms']} ms")

    assert res_ideal["found_indices"] == classical_indices
    assert p_ideal > 0.95
    assert sum(res_noisy["counts"].values()) == 1024
    assert math.isclose(sum(res_noisy["probabilities"].values()), 1.0, abs_tol=1e-9)
    assert res_noisy["depth"] > 0
    assert res_noisy["cx_count"] > 0