from __future__ import annotations
from qiskit.quantum_info import Statevector
from src.dna_core import process_dna
from src.q_grover import build_grover_circuit
import math

def _extract_index_probabilities(probs: dict[str, float], n_data: int) -> dict[int, float]:
    idx_probs: dict[int, float] = {}
    for bitstream, prob in probs.items():
        assert bitstream[0] == "0"
        assert bitstream[1: 1 + n_data] == "0" * n_data
        idx_val = int(bitstream[1 + n_data:], 2)
        idx_probs[idx_val] = prob
    return idx_probs

def test_dna_to_grover_single_match() -> None:
    payload = process_dna("ACGTTAG", "TT")
    classical_indices = payload["classical_res"]["indices"]
    assert classical_indices == [3]

    qc, k = build_grover_circuit(payload, num_solutions=len(classical_indices), measure=False)
    assert k == 2

    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    idx_probs = _extract_index_probabilities(probs, payload["n_data"])

    assert idx_probs[3] > 0.94
    for idx, prob in idx_probs.items():
        if idx != 3:
            assert prob < 0.02

def test_dna_to_grover_aa_pattern_ignores_dna_core_padding() -> None:
    payload = process_dna("CGAATC", "AA")
    classical_indices = payload["classical_res"]["indices"]
    assert classical_indices == [2]
    assert payload["valid_len"] == 5

    qc, k = build_grover_circuit(payload, num_solutions=len(classical_indices), measure=False)
    assert k == 2

    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    idx_probs = _extract_index_probabilities(probs, payload["n_data"])

    assert idx_probs[2] > 0.94
    for padded_idx in (5, 6, 7):
        assert idx_probs[padded_idx] < 0.02

def test_dna_to_grover_multiple_matches() -> None:
    payload = process_dna("ACGTACGT", "CG")
    classical_indices = payload["classical_res"]["indices"]
    assert classical_indices == [1, 5]

    qc, k = build_grover_circuit(payload, num_solutions=len(classical_indices), measure=False)
    assert k == 1

    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    idx_probs = _extract_index_probabilities(probs, payload["n_data"])

    assert math.isclose(idx_probs[1], 0.5, abs_tol=1e-10)
    assert math.isclose(idx_probs[5], 0.5, abs_tol=1e-10)
    assert math.isclose(idx_probs[1] + idx_probs[5], 1.0, abs_tol=1e-10)

def test_dna_to_grover_absent_pattern_stays_uniform() -> None:
    payload = process_dna("ACGTAC", "TT")
    assert payload["classical_res"]["indices"] == []

    qc, _ = build_grover_circuit(payload, num_solutions=1, measure=False)
    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    idx_probs = _extract_index_probabilities(probs, payload["n_data"])

    expected_uniform = 1.0 / (1 << payload["n_idx"])
    for prob in idx_probs.values():
        assert math.isclose(prob, expected_uniform, abs_tol=1e-10)