from __future__ import annotations
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import Statevector
from src.q_grover import (build_diffuser, build_grover_circuit, build_pattern_oracle, compute_optimal_iterations)
from src.q_loader import build_dna_loader
import math
import pytest

# Test optimal iterations
@pytest.mark.parametrize(
    "n_idx, num_solutions, expected",
    [
        (2, 1, 1),
        (3, 1, 2),
        (3, 2, 1),
        (4, 1, 3),
        (2, 4, 0),
    ],
)
def test_compute_optimal_iterations(n_idx: int, num_solutions: int, expected: int) -> None:
    assert compute_optimal_iterations(n_idx, num_solutions) == expected

# Test optimal iterations with invalid inputs
@pytest.mark.parametrize(
    "n_idx, num_solutions",
    [
        (0, 1),
        (2, 0),
        (2, 5),
    ],
)
def test_compute_optimal_iterations_invalid(n_idx: int, num_solutions: int) -> None:
    with pytest.raises(ValueError):
        compute_optimal_iterations(n_idx, num_solutions)

# Test phase kickback and uncomputation
def test_pattern_oracle_phase_kickback_and_uncomputation() -> None:
    windows_bin = ["0000", "0110", "1101", "1010", "0011"]
    pattern_bin = "0000"
    n_idx = 3
    n_data = 4
    valid_len = len(windows_bin)

    loader_gate = build_dna_loader(windows_bin, n_idx=n_idx, n_data=n_data, valid_len=valid_len)
    oracle_gate = build_pattern_oracle(loader_gate, pattern_bin, n_idx=n_idx, n_data=n_data)

    idx_reg = QuantumRegister(n_idx, "idx")
    data_reg = QuantumRegister(n_data, "data")
    valid_reg = QuantumRegister(1, "valid")
    qc = QuantumCircuit(idx_reg, data_reg, valid_reg)

    qc.h(idx_reg)
    qc.append(oracle_gate, qc.qubits)

    sv = Statevector.from_instruction(qc)
    expected_amp = 1.0 / math.sqrt(1 << n_idx)

    for state_int, amp in enumerate(sv.data):
        if abs(amp) < 1e-10:
            continue
        bitstream = format(state_int, f"0{n_idx + n_data + 1}b")
        valid_bit = bitstream[0]
        data_bits = bitstream[1: 1 + n_data]
        idx_val = int(bitstream[1 + n_data:], 2)

        assert valid_bit == "0"
        assert data_bits == "0" * n_data

        if idx_val == 0:
            assert math.isclose(amp.real, -expected_amp, abs_tol=1e-10)
        else:
            assert math.isclose(amp.real, expected_amp, abs_tol=1e-10)
        assert math.isclose(amp.imag, 0.0, abs_tol=1e-10)

# Test single match
def test_grover_single_match_high_probability() -> None:
    payload = {
        "windows_bin": ["0001", "0110", "1011", "1100", "0010", "0101"],
        "pattern_bin": "1011",
        "n_idx": 3,
        "n_data": 4,
        "valid_len": 6,
    }
    qc, k = build_grover_circuit(payload, num_solutions=1, measure=False)
    assert k == 2

    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    idx_probs: dict[int, float] = {}
    for bitstream, prob in probs.items():
        assert bitstream[0] == "0"
        assert bitstream[1:5] == "0000"
        idx_val = int(bitstream[5:], 2)
        idx_probs[idx_val] = prob

    assert idx_probs[2] > 0.94
    for i in range(8):
        if i != 2:
            assert idx_probs[i] < 0.02

# Test multiple matches
def test_grover_multiple_matches_exact_amplification() -> None:
    payload = {
        "windows_bin": ["0011", "1010", "0101", "1010", "1111", "0001"],
        "pattern_bin": "1010",
        "n_idx": 3,
        "n_data": 4,
        "valid_len": 6,
    }
    qc, k = build_grover_circuit(payload, num_solutions=2, measure=False)
    assert k == 1

    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    idx_probs = {int(bitstream[5:], 2): prob for bitstream, prob in probs.items()}

    assert math.isclose(idx_probs[1], 0.5, abs_tol=1e-10)
    assert math.isclose(idx_probs[3], 0.5, abs_tol=1e-10)
    assert math.isclose(idx_probs[1] + idx_probs[3], 1.0, abs_tol=1e-10)

# Test AA ignoring padding
def test_grover_aa_pattern_ignores_padding() -> None:
    payload = {
        "windows_bin": ["0001", "0110", "1011", "1100", "0010", "0101"], "pattern_bin": "1011", "n_idx": 3, "n_data": 4, "valid_len": 6,
    }
    qc, k = build_grover_circuit(payload, num_solutions=1, measure=False)
    assert k == 2

    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    idx_probs = {int(bitstream[5:], 2): prob for bitstream, prob in probs.items()}

    assert idx_probs[2] > 0.94
    for padded_idx in (5, 6, 7):
        assert idx_probs[padded_idx] < 0.02

# Test circuit with measurements
def test_grover_circuit_with_measurements_contract() -> None:
    payload = {
        "windows_bin": ["01", "10", "11"],
        "pattern_bin": "10",
        "n_idx": 2,
        "n_data": 2,
        "valid_len": 3,
    }

    qc, k = build_grover_circuit(payload, num_solutions=1, measure=True)
    assert isinstance(qc, QuantumCircuit)
    assert k == 1
    assert qc.num_qubits == 2 + 2 + 1
    assert qc.num_clbits == 2

    measure_ops = [inst for inst in qc.data if inst.operation.name == "measure"]
    assert len(measure_ops) == 2

def test_diffuser_invalid_n_idx() -> None:
    with pytest.raises(ValueError):
        build_diffuser(0)
