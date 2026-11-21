from __future__ import annotations
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.quantum_info import Statevector
from src.q_loader import build_dna_loader
import math
import pytest

# Test on 5 real DNA windows
def test_loader_exact_superposition_and_padding() -> None:
    windows_bin = ["0000", "0110", "1101", "1010", "0011", ]
    n_idx = 3
    n_data = 4
    valid_len = len(windows_bin)

    loader_gate = build_dna_loader(windows_bin, n_idx=n_idx, n_data=n_data, valid_len=valid_len)

    idx_reg = QuantumRegister(n_idx, "idx")
    data_reg = QuantumRegister(n_data, "data")
    valid_reg = QuantumRegister(1, "valid")
    qc = QuantumCircuit(idx_reg, data_reg, valid_reg)

    qc.h(idx_reg)
    qc.append(loader_gate, qc.qubits)

    sv = Statevector.from_instruction(qc)
    probs = sv.probabilities_dict(decimals=12)

    expected_prob = 1.0 / (1 << n_idx)
    assert len(probs) == (1 << n_idx)

    for bitstream, prob in probs.items():
        assert math.isclose(prob, expected_prob, abs_tol=1e-10)

        valid_bit = bitstream[0]
        data_bits = bitstream[1: 1 + n_data][::-1]
        idx_val = int(bitstream[1 + n_data:], 2)

        if idx_val < valid_len:
            assert valid_bit == "1", f"Indice valido {idx_val} dovrebbe avere valid=1"
            assert data_bits == windows_bin[idx_val], f"Mismatch all'indice {idx_val}: atteso {windows_bin[idx_val]}, ottenuto {data_bits}"
        else:
            assert valid_bit == "0", f"Indice padding {idx_val} dovrebbe avere valid=0"
            assert data_bits == "0" * n_data, f"Indice padding {idx_val} dovrebbe avere data=0...0"

# Test that uncomputation erase entanglement and restore the initial state
def test_loader_uncomputation_is_exact_identity() -> None:
    windows_bin = ["0101", "1110", "0011"]
    n_idx = 2
    n_data = 4

    loader_gate = build_dna_loader(windows_bin, n_idx=n_idx, n_data=n_data, valid_len=3)

    idx_reg = QuantumRegister(n_idx, "idx")
    data_reg = QuantumRegister(n_data, "data")
    valid_reg = QuantumRegister(1, "valid")
    qc = QuantumCircuit(idx_reg, data_reg, valid_reg)

    qc.h(idx_reg)
    qc.append(loader_gate, qc.qubits)
    qc.append(loader_gate.inverse(), qc.qubits)

    sv_after = Statevector.from_instruction(qc)

    qc_ref = QuantumCircuit(idx_reg, data_reg, valid_reg)
    qc_ref.h(idx_reg)
    sv_ref = Statevector.from_instruction(qc_ref)

    assert sv_after.equiv(sv_ref), "L'uncomputation non ha ripristinato lo stato iniziale!"

# Test that the loader ignores padded entries beyond valid length
def test_loader_ignores_padded_entries_beyond_valid_len() -> None:
    windows_bin = ["11", "01", "11", "11"]
    loader_gate = build_dna_loader(windows_bin, n_idx=2, n_data=2, valid_len=2)

    idx_reg = QuantumRegister(2, "idx")
    data_reg = QuantumRegister(2, "data")
    valid_reg = QuantumRegister(1, "valid")
    qc = QuantumCircuit(idx_reg, data_reg, valid_reg)

    qc.h(idx_reg)
    qc.append(loader_gate, qc.qubits)

    probs = Statevector.from_instruction(qc).probabilities_dict(decimals=12)
    for bitstream in probs:
        valid_bit = bitstream[0]
        data_bits = bitstream[1:3][::-1]
        idx_val = int(bitstream[3:], 2)
        if idx_val >= 2:
            assert valid_bit == "0"
            assert data_bits == "00"

# Test input validation
@pytest.mark.parametrize(
    "windows, n_idx, n_data, valid_len",
    [
        (["01"], 0, 2, 1),
        (["01"], 1, 0, 1),
        (["01", "10", "11"], 1, 2, 3),
        (["01"], 2, 2, 2),
        (["02"], 1, 2, 1),
        (["010"], 1, 2, 1),
    ],
)
def test_loader_invalid_inputs_raise_value_error(windows: list[str], n_idx: int, n_data: int, valid_len: int) -> None:
    with pytest.raises(ValueError):
        build_dna_loader(windows, n_idx=n_idx, n_data=n_data, valid_len=valid_len)
