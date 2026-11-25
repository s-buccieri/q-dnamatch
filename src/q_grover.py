from __future__ import annotations
from typing import Any
from qiskit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.circuit import Gate
from src.q_loader import build_dna_loader
import math

def compute_optimal_iterations(n_idx: int, num_solutions: int = 1) -> int:
    if n_idx < 1:
        raise ValueError(f"n_idx deve essere >= 1, ricevuto {n_idx}.")

    total_states = 1 << n_idx
    if num_solutions < 1 or num_solutions > total_states:
        raise ValueError(f"num_solutions ({num_solutions}) deve essere compreso tra 1 e {total_states}.")

    if num_solutions == total_states:
        return 0
    return math.floor((math.pi / 4.0) * math.sqrt(total_states / num_solutions))

def build_pattern_oracle(loader_gate: Gate, pattern_bin: str, n_idx: int, n_data: int, ) -> Gate:
    if n_idx < 1:
        raise ValueError(f"n_idx deve essere >= 1, ricevuto {n_idx}.")

    if n_data < 1:
        raise ValueError(f"n_data deve essere >= 1, ricevuto {n_data}.")

    if len(pattern_bin) != n_data or any(ch not in ("0", "1") for ch in pattern_bin):
        raise ValueError(f"pattern_bin non valido: '{pattern_bin}'. Attesa stringa binaria di lunghezza {n_data}.")

    expected_qubits = n_idx + n_data + 1
    if loader_gate.num_qubits != expected_qubits:
        raise ValueError(f"loader_gate opera su {loader_gate.num_qubits} qubit, attesi {expected_qubits}.")

    idx_reg = QuantumRegister(n_idx, "idx")
    data_reg = QuantumRegister(n_data, "data")
    valid_reg = QuantumRegister(1, "valid")
    qc = QuantumCircuit(idx_reg, data_reg, valid_reg, name="Pattern_Oracle")

    qc.append(loader_gate, qc.qubits)

    zero_bits = [data_reg[pos] for pos, bit in enumerate(pattern_bin) if bit == "0"]
    for q in zero_bits:
        qc.x(q)

    data_controls = [data_reg[i] for i in range(n_data)]
    qc.h(valid_reg[0])
    qc.mcx(data_controls, valid_reg[0])
    qc.h(valid_reg[0])

    for q in zero_bits:
        qc.x(q)

    qc.append(loader_gate.inverse(), qc.qubits)

    return qc.to_gate(label="Pattern_Oracle")

def build_diffuser(n_idx: int) -> Gate:
    if n_idx < 1:
        raise ValueError(f"n_idx deve essere >= 1, ricevuto {n_idx}.")

    idx_reg = QuantumRegister(n_idx, "idx")
    qc = QuantumCircuit(idx_reg, name="Diffuser")

    qc.h(idx_reg)
    qc.x(idx_reg)

    if n_idx == 1:
        qc.z(idx_reg[0])
    else:
        controls = [idx_reg[i] for i in range(n_idx - 1)]
        target = idx_reg[n_idx - 1]
        qc.h(target)
        qc.mcx(controls, target)
        qc.h(target)

    qc.x(idx_reg)
    qc.h(idx_reg)

    return qc.to_gate(label="Diffuser")

def build_grover_circuit(payload: dict[str, Any], num_solutions: int = 1, iterations: int | None = None, measure: bool = True, ) -> tuple[QuantumCircuit, int]:
    windows_bin: list[str] = payload["windows_bin"]
    pattern_bin: str = payload["pattern_bin"]
    n_idx: int = payload["n_idx"]
    n_data: int = payload["n_data"]
    valid_len: int = payload.get("valid_len", len(windows_bin))

    k = (compute_optimal_iterations(n_idx=n_idx, num_solutions=num_solutions) if iterations is None else iterations)
    if k < 0:
        raise ValueError(f"iterations deve essere >= 0, ricevuto {k}.")

    loader_gate = build_dna_loader(windows_bin=windows_bin, n_idx=n_idx, n_data=n_data, valid_len=valid_len, )
    oracle_gate = build_pattern_oracle(loader_gate=loader_gate, pattern_bin=pattern_bin, n_idx=n_idx, n_data=n_data, )
    diffuser_gate = build_diffuser(n_idx=n_idx)

    idx_reg = QuantumRegister(n_idx, "idx")
    data_reg = QuantumRegister(n_data, "data")
    valid_reg = QuantumRegister(1, "valid")

    if measure:
        c_idx = ClassicalRegister(n_idx, "c_idx")
        qc = QuantumCircuit(idx_reg, data_reg, valid_reg, c_idx, name="Q_DNAMatch_Grover")
    else:
        qc = QuantumCircuit(idx_reg, data_reg, valid_reg, name="Q_DNAMatch_Grover")

    qc.h(idx_reg)

    idx_qubits = [idx_reg[i] for i in range(n_idx)]
    all_qubits = qc.qubits

    for _ in range(k):
        qc.append(oracle_gate, all_qubits)
        qc.append(diffuser_gate, idx_qubits)

    if measure:
        qc.measure(idx_reg, c_idx)

    return qc, k
