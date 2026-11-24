from __future__ import annotations
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit import Gate

def input_validation(windows_bin: list[str], n_idx: int, n_data: int, valid_len: int | None = None) -> None:
    if n_idx < 1:
        raise ValueError(f"n_idx deve essere >= 1, ricevuto {n_idx}.")
    if n_data < 1:
        raise ValueError(f"n_data deve essere >= 1, ricevuto {n_data}.")

    max_states = 1 << n_idx
    effective_valid_len = len(windows_bin) if valid_len is None else valid_len

    if effective_valid_len < 0 or effective_valid_len > max_states:
        raise ValueError(f"valid_len ({effective_valid_len}) deve essere compreso tra 0 e 2**n_idx ({max_states}).")
    if effective_valid_len > len(windows_bin):
        raise ValueError(f"valid_len ({effective_valid_len}) supera len(windows_bin) ({len(windows_bin)}).")

    for i in range(effective_valid_len):
        w = windows_bin[i]
        if len(w) != n_data or any(ch not in ("0", "1") for ch in w):
            raise ValueError(f"Finestra all'indice {i} non valida: '{w}'. " f"Attesa stringa binaria di lunghezza {n_data}.")

def build_qrom_dna_gate(windows_bin: list[str], n_idx: int, n_data: int, valid_len: int | None = None) -> Gate:
    idx_reg = QuantumRegister(n_idx, "idx")
    data_reg = QuantumRegister(n_data, "data")
    valid_reg = QuantumRegister(1, "valid")
    qc = QuantumCircuit(idx_reg, data_reg, valid_reg, name="QROM_DNA")

    ctrl_qubits = list(idx_reg[i] for i in range(n_idx))
    effective_valid_len = len(windows_bin) if valid_len is None else valid_len

    for i in range(effective_valid_len):
        window = windows_bin[i]
        idx_bits_lsb = format(i, f"0{n_idx}b")[::-1]
        zero_controls = [idx_reg[q] for q, bit in enumerate(idx_bits_lsb) if bit == "0"]

        for q in zero_controls:
            qc.x(q)

        for bit_pos, bit_val in enumerate(window):
            if bit_val == "1":
                qc.mcx(ctrl_qubits, data_reg[bit_pos])

        qc.mcx(ctrl_qubits, valid_reg[0])

        for q in zero_controls:
            qc.x(q)

    return qc.to_gate(label="QROM_DNA")

def build_dna_loader(windows_bin: list[str], n_idx: int, n_data: int, valid_len: int | None = None, ) -> Gate:
    input_validation(windows_bin, n_idx, n_data, valid_len)
    return build_qrom_dna_gate(windows_bin, n_idx, n_data, valid_len)