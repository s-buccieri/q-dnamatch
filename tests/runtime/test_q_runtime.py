import pytest
from qiskit import QuantumCircuit

from src.q_runtime import (
    count_two_qubit_gates,
    decode_results,
    run_circuit,
)


def _build_mock_target_circuit(n_idx: int, target_index: int) -> QuantumCircuit:
    """
    Crea un circuito mock che prepara in modo deterministico lo stato |target_index>
    su n_idx qubit (rispettando il little-endian di Qiskit: q0 = bit meno significativo)
    e misura tutti i qubit di indice.
    """
    qc = QuantumCircuit(n_idx, n_idx)
    # Converte target_index in binario a n_idx bit (es. 1 su 3 bit -> '001')
    bitstring = format(target_index, f"0{n_idx}b")
    # In Qiskit q0 è il bit più a destra: invertiamo la stringa per iterare su q0, q1, ...
    for qubit_idx, bit in enumerate(reversed(bitstring)):
        if bit == "1":
            qc.x(qubit_idx)

    qc.measure(range(n_idx), range(n_idx))
    return qc


class TestDecodeResults:
    """Suite JUnit-style per la decodifica little-endian e il filtro del rumore."""

    def test_endianness_little_endian(self):
        # '001' deve corrispondere all'indice intero 1, '100' a 4
        counts = {"001": 900, "100": 100}
        found_indices, probabilities = decode_results(counts, threshold=0.5)

        assert found_indices == [1]
        assert probabilities[1] == pytest.approx(0.9)
        assert probabilities[4] == pytest.approx(0.1)

    def test_noise_filtering_multiple_peaks(self):
        # Due picchi reali (indici 2 e 5) + rumore di fondo (indici 0 e 7)
        counts = {"010": 450, "101": 450, "000": 50, "111": 50}
        found_indices, probabilities = decode_results(counts, threshold=0.2)

        assert found_indices == [2, 5]
        assert sum(probabilities.values()) == pytest.approx(1.0)

    def test_empty_counts_edge_case(self):
        found_indices, probabilities = decode_results({})
        assert found_indices == []
        assert probabilities == {}

    def test_default_dynamic_threshold(self):
        # Caso 1: 3 bit -> 8 stati -> soglia dinamica = max(0.10, 2/8) = 0.25
        counts_3bit = {"001": 600, "010": 200, "100": 200}
        found_3bit, _ = decode_results(counts_3bit, threshold=None)
        assert found_3bit == [1]

        # Caso 2: 5 bit -> 32 stati -> 2/32 = 0.0625 -> scatta il minimo 0.10
        counts_5bit = {"00000": 770, "00001": 150, "00010": 80}
        found_5bit, _ = decode_results(counts_5bit, threshold=None)
        assert found_5bit == [0, 1]

    def test_bitstrings_with_spaces(self):
        # Qiskit separa con spazi registri classici multipli: "01 0" -> "010" (2)
        counts = {"01 0": 800, "10 1": 200}
        found_indices, probabilities = decode_results(counts, threshold=0.5)

        assert found_indices == [2]
        assert probabilities[2] == pytest.approx(0.8)
        assert probabilities[5] == pytest.approx(0.2)


class TestQRuntimeExecution:
    """Suite JUnit-style per il motore di esecuzione (ideale, rumoroso e scalabilità)."""

    def _assert_contract(self, res: dict):
        """Verifica che il dizionario rispetti rigorosamente il contratto ExecutionResult."""
        expected_keys = {
            "counts",
            "found_indices",
            "probabilities",
            "depth",
            "cx_count",
            "exec_time_ms",
        }
        assert set(res.keys()) == expected_keys
        assert isinstance(res["counts"], dict)
        assert isinstance(res["found_indices"], list)
        assert isinstance(res["probabilities"], dict)
        assert isinstance(res["depth"], int) and res["depth"] >= 0
        assert isinstance(res["cx_count"], int) and res["cx_count"] >= 0
        assert isinstance(res["exec_time_ms"], float) and res["exec_time_ms"] >= 0.0

    def test_execution_ideal_contract_and_endianness(self):
        # Circuito mock 3 qubit con target = 1 ('001')
        qc = _build_mock_target_circuit(n_idx=3, target_index=1)
        res = run_circuit(qc, shots=512, use_noise=False, seed=42)

        self._assert_contract(res)
        assert res["found_indices"] == [1]
        assert res["probabilities"][1] == pytest.approx(1.0)

    def test_execution_noisy_fake_brisbane(self):
        # Circuito mock con un CNOT per verificare anche il conteggio dei gate a 2 qubit
        # |00> -> X su q0 -> |01> -> CX(0, 1) -> |11> (indice 3)
        qc = QuantumCircuit(2, 2)
        qc.x(0)
        qc.cx(0, 1)
        qc.measure([0, 1], [0, 1])

        res = run_circuit(qc, shots=512, use_noise=True, seed=42)

        self._assert_contract(res)
        # Su FakeBrisbane il gate a 2 qubit (cx -> ecr) deve essere contato (>= 1)
        assert res["cx_count"] >= 1
        # Anche con il rumore, su un circuito così corto l'indice 3 deve dominare
        assert 3 in res["found_indices"]
        assert res["probabilities"][3] > 0.7

    @pytest.mark.parametrize("n_idx", [2, 3, 4])
    def test_scaling_execution_times_n_idx(self, n_idx):
        # Verifica esecuzione e tempi per n_idx in {2, 3, 4} come da specifica
        target_index = 1
        qc = _build_mock_target_circuit(n_idx=n_idx, target_index=target_index)
        res = run_circuit(qc, shots=512, use_noise=False, seed=42)

        self._assert_contract(res)
        assert res["found_indices"] == [target_index]
        assert res["exec_time_ms"] > 0.0

    def test_count_two_qubit_gates_helper(self):
        qc = QuantumCircuit(3)
        qc.h(0)
        qc.cx(0, 1)
        qc.cx(1, 2)
        qc.barrier()
        assert count_two_qubit_gates(qc) == 2

    def test_transpiler_optimization_level_3_cancels_redundant_cx(self):
        # Due CX consecutivi identici si annullano: opt level 3 deve rimuoverli
        qc = QuantumCircuit(2, 2)
        qc.x(0)
        qc.cx(0, 1)
        qc.cx(0, 1)
        qc.measure([0, 1], [0, 1])

        res = run_circuit(qc, shots=512, use_noise=False, seed=42)

        self._assert_contract(res)
        assert res["cx_count"] == 0
        assert res["found_indices"] == [1]

    def test_circuit_with_unmeasured_data_and_valid_qubits(self):
        # Simula q_grover.py: 5 qubit totali (2 idx + 2 data + 1 valid),
        # ma misura solo i 2 qubit di indice (q0, q1). Target = 2 ('10' -> X su q1).
        qc = QuantumCircuit(5, 2)
        qc.x(1)  # idx = '10' (2)
        qc.x(2)  # data non misurato
        qc.x(4)  # valid non misurato
        qc.measure([0, 1], [0, 1])

        res = run_circuit(qc, shots=512, use_noise=False, seed=42)

        self._assert_contract(res)
        assert res["found_indices"] == [2]
        assert res["probabilities"][2] == pytest.approx(1.0)
        assert all(len(bits) == 2 for bits in res["counts"].keys())