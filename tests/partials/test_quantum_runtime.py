from __future__ import annotations
from src.q_grover import build_grover_circuit
from src.q_runtime import run_circuit
import math

def test_grover_runtime_ideal_single_match_and_metrics() -> None:
    mock_payload = {
        "windows_bin": ["0001", "0110", "1011", "1100", "0010", "0101", "0000", "0000"], "pattern_bin": "1011", "n_idx": 3, "n_data": 4, "valid_len": 6,
    }

    qc, k = build_grover_circuit(mock_payload, num_solutions=1, measure=True)
    assert k == 2

    shots = 1024
    res = run_circuit(qc, shots=shots, use_noise=False, seed=42)

    assert set(res.keys()) == {"counts", "found_indices", "probabilities", "depth", "cx_count", "exec_time_ms", }
    assert sum(res["counts"].values()) == shots
    assert math.isclose(sum(res["probabilities"].values()), 1.0, abs_tol=1e-9)
    assert res["found_indices"] == [2]
    assert res["probabilities"][2] > 0.90
    assert res["depth"] > 0
    assert res["cx_count"] > 0
    assert res["exec_time_ms"] >= 0.0

def test_grover_runtime_ideal_multiple_matches() -> None:
    mock_payload = {
        "windows_bin": ["0011", "1010", "0101", "1010", "1111", "0001", "0000", "0000"], "pattern_bin": "1010", "n_idx": 3, "n_data": 4, "valid_len": 6,
    }

    qc, k = build_grover_circuit(mock_payload, num_solutions=2, measure=True)
    assert k == 1

    res = run_circuit(qc, shots=1024, use_noise=False, seed=42)

    assert res["found_indices"] == [1, 3]
    assert res["probabilities"][1] + res["probabilities"][3] > 0.99
    assert res["depth"] > 0
    assert res["cx_count"] > 0

def test_grover_runtime_noisy_fake_brisbane_metrics() -> None:
    mock_payload = {
        "windows_bin": ["00", "01", "10", "00"], "pattern_bin": "10", "n_idx": 2, "n_data": 2, "valid_len": 3,
    }

    qc, k = build_grover_circuit(mock_payload, num_solutions=1, measure=True)
    assert k == 1

    shots = 512
    res_ideal = run_circuit(qc, shots=shots, use_noise=False, seed=42)
    res_noisy = run_circuit(qc, shots=shots, use_noise=True, seed=42)

    assert res_ideal["found_indices"] == [2]
    assert sum(res_noisy["counts"].values()) == shots
    assert math.isclose(sum(res_noisy["probabilities"].values()), 1.0, abs_tol=1e-9)
    assert isinstance(res_noisy["found_indices"], list)
    assert res_noisy["depth"] > 0
    assert res_noisy["cx_count"] > 0
    assert res_noisy["exec_time_ms"] >= 0.0