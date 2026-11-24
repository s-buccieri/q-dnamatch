from __future__ import annotations
from typing import Any
from src.analytics import (build_scaling_dataframe, plot_histogram, plot_noise_vs_ideal, plot_scaling, )
from src.dna_core import process_dna
from src.q_grover import build_grover_circuit
from src.q_runtime import run_circuit
import matplotlib.pyplot as plt
import streamlit as st

st.set_page_config(page_title="Q-DNAMatch", layout="wide", initial_sidebar_state="expanded")
st.markdown("""
    <style>
        .block-container { padding-top: 2rem; padding-bottom: 2rem; max-width: 1120px; }
        .dna-box {
            font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
            font-size: 1.35rem;
            letter-spacing: 0.18rem;
            padding: 0.9rem 1.2rem;
            background-color: #F7F7F8;
            border: 1px solid #E2E4E9;
            border-radius: 6px;
            overflow-x: auto;
        }
        .dna-match {
            background-color: #D4EDDA;
            color: #155724;
            font-weight: 700;
            padding: 0.1rem 0.25rem;
            border-radius: 3px;
        }
        .dna-base { color: #2B2D33; }
        .subtle-caption { color: #6B7280; font-size: 0.9rem; margin-bottom: 1.2rem; }
    </style>
    """, unsafe_allow_html=True, )

PRESETS: dict[str, tuple[str, str]] = {
    "Singolo match (ACGTTAGC / TT)": ("ACGTTAGC", "TT"), "Protezione padding (CGAATC / AA)": ("CGAATC", "AA"), "Match multipli (ACGTACGT / CG)": ("ACGTACGT", "CG"), "Compatto NISQ (ACGT / G)": ("ACGT", "G"), "Pattern assente (ACGTAC / TT)": ("ACGTAC", "TT"),
}


def render_highlighted_dna(seq: str, pattern_len: int, indices: list[int]) -> str:
    seq_u = seq.upper()
    matched_positions: set[int] = set()
    for start in indices:
        for offset in range(pattern_len):
            if start + offset < len(seq_u):
                matched_positions.add(start + offset)

    spans = [f'<span class="dna-match">{base}</span>' if i in matched_positions else f'<span class="dna-base">{base}</span>' for i, base in enumerate(seq_u)]
    return f'<div class="dna-box">{"".join(spans)}</div>'


def _pipeline_runner_for_scaling(seq_s: str, pat_s: str) -> dict[str, Any]:
    payload_s = process_dna(seq_s, pat_s)
    s_count = max(1, len(payload_s["classical_res"]["indices"]))
    _, k_val = build_grover_circuit(payload_s, num_solutions=s_count, measure=True)
    return {
        "classical_queries": payload_s["classical_res"]["queries"], "grover_iterations": k_val,
    }


st.title("Q-DNAMatch")
st.markdown('<div class="subtle-caption">Quantum Pattern Matching su sequenze DNA tramite QROM e algoritmo di Grover — Qiskit Fall Fest 2026</div>', unsafe_allow_html=True, )

with st.sidebar:
    st.subheader("Configurazione")
    preset_label = st.selectbox("Esempio rapido", options=["Personalizzato"] + list(PRESETS.keys()), index=1, )

    default_seq, default_pat = PRESETS.get(preset_label, ("ACGTTAGC", "TT"))
    seq_input = st.text_input("Sequenza DNA (A, C, G, T)", value=default_seq).strip().upper()
    pat_input = st.text_input("Pattern da cercare", value=default_pat).strip().upper()

    st.divider()
    shots = st.select_slider("Shots simulazione", options=[256, 512, 1024, 2048, 4096], value=1024)
    seed = int(st.number_input("Seed", min_value=0, max_value=99999, value=42, step=1))

    custom_iter = st.checkbox("Forza iterazioni Grover", value=False)
    manual_k = (int(st.number_input("Numero iterazioni (k)", min_value=0, max_value=10, value=1, step=1)) if custom_iter else None)

tab_demo, tab_scaling, tab_noise = st.tabs(["Demo Live", "Superiorità Quantistica", "Hardware & Noise (NISQ)"])

payload: dict[str, Any] | None = None
grover_qc = None
grover_k = 0

try:
    payload = process_dna(seq_input, pat_input)
    classical_indices: list[int] = payload["classical_res"]["indices"]
    classical_queries: int = payload["classical_res"]["queries"]
    num_solutions_param = max(1, len(classical_indices))

    total_qubits = payload["n_idx"] + payload["n_data"] + 1
    if total_qubits > 14:
        st.warning(f"Circuito troppo grande per la simulazione live ({total_qubits} qubit). "
                   "Usa una sequenza più corta (es. fino a 16 basi) o un pattern da 1–3 basi.")
        payload = None
    else:
        grover_qc, grover_k = build_grover_circuit(payload, num_solutions=num_solutions_param, iterations=manual_k, measure=True, )
except ValueError as err:
    st.error(str(err))

with tab_demo:
    if payload is not None and grover_qc is not None:
        with st.spinner("Esecuzione circuito su AerSimulator (ideale)..."):
            ideal_res = run_circuit(grover_qc, shots=shots, use_noise=False, seed=seed)

        if not classical_indices:
            st.info("Il pattern non compare nella sequenza: l'oracolo non marca alcuno stato e la distribuzione resta uniforme.")

        st.markdown("**Sequenza analizzata (match evidenziati)**")
        st.markdown(render_highlighted_dna(seq_input, len(pat_input), ideal_res["found_indices"]), unsafe_allow_html=True, )
        st.write("")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Indici Quantum", str(ideal_res["found_indices"]))
        c2.metric("Indici Classici", str(classical_indices))
        c3.metric("Query O(N) vs Grover O(√N)", f"{classical_queries} vs {grover_k}")
        c4.metric("Qubit / Depth / 2Q-Gates", f"{grover_qc.num_qubits} / {ideal_res['depth']} / {ideal_res['cx_count']}")

        fig_hist = plot_histogram(ideal_res, target_indices=classical_indices if classical_indices else None, )
        st.pyplot(fig_hist, width="stretch")
        plt.close(fig_hist)

        with st.expander("Circuito quantistico e dettagli tecnici"):
            st.text(f"Tempo simulazione: {ideal_res['exec_time_ms']} ms | Finestre valide: {payload['valid_len']} / {1 << payload['n_idx']}")
            st.code(str(grover_qc.draw(output="text")), language="text")
            col_p, col_c = st.columns(2)
            with col_p:
                st.caption("DNAPayload")
                st.json(payload)
            with col_c:
                st.caption("ExecutionResult")
                st.json(ideal_res)

with tab_scaling:
    col_a, col_b, col_c = st.columns([1, 1, 2])
    with col_a:
        bench_pat_len = int(st.number_input("Lunghezza pattern (M)", min_value=1, max_value=6, value=2, step=1))
    with col_b:
        bench_solutions = int(st.number_input("Occorrenze attese (S)", min_value=1, max_value=4, value=1, step=1))
    with col_c:
        use_real_sim = st.checkbox("Usa la pipeline reale sui punti simulati (N = 4, 8)", value=False, )

    runner = _pipeline_runner_for_scaling if use_real_sim else None
    sim_lengths = [4, 8] if use_real_sim else [4, 8, 16]

    with st.spinner("Calcolo benchmark di scalabilità..."):
        df_scaling = build_scaling_dataframe(simulated_lengths=sim_lengths, pattern_len=bench_pat_len, num_solutions=bench_solutions, pipeline_runner=runner, seed=seed, )

    fig_scale = plot_scaling(df_scaling)
    st.pyplot(fig_scale, width="stretch")
    plt.close(fig_scale)

    st.dataframe(df_scaling, width="stretch", hide_index=True)

with tab_noise:
    if payload is not None and grover_qc is not None:
        st.caption(f"Confronto diretto tra simulatore ideale e modello di rumore hardware IBM FakeBrisbane "
                   f"({grover_qc.num_qubits} qubit, {shots} shots).")

        if grover_qc.num_qubits > 9:
            st.warning("Con più di 9 qubit la simulazione rumorosa su FakeBrisbane può richiedere diversi secondi. "
                       "Per una demo rapida seleziona il preset 'Compatto NISQ (ACGT / G)'.")

        if st.button("Esegui confronto Ideale vs FakeBrisbane", type="primary"):
            with st.spinner("Transpilazione e simulazione su FakeBrisbane in corso..."):
                res_ideal_cmp = run_circuit(grover_qc, shots=shots, use_noise=False, seed=seed)
                res_noisy_cmp = run_circuit(grover_qc, shots=shots, use_noise=True, seed=seed)

            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Depth (Ideale → Noisy)", f"{res_ideal_cmp['depth']} → {res_noisy_cmp['depth']}")
            m2.metric("2Q-Gates (Ideale → Noisy)", f"{res_ideal_cmp['cx_count']} → {res_noisy_cmp['cx_count']}")
            m3.metric("Indici Noisy", str(res_noisy_cmp["found_indices"]))
            m4.metric("Tempo Noisy", f"{res_noisy_cmp['exec_time_ms']} ms")

            fig_cmp = plot_noise_vs_ideal(res_ideal_cmp, res_noisy_cmp, target_indices=classical_indices if classical_indices else None, )
            st.pyplot(fig_cmp, width="stretch")
            plt.close(fig_cmp)
