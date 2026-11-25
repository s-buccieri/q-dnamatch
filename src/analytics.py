import math
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

def compute_grover_iterations(num_windows: int, num_solutions: int = 1) -> int:
    if num_windows <= 0 or num_solutions <= 0:
        return 0
    n_idx = max(1, math.ceil(math.log2(num_windows)))
    space_size = 2 ** n_idx
    if num_solutions >= space_size:
        return 0
    return max(1, math.floor((math.pi / 4.0) * math.sqrt(space_size / num_solutions)))


def build_scaling_dataframe(simulated_lengths: list[int] | None = None, analytical_lengths: list[int] | None = None, pattern_len: int = 2, num_solutions: int = 1, pipeline_runner: Callable[[str, str], dict[str, Any]] | None = None, seed: int = 42, ) -> pd.DataFrame:
    if simulated_lengths is None:
        simulated_lengths = [4, 8, 16]
    if analytical_lengths is None:
        analytical_lengths = [32, 64, 128, 256, 512, 1024, 4096, 16384, 65536, 262144, 1_000_000]

    rng = np.random.default_rng(seed)
    bases = np.array(["A", "C", "G", "T"])
    rows: list[dict[str, Any]] = []

    for n in simulated_lengths:
        valid_windows = max(1, n - pattern_len + 1)
        if pipeline_runner is not None:
            seq = "".join(rng.choice(bases, size=n))
            pattern = seq[:pattern_len]
            res = pipeline_runner(seq, pattern)
            classical_q = int(res["classical_queries"])
            grover_k = int(res["grover_iterations"])
        else:
            classical_q = valid_windows
            grover_k = compute_grover_iterations(valid_windows, num_solutions=num_solutions)

        speedup = classical_q / grover_k if grover_k > 0 else 1.0
        rows.append({
            "N": n, "valid_windows": valid_windows, "classical_queries": classical_q, "grover_iterations": grover_k, "speedup": round(speedup, 2), "mode": "simulated",
        })

    for n in analytical_lengths:
        valid_windows = max(1, n - pattern_len + 1)
        classical_q = valid_windows
        grover_k = compute_grover_iterations(valid_windows, num_solutions=num_solutions)
        speedup = classical_q / grover_k if grover_k > 0 else 1.0
        rows.append({
            "N": n, "valid_windows": valid_windows, "classical_queries": classical_q, "grover_iterations": grover_k, "speedup": round(speedup, 2), "mode": "analytical",
        })

    df = pd.DataFrame(rows).sort_values("N").reset_index(drop=True)
    return df


def plot_scaling(df: pd.DataFrame) -> Figure:
    fig, ax = plt.subplots(figsize=(8, 5), dpi=120)

    ax.loglog(df["N"], df["classical_queries"], color="#D62728", linewidth=2.2, marker="o", markersize=5, label="Ricerca Classica — $O(N)$", )
    ax.loglog(df["N"], df["grover_iterations"], color="#1F77B4", linewidth=2.2, marker="s", markersize=5, label="Grover Quantistico — $O(\\sqrt{N})$", )

    sim_df = df[df["mode"] == "simulated"]
    if not sim_df.empty:
        ax.axvspan(sim_df["N"].min() * 0.8, sim_df["N"].max() * 1.25, color="#2CA02C", alpha=0.10, label="Zona simulata (Qiskit Aer)", )

    last_row = df.iloc[-1]
    ax.annotate(f"Speedup ~{int(last_row['speedup'])}x\n(N = {int(last_row['N']):,})", xy=(last_row["N"], last_row["grover_iterations"]), xytext=(last_row["N"] / 25, last_row["grover_iterations"] * 12), arrowprops=dict(arrowstyle="->", color="#1F77B4", lw=1.2), fontsize=9, fontweight="bold", color="#1F77B4", bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#1F77B4", alpha=0.9), )

    ax.set_title("Scalabilità: Ricerca Classica $O(N)$ vs Grover $O(\\sqrt{N})$", fontsize=12, fontweight="bold")
    ax.set_xlabel("Lunghezza sequenza DNA ($N$ basi) — scala log", fontsize=10)
    ax.set_ylabel("Numero di query / iterazioni — scala log", fontsize=10)
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", frameon=True)
    fig.tight_layout()
    return fig


def plot_histogram(exec_result: dict[str, Any], target_indices: list[int] | None = None, title: str = "Distribuzione delle Probabilità (Misura Indici DNA)", ) -> Figure:
    probabilities: dict[int, float] = exec_result.get("probabilities", {})
    found_indices: list[int] = exec_result.get("found_indices", [])
    highlight = set(target_indices) if target_indices is not None else set(found_indices)

    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)

    if not probabilities:
        ax.text(0.5, 0.5, "Nessun dato disponibile", ha="center", va="center", fontsize=12)
        ax.set_title(title, fontsize=12, fontweight="bold")
        fig.tight_layout()
        return fig

    indices = sorted(probabilities.keys())
    probs = [probabilities[i] for i in indices]
    colors = ["#2CA02C" if i in highlight else "#7F7F7F" for i in indices]

    bars = ax.bar([str(i) for i in indices], probs, color=colors, edgecolor="black", alpha=0.88)

    for bar, prob in zip(bars, probs):
        if prob >= 0.02:
            ax.text(bar.get_x() + bar.get_width() / 2.0, prob + 0.015, f"{prob * 100:.1f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", )

    ax.set_ylim(0, min(1.12, max(probs) + 0.15))
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Indice finestra DNA", fontsize=10)
    ax.set_ylabel("Probabilità", fontsize=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    from matplotlib.patches import Patch

    legend_handles = [
        Patch(facecolor="#2CA02C", edgecolor="black", label="Match trovato (sopra soglia)"), Patch(facecolor="#7F7F7F", edgecolor="black", label="Rumore di fondo"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=True)
    fig.tight_layout()
    return fig


def plot_noise_vs_ideal(ideal_result: dict[str, Any], noisy_result: dict[str, Any], target_indices: list[int] | None = None, ) -> Figure:
    ideal_probs: dict[int, float] = ideal_result.get("probabilities", {})
    noisy_probs: dict[int, float] = noisy_result.get("probabilities", {})
    all_indices = sorted(set(ideal_probs.keys()) | set(noisy_probs.keys()))

    fig, ax = plt.subplots(figsize=(9, 4.8), dpi=120)

    if not all_indices:
        ax.text(0.5, 0.5, "Nessun dato disponibile", ha="center", va="center", fontsize=12)
        ax.set_title("Confronto Ideale vs Rumoroso (NISQ)", fontsize=12, fontweight="bold")
        fig.tight_layout()
        return fig

    x = np.arange(len(all_indices))
    width = 0.38

    y_ideal = [ideal_probs.get(i, 0.0) for i in all_indices]
    y_noisy = [noisy_probs.get(i, 0.0) for i in all_indices]

    ax.bar(x - width / 2, y_ideal, width, label="Ideale (AerSimulator)", color="#1F77B4", edgecolor="black")
    ax.bar(x + width / 2, y_noisy, width, label="Rumoroso (FakeBrisbane)", color="#FF7F0E", edgecolor="black")

    info_lines = [
        f"Ideale   — depth: {ideal_result.get('depth', '-')}, 2Q-gates: {ideal_result.get('cx_count', '-')}", f"Rumoroso — depth: {noisy_result.get('depth', '-')}, 2Q-gates: {noisy_result.get('cx_count', '-')}",
    ]
    if target_indices:
        ideal_fid = sum(ideal_probs.get(i, 0.0) for i in target_indices)
        noisy_fid = sum(noisy_probs.get(i, 0.0) for i in target_indices)
        info_lines.append(f"Prob. target — Ideale: {ideal_fid * 100:.1f}% | Rumoroso: {noisy_fid * 100:.1f}%")

    ax.text(0.98, 0.95, "\n".join(info_lines), transform=ax.transAxes, ha="right", va="top", fontsize=8.5, family="monospace", bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#444444", alpha=0.92), )

    ax.set_xticks(x)
    ax.set_xticklabels([str(i) for i in all_indices])
    ax.set_ylim(0, min(1.15, max(max(y_ideal), max(y_noisy)) + 0.22))
    ax.set_title("Analisi NISQ: Simulatore Ideale vs Rumoroso (FakeBrisbane)", fontsize=12, fontweight="bold")
    ax.set_xlabel("Indice finestra DNA", fontsize=10)
    ax.set_ylabel("Probabilità", fontsize=10)
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.legend(loc="upper left", frameon=True)
    fig.tight_layout()
    return fig


def export_figure(fig: Figure, output_path: str | Path, dpi: int = 300) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    return path
