import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

# Backend non interattivo: evita che pytest apra finestre grafiche su Windows
matplotlib.use("Agg")

from src.analytics import (build_scaling_dataframe, compute_grover_iterations, export_figure, plot_histogram, plot_noise_vs_ideal, plot_scaling, )


class TestScalingBenchmark:
    """Suite JUnit-style per formula di Grover, DataFrame di benchmark e casi limite."""

    def test_compute_grover_iterations_exact_and_padding(self):
        # Potenze di 2 esatte
        assert compute_grover_iterations(num_windows=8, num_solutions=1) == 2
        assert compute_grover_iterations(num_windows=8, num_solutions=2) == 1
        # Non potenza di 2: 5 finestre -> padding a 2^3 = 8 stati -> k = 2
        assert compute_grover_iterations(num_windows=5, num_solutions=1) == 2
        # Caso minimo valido: 1 finestra -> n_idx=1 -> space_size=2, S=1 -> k = 1
        assert compute_grover_iterations(num_windows=1, num_solutions=1) == 1

    def test_compute_grover_iterations_edge_cases(self):
        # Input nulli o negativi
        assert compute_grover_iterations(num_windows=0, num_solutions=1) == 0
        assert compute_grover_iterations(num_windows=-4, num_solutions=1) == 0
        assert compute_grover_iterations(num_windows=8, num_solutions=0) == 0
        assert compute_grover_iterations(num_windows=8, num_solutions=-2) == 0
        # Soluzioni >= spazio totale: nessuna iterazione necessaria
        assert compute_grover_iterations(num_windows=4, num_solutions=4) == 0
        assert compute_grover_iterations(num_windows=4, num_solutions=10) == 0

    def test_scaling_dataframe_reproducibility_and_seed_effect(self):
        df1 = build_scaling_dataframe(seed=42)
        df2 = build_scaling_dataframe(seed=42)
        pd.testing.assert_frame_equal(df1, df2)

        # Verifica che il seed controlli davvero le sequenze DNA passate al runner
        seqs_seed_42, seqs_seed_99 = [], []

        def capture(store):
            def _runner(seq: str, pattern: str) -> dict:
                store.append((seq, pattern))
                return {"classical_queries": len(seq) - len(pattern) + 1, "grover_iterations": 1}

            return _runner

        build_scaling_dataframe(simulated_lengths=[16], analytical_lengths=[], pipeline_runner=capture(seqs_seed_42), seed=42)
        build_scaling_dataframe(simulated_lengths=[16], analytical_lengths=[], pipeline_runner=capture(seqs_seed_99), seed=99)

        assert seqs_seed_42 != seqs_seed_99
        # Verifica alfabeto DNA valido
        assert set(seqs_seed_42[0][0]).issubset({"A", "C", "G", "T"})

    def test_scaling_dataframe_math_and_no_scale_errors(self):
        df = build_scaling_dataframe(seed=42)
        expected_cols = {"N", "valid_windows", "classical_queries", "grover_iterations", "speedup", "mode", }
        assert set(df.columns) == expected_cols
        assert df["N"].is_monotonic_increasing
        assert df["N"].iloc[-1] == 1_000_000

        # Controllo matematico riga per riga dello speedup
        for _, row in df.iterrows():
            expected_speedup = round(row["classical_queries"] / row["grover_iterations"], 2)
            assert row["speedup"] == pytest.approx(expected_speedup)
            assert row["grover_iterations"] >= 1

        # A N = 10^6 il vantaggio quadratico deve superare 1000x
        last = df.iloc[-1]
        assert last["classical_queries"] == 999_999
        assert last["speedup"] > 1000

    def test_scaling_with_mock_pipeline_runner_and_pattern_longer_than_n(self):
        called = []

        def mock_pipeline(seq: str, pattern: str) -> dict:
            called.append((len(seq), len(pattern)))
            return {"classical_queries": 7, "grover_iterations": 2}

        df = build_scaling_dataframe(simulated_lengths=[4, 8], analytical_lengths=[1024], pattern_len=2, pipeline_runner=mock_pipeline, seed=42, )
        assert called == [(4, 2), (8, 2)]
        sim_rows = df[df["mode"] == "simulated"]
        assert list(sim_rows["classical_queries"]) == [7, 7]
        assert list(sim_rows["grover_iterations"]) == [2, 2]
        assert list(sim_rows["speedup"]) == [3.5, 3.5]


class TestAnalyticsPlots:
    """Suite JUnit-style che verifica contenuto reale, scale, colori e testi dei 3 grafici."""

    def test_plot_scaling_axes_data_and_no_simulated_mode(self):
        df = build_scaling_dataframe(seed=42)
        fig = plot_scaling(df)
        try:
            assert isinstance(fig, Figure)
            ax = fig.axes[0]
            assert ax.get_xscale() == "log"
            assert ax.get_yscale() == "log"
            assert len(ax.lines) == 2

            # Verifica che i punti plottati corrispondano ai dati del DataFrame
            np.testing.assert_array_equal(ax.lines[0].get_xdata(), df["N"].to_numpy())
            np.testing.assert_array_equal(ax.lines[0].get_ydata(), df["classical_queries"].to_numpy())
            np.testing.assert_array_equal(ax.lines[1].get_ydata(), df["grover_iterations"].to_numpy())
        finally:
            plt.close(fig)

        # Caso limite: DataFrame senza punti simulati (solo analytical)
        df_only_analytical = df[df["mode"] == "analytical"].reset_index(drop=True)
        fig_an = plot_scaling(df_only_analytical)
        try:
            assert isinstance(fig_an, Figure)
            assert len(fig_an.axes[0].lines) == 2
        finally:
            plt.close(fig_an)

    def test_plot_histogram_bars_heights_and_colors(self):
        mock_exec_res = {
            "found_indices": [2], "probabilities": {0: 0.05, 1: 0.05, 2: 0.85, 3: 0.05},
        }
        # Test senza target_indices esplicito: deve usare found_indices ([2])
        fig = plot_histogram(mock_exec_res)
        try:
            ax = fig.axes[0]
            patches = ax.patches
            assert len(patches) == 4
            heights = [p.get_height() for p in patches]
            assert heights == pytest.approx([0.05, 0.05, 0.85, 0.05])

            # Barra dell'indice 2 deve essere verde (#2CA02C), le altre grigie (#7F7F7F)
            green_rgb = mcolors.to_rgb("#2CA02C")
            gray_rgb = mcolors.to_rgb("#7F7F7F")
            assert patches[2].get_facecolor()[:3] == pytest.approx(green_rgb)
            assert patches[0].get_facecolor()[:3] == pytest.approx(gray_rgb)
        finally:
            plt.close(fig)

    def test_plot_histogram_empty_input(self):
        fig_empty = plot_histogram({})
        try:
            ax = fig_empty.axes[0]
            assert len(ax.patches) == 0
            assert any("Nessun dato" in t.get_text() for t in ax.texts)
        finally:
            plt.close(fig_empty)

    def test_plot_noise_vs_ideal_union_indices_and_metrics_box(self):
        # Caso reale NISQ: il rumoroso ha indici spuri (1 e 3) assenti nell'ideale
        mock_ideal = {
            "probabilities": {2: 0.98, 0: 0.02}, "depth": 24, "cx_count": 12,
        }
        mock_noisy = {
            "probabilities": {0: 0.10, 1: 0.15, 2: 0.60, 3: 0.15}, "depth": 68, "cx_count": 34,
        }
        fig = plot_noise_vs_ideal(mock_ideal, mock_noisy, target_indices=[2])
        try:
            ax = fig.axes[0]
            # Unione degli indici: {0, 1, 2, 3} -> 4 coppie di barre (2 container da 4 barre)
            assert len(ax.containers) == 2
            ideal_heights = [b.get_height() for b in ax.containers[0]]
            noisy_heights = [b.get_height() for b in ax.containers[1]]
            assert ideal_heights == pytest.approx([0.02, 0.0, 0.98, 0.0])
            assert noisy_heights == pytest.approx([0.10, 0.15, 0.60, 0.15])

            # Verifica che il box informativo contenga depth, 2Q-gates e fedeltà target
            all_text = "\n".join(t.get_text() for t in ax.texts)
            assert "depth: 24" in all_text and "2Q-gates: 12" in all_text
            assert "depth: 68" in all_text and "2Q-gates: 34" in all_text
            assert "98.0%" in all_text and "60.0%" in all_text
        finally:
            plt.close(fig)

    def test_plot_noise_vs_ideal_empty_input(self):
        fig_empty = plot_noise_vs_ideal({}, {})
        try:
            ax = fig_empty.axes[0]
            assert len(ax.containers) == 0
            assert any("Nessun dato" in t.get_text() for t in ax.texts)
        finally:
            plt.close(fig_empty)

    def test_export_figure_png_svg_and_nested_dirs(self, tmp_path):
        df = build_scaling_dataframe(seed=42)
        fig = plot_scaling(df)
        png_path = tmp_path / "slides" / "nested" / "fig_scaling.png"
        svg_path = tmp_path / "slides" / "fig_scaling.svg"
        try:
            saved_png = export_figure(fig, png_path, dpi=150)
            saved_svg = export_figure(fig, svg_path)
            assert saved_png.exists() and saved_png.stat().st_size > 0
            assert saved_svg.exists() and saved_svg.stat().st_size > 0
        finally:
            plt.close(fig)
