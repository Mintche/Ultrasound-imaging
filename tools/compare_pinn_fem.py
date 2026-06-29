#!/usr/bin/env python3
"""Compare a PINN complex field exported as CSV with a FEM reference grid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare un champ complexe PINN avec une grille FEM exportee par generate_pinn_data.x."
    )
    parser.add_argument("--fem", required=True, type=Path, help="fem_grid_*.csv")
    parser.add_argument(
        "--pinn",
        required=True,
        type=Path,
        help="CSV PINN contenant x,y,Re_U,Im_U (et facultativement f,mode)",
    )
    parser.add_argument("--frequency", type=float, help="Frequence a comparer si le CSV en contient plusieurs")
    parser.add_argument("--mode", type=int, help="Mode a comparer si le CSV en contient plusieurs")
    parser.add_argument(
        "--normalization",
        choices=("none", "max", "l2"),
        default="none",
        help="Normaliser chaque champ avant comparaison (defaut: none)",
    )
    parser.add_argument(
        "--pinn-scale",
        type=float,
        default=1.0,
        help="Facteur multiplicatif du champ PINN, par exemple U_norm (defaut: 1)",
    )
    parser.add_argument("--output", type=Path, help="Figure PNG/PDF de comparaison")
    parser.add_argument("--metrics-json", type=Path, help="Ecrire les metriques dans un JSON")
    parser.add_argument("--error-csv", type=Path, help="Ecrire le champ d'erreur point par point")
    parser.add_argument(
        "--interpolate-pinn",
        action="store_true",
        help="Autoriser une interpolation lineaire si la grille PINN differe (desactive par defaut)",
    )
    parser.add_argument("--no-plot", action="store_true", help="Calculer uniquement les metriques")
    return parser.parse_args()


def read_csv(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    if data.dtype.names is None:
        raise ValueError(f"{path}: en-tete CSV absent ou invalide")
    return np.atleast_1d(data)


def column(data: np.ndarray, aliases: Iterable[str]) -> np.ndarray:
    names = data.dtype.names or ()
    lower_to_name = {name.lower(): name for name in names}
    for alias in aliases:
        if alias.lower() in lower_to_name:
            return np.asarray(data[lower_to_name[alias.lower()]], dtype=float)
    raise ValueError(f"Colonne manquante parmi {tuple(aliases)}; colonnes disponibles: {names}")


def select_case(data: np.ndarray, frequency: float | None, mode: int | None, label: str) -> np.ndarray:
    selected = data
    names = {name.lower(): name for name in (data.dtype.names or ())}

    if "f" in names:
        values = np.unique(np.asarray(selected[names["f"]], dtype=float))
        if frequency is None:
            if len(values) != 1:
                raise ValueError(f"{label}: plusieurs frequences {values}; utilisez --frequency")
            frequency = float(values[0])
        mask = np.isclose(np.asarray(selected[names["f"]], dtype=float), frequency)
        selected = selected[mask]

    if "mode" in names:
        values = np.unique(np.asarray(selected[names["mode"]], dtype=int))
        if mode is None:
            if len(values) != 1:
                raise ValueError(f"{label}: plusieurs modes {values}; utilisez --mode")
            mode = int(values[0])
        selected = selected[np.asarray(selected[names["mode"]], dtype=int) == mode]

    if len(selected) == 0:
        raise ValueError(f"{label}: aucun point ne correspond au cas demande")
    return selected


def complex_field(data: np.ndarray) -> np.ndarray:
    real = column(data, ("Re_U", "u_re", "real", "ReU"))
    imag = column(data, ("Im_U", "u_im", "imag", "ImU"))
    return real + 1j * imag


def average_duplicate_points(x: np.ndarray, y: np.ndarray, values: np.ndarray):
    points = np.column_stack((x, y))
    unique_points, inverse = np.unique(points, axis=0, return_inverse=True)
    if len(unique_points) == len(points):
        return x, y, values
    sums = np.zeros(len(unique_points), dtype=complex)
    counts = np.zeros(len(unique_points), dtype=int)
    np.add.at(sums, inverse, values)
    np.add.at(counts, inverse, 1)
    return unique_points[:, 0], unique_points[:, 1], sums / counts


def interpolate_to_reference(
    x_source: np.ndarray,
    y_source: np.ndarray,
    values: np.ndarray,
    x_target: np.ndarray,
    y_target: np.ndarray,
    allow_interpolation: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    x_source, y_source, values = average_duplicate_points(x_source, y_source, values)
    lookup = {
        (round(float(x), 12), round(float(y), 12)): value
        for x, y, value in zip(x_source, y_source, values)
    }
    keys = [(round(float(x), 12), round(float(y), 12)) for x, y in zip(x_target, y_target)]
    if all(key in lookup for key in keys):
        return np.asarray([lookup[key] for key in keys]), np.ones(len(keys), dtype=bool)

    if not allow_interpolation:
        missing = sum(key not in lookup for key in keys)
        raise ValueError(
            f"Grilles FEM et PINN differentes ({missing} points FEM absents du CSV PINN). "
            "Evaluez le PINN sur fem_evaluation_grid_*.csv, ou utilisez --interpolate-pinn."
        )

    if len(values) < 3:
        raise ValueError("Il faut au moins trois points PINN non alignes pour interpoler le champ")

    import matplotlib.tri as mtri

    triangulation = mtri.Triangulation(x_source, y_source)
    real_interpolator = mtri.LinearTriInterpolator(triangulation, values.real)
    imag_interpolator = mtri.LinearTriInterpolator(triangulation, values.imag)
    real = real_interpolator(x_target, y_target)
    imag = imag_interpolator(x_target, y_target)
    valid = ~(np.ma.getmaskarray(real) | np.ma.getmaskarray(imag))
    interpolated = np.asarray(np.ma.filled(real, np.nan)) + 1j * np.asarray(
        np.ma.filled(imag, np.nan)
    )
    return interpolated, valid


def normalize(values: np.ndarray, method: str) -> np.ndarray:
    if method == "none":
        return values
    norm = np.max(np.abs(values)) if method == "max" else np.linalg.norm(values)
    if norm == 0.0:
        raise ValueError("Impossible de normaliser un champ nul")
    return values / norm


def compute_metrics(reference: np.ndarray, prediction: np.ndarray) -> dict[str, float | int]:
    error = prediction - reference
    reference_norm = np.linalg.norm(reference)
    return {
        "n_points": int(len(reference)),
        "relative_l2_complex": float(np.linalg.norm(error) / max(reference_norm, 1e-30)),
        "relative_l2_real": float(
            np.linalg.norm(error.real) / max(np.linalg.norm(reference.real), 1e-30)
        ),
        "relative_l2_imag": float(
            np.linalg.norm(error.imag) / max(np.linalg.norm(reference.imag), 1e-30)
        ),
        "rmse_complex": float(np.sqrt(np.mean(np.abs(error) ** 2))),
        "mae_magnitude": float(np.mean(np.abs(np.abs(prediction) - np.abs(reference)))),
        "max_abs_error": float(np.max(np.abs(error))),
    }


def plot_comparison(
    path: Path,
    x: np.ndarray,
    y: np.ndarray,
    reference: np.ndarray,
    prediction: np.ndarray,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    triangulation = mtri.Triangulation(x, y)
    error = prediction - reference
    fields = (
        (reference.real, "Re(U) FEM"),
        (prediction.real, "Re(U) PINN"),
        (error.real, "Re(erreur)"),
        (reference.imag, "Im(U) FEM"),
        (prediction.imag, "Im(U) PINN"),
        (np.abs(error), "|erreur complexe|"),
    )

    figure, axes = plt.subplots(2, 3, figsize=(14, 6), constrained_layout=True)
    for axis, (values, title) in zip(axes.flat, fields):
        image = axis.tricontourf(triangulation, values, levels=50, cmap="viridis")
        axis.set_title(title)
        axis.set_aspect("equal")
        axis.set_xlabel("x")
        axis.set_ylabel("y")
        figure.colorbar(image, ax=axis)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def main() -> None:
    args = parse_args()
    fem = select_case(read_csv(args.fem), args.frequency, args.mode, "FEM")
    pinn = select_case(read_csv(args.pinn), args.frequency, args.mode, "PINN")

    x_fem = column(fem, ("x",))
    y_fem = column(fem, ("y",))
    u_fem = complex_field(fem)
    x_pinn = column(pinn, ("x",))
    y_pinn = column(pinn, ("y",))
    u_pinn = complex_field(pinn) * args.pinn_scale

    u_pinn_on_fem, valid = interpolate_to_reference(
        x_pinn, y_pinn, u_pinn, x_fem, y_fem,
        allow_interpolation=args.interpolate_pinn,
    )
    if not np.any(valid):
        raise ValueError("Aucun point FEM ne se trouve dans le domaine couvert par les predictions PINN")

    x = x_fem[valid]
    y = y_fem[valid]
    u_fem = normalize(u_fem[valid], args.normalization)
    u_pinn_on_fem = normalize(u_pinn_on_fem[valid], args.normalization)
    metrics = compute_metrics(u_fem, u_pinn_on_fem)
    metrics["coverage"] = float(np.mean(valid))
    metrics["normalization"] = args.normalization
    metrics["pinn_scale"] = float(args.pinn_scale)

    print(json.dumps(metrics, indent=2))

    if args.metrics_json:
        args.metrics_json.parent.mkdir(parents=True, exist_ok=True)
        args.metrics_json.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")

    if args.error_csv:
        args.error_csv.parent.mkdir(parents=True, exist_ok=True)
        error = u_pinn_on_fem - u_fem
        output = np.column_stack((x, y, u_fem.real, u_fem.imag,
                                  u_pinn_on_fem.real, u_pinn_on_fem.imag,
                                  error.real, error.imag, np.abs(error)))
        np.savetxt(
            args.error_csv,
            output,
            delimiter=",",
            header="x,y,Re_U_FEM,Im_U_FEM,Re_U_PINN,Im_U_PINN,Re_error,Im_error,abs_error",
            comments="",
        )

    if not args.no_plot:
        output_path = args.output or Path("pinn_fem_comparison.png")
        plot_comparison(output_path, x, y, u_fem, u_pinn_on_fem)


if __name__ == "__main__":
    main()
