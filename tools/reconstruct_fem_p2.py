#!/usr/bin/env python3
"""Reconstruct a complex FEM P2 field on the fixed evaluation grid."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def read_csv(path: Path) -> np.ndarray:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    if data.dtype.names is None:
        raise ValueError(f"{path}: en-tete CSV absent ou invalide")
    return np.atleast_1d(data)


def column(data: np.ndarray, name: str, dtype=float) -> np.ndarray:
    names = {candidate.lower(): candidate for candidate in (data.dtype.names or ())}
    if name.lower() not in names:
        raise ValueError(f"Colonne {name!r} absente; colonnes disponibles: {data.dtype.names}")
    return np.asarray(data[names[name.lower()]], dtype=dtype)


@dataclass(frozen=True)
class P2EvaluationGrid:
    grid_i: np.ndarray
    x: np.ndarray
    y: np.ndarray
    x_norm: np.ndarray
    y_norm: np.ndarray
    element_id: np.ndarray
    node_ids: np.ndarray
    weights: np.ndarray
    region_tag: np.ndarray
    is_defect: np.ndarray
    sound_speed: np.ndarray

    @classmethod
    def from_csv(cls, path: Path | str) -> "P2EvaluationGrid":
        data = read_csv(Path(path))
        order = np.argsort(column(data, "grid_i", int))

        def ordered(name: str, dtype=float) -> np.ndarray:
            return column(data, name, dtype)[order]

        node_ids = np.column_stack([ordered(f"node{i}", int) for i in range(6)])
        weights = np.column_stack([ordered(f"phi{i}") for i in range(6)])
        if not np.allclose(np.sum(weights, axis=1), 1.0, rtol=0.0, atol=5e-13):
            raise ValueError("Les poids P2 ne verifient pas la partition de l'unite")

        return cls(
            grid_i=ordered("grid_i", int),
            x=ordered("x"),
            y=ordered("y"),
            x_norm=ordered("x_norm"),
            y_norm=ordered("y_norm"),
            element_id=ordered("element_id", int),
            node_ids=node_ids,
            weights=weights,
            region_tag=ordered("region_tag", int),
            is_defect=ordered("is_defect", int),
            sound_speed=ordered("c"),
        )

    def evaluate(self, nodal_values: np.ndarray) -> np.ndarray:
        """Evaluate U_h(x)=sum_i phi_i(x) U_i at every fixed-grid point."""
        nodal_values = np.asarray(nodal_values, dtype=complex)
        if nodal_values.ndim != 1:
            raise ValueError("nodal_values doit etre un vecteur complexe 1D")
        if np.min(self.node_ids) < 0 or np.max(self.node_ids) >= len(nodal_values):
            raise ValueError("La connectivite de la grille reference un node_id absent")
        return np.sum(self.weights * nodal_values[self.node_ids], axis=1)


@dataclass(frozen=True)
class NodalP2Field:
    frequency: float
    k0: float
    mode: int
    values: np.ndarray

    @classmethod
    def from_csv(
        cls,
        path: Path | str,
        frequency: float | None = None,
        mode: int | None = None,
    ) -> "NodalP2Field":
        data = read_csv(Path(path))
        frequencies = np.unique(column(data, "f"))
        modes = np.unique(column(data, "mode", int))

        if frequency is None:
            if len(frequencies) != 1:
                raise ValueError(f"Plusieurs frequences {frequencies}; precisez --frequency")
            frequency = float(frequencies[0])
        if mode is None:
            if len(modes) != 1:
                raise ValueError(f"Plusieurs modes {modes}; precisez --mode")
            mode = int(modes[0])

        selected = data[
            np.isclose(column(data, "f"), frequency)
            & (column(data, "mode", int) == mode)
        ]
        if len(selected) == 0:
            raise ValueError("Aucun coefficient nodal ne correspond au cas demande")

        node_ids = column(selected, "node_id", int)
        if len(np.unique(node_ids)) != len(node_ids):
            raise ValueError("node_id duplique dans le champ nodal")
        values = np.full(int(np.max(node_ids)) + 1, np.nan + 1j * np.nan, dtype=complex)
        values[node_ids] = column(selected, "Re_U") + 1j * column(selected, "Im_U")
        if np.any(~np.isfinite(values)):
            raise ValueError("La numerotation nodale n'est pas complete ou U contient NaN/Inf")

        k0_values = np.unique(column(selected, "k0"))
        if len(k0_values) != 1:
            raise ValueError("Plusieurs valeurs de k0 pour un meme cas")
        return cls(float(frequency), float(k0_values[0]), int(mode), values)


def write_reconstructed_csv(
    path: Path,
    grid: P2EvaluationGrid,
    field: NodalP2Field,
    values: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "f", "k0", "mode", "grid_i", "x", "y", "x_norm", "y_norm",
        "element_id", "region_tag", "is_defect", "c", "Re_U", "Im_U", "abs_U",
    )
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(header)
        for index, value in enumerate(values):
            writer.writerow((
                f"{field.frequency:.17g}", f"{field.k0:.17g}", field.mode,
                int(grid.grid_i[index]), f"{grid.x[index]:.17g}", f"{grid.y[index]:.17g}",
                f"{grid.x_norm[index]:.17g}", f"{grid.y_norm[index]:.17g}",
                int(grid.element_id[index]), int(grid.region_tag[index]),
                int(grid.is_defect[index]), f"{grid.sound_speed[index]:.17g}",
                f"{value.real:.17g}", f"{value.imag:.17g}", f"{abs(value):.17g}",
            ))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reconstruit exactement l'interpolant FEM P2 sur la grille fixe exportee."
    )
    parser.add_argument("--grid-map", required=True, type=Path, help="fem_evaluation_grid_*.csv")
    parser.add_argument("--field", required=True, type=Path, help="fem_field_*_modeN.csv")
    parser.add_argument("--frequency", type=float)
    parser.add_argument("--mode", type=int)
    parser.add_argument("--output", required=True, type=Path, help="CSV reconstruit au format fem_grid_*")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grid = P2EvaluationGrid.from_csv(args.grid_map)
    field = NodalP2Field.from_csv(args.field, args.frequency, args.mode)
    values = grid.evaluate(field.values)
    write_reconstructed_csv(args.output, grid, field, values)
    print(
        f"{len(values)} valeurs P2 reconstruites pour f={field.frequency:g} Hz, "
        f"mode={field.mode}, dans {args.output}"
    )


if __name__ == "__main__":
    main()
