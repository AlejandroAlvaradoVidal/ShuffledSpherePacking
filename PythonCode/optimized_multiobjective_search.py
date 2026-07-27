from __future__ import annotations

import datetime
import math
import os
import random
from dataclasses import dataclass
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import torch


# ============================================================
# Configuration
# ============================================================

@dataclass
class Config:
    NF: int = 31
    R: int = 1
    mux_values: Tuple[int, ...] = tuple(range(1, 17))

    max_iterations: int = 4000
    patience: int = 700
    burn_in: int = 50
    print_every: int = 50

    # Multi-objective weights used by the simulated-annealing energy.
    # Density is maximized; condition number and uniformity error are minimized.
    weight_density: float = 0.55
    weight_condition: float = 0.25
    weight_uniformity: float = 0.20

    # Constraint penalty. Large values make the search strongly prefer feasible points.
    constraint_penalty: float = 20.0

    # Maximum admissible condition number relative to the initial condition number.
    # Example: 1.05 means at most 5% above the initial value.
    condition_limit_factor: float = 1.05

    # Progressive normalized-uniformity constraint.
    # The limit decreases linearly from max(initial error, start) to final.
    uniformity_limit_start: float = 0.20
    uniformity_limit_final: float = 0.01

    # Simulated annealing prevents a completely greedy search.
    initial_temperature: float = 0.03
    final_temperature: float = 1e-4

    # Candidate generation.
    critical_move_probability: float = 0.85
    extra_random_moves_probability: float = 0.10
    max_extra_random_moves: int = 2

    # Condition-number definition:
    # False -> cond(sigma), numerically safer.
    # True  -> cond(sigma.T @ sigma), compatible with the original code.
    use_gram_condition: bool = True

    # Optional post-processing. Keeping this enabled makes the output comparable
    # with the original workflow, although uniformity is already used during search.
    apply_final_uniformization: bool = True
    uniformizer_max_iter_factor: int = 200

    random_seed: int = 1234
    plot_graphs: bool = True
    graph_export_path: str = "./PythonCode/graphs_multiobjective"
    output_folder: str = "./PythonCode/Patrones/SSPv5_multiobjective"
    log_filename: str = "resumen_metricas_multiobjetivo.txt"


# ============================================================
# User-specific initialization functions
# ============================================================

def run_eq10mux5(NF: int, M: int, device: torch.device) -> torch.Tensor:
    """
    Placeholder for the user's original RunEQ10MUX5 implementation.

    Replace this function with the real generator when available.
    The expected output shape is (NF, NF, M), with integer values in [1, NF].
    """
    return torch.randint(1, NF + 1, (NF, NF, M), device=device)


def generate_coded_aperture(
    G: torch.Tensor,
    NF: int,
    B: int,
    N: int,
    device: torch.device,
) -> Tuple[torch.Tensor, torch.Tensor]:
    C = torch.zeros((NF, NF, NF), device=device)
    T = torch.zeros((N, N, NF), device=device)

    for j in range(1, NF + 1):
        for i in range(B):
            C[:, :, j - 1] += (G[:, :, i] == j).float()

    C.clamp_(0, 1)
    return C, T


# ============================================================
# Metrics
# ============================================================

def get_cond(C: torch.Tensor, use_gram: bool = True) -> float:
    """Compute cond(sigma) or cond(sigma.T @ sigma)."""
    NF = C.shape[2]
    sigma = C.reshape(-1, NF).float()

    max_val = sigma.abs().max()
    if max_val > 0:
        sigma = sigma / max_val

    matrix = sigma.T @ sigma if use_gram else sigma

    try:
        singular_values = torch.linalg.svdvals(matrix)
        if singular_values.numel() == 0:
            return float("inf")

        s_max = singular_values.max()
        s_min = singular_values.min()
        eps = torch.finfo(singular_values.dtype).eps

        if not torch.isfinite(s_max) or not torch.isfinite(s_min):
            return float("inf")
        if s_min <= eps * max(float(max(matrix.shape)), 1.0) * s_max:
            return float("inf")

        return (s_max / s_min).item()
    except RuntimeError:
        return float("inf")


def compute_density_irregular_sp(T: torch.Tensor) -> Tuple[float, float]:
    """Original irregular sphere-packing density metric, computed on GPU."""
    NF = T.shape[2]
    coords = torch.nonzero(T > 0)

    if coords.shape[0] < 2:
        return 0.0, 0.0

    points = coords.float()
    distances = torch.cdist(points, points, p=2.0)
    distances.fill_diagonal_(float("inf"))

    min_distances = torch.min(distances, dim=1).values
    valid = torch.isfinite(min_distances) & (min_distances > 0)

    radii = torch.zeros_like(min_distances)
    radii[valid] = min_distances[valid] / 2.0

    cube_volume = float((NF + 1) ** 3)
    sphere_volume = torch.sum((4.0 / 3.0) * math.pi * radii.pow(3))
    density = sphere_volume / cube_volume
    mean_diameter = 2.0 * torch.mean(radii)

    return density.item(), mean_diameter.item()


def compute_density_regular_sp(T: torch.Tensor) -> Tuple[float, float]:
    return compute_density_irregular_sp(T)


def compute_uniformity_error(
    C: torch.Tensor,
    M: int,
    normalized: bool = True,
) -> float:
    """
    RMSE of the three marginal sums relative to the desired value M.

    A value of 0 means perfect uniformity.
    """
    sum_x = torch.sum(C, dim=(1, 2))
    sum_y = torch.sum(C, dim=(0, 2))
    sum_z = torch.sum(C, dim=(0, 1))

    errors = torch.cat((sum_x - M, sum_y - M, sum_z - M)).float()
    rmse = torch.sqrt(torch.mean(errors.pow(2)))

    if normalized:
        rmse = rmse / max(float(M), 1.0)

    return rmse.item()


def evaluate(
    C: torch.Tensor,
    M: int,
    R: int,
    use_gram_condition: bool,
) -> Dict[str, float]:
    condition = get_cond(C, use_gram=use_gram_condition)
    if M == 1 and math.isfinite(condition):
        condition = 1.0

    if R == 1:
        density, diameter = compute_density_irregular_sp(C)
    else:
        density, diameter = compute_density_regular_sp(C)

    uniformity = compute_uniformity_error(C, M, normalized=True)

    return {
        "condition": condition,
        "density": density,
        "diameter": diameter,
        "uniformity": uniformity,
    }


# ============================================================
# Candidate moves
# ============================================================

def unique_distance_irregular_sp(C: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
    coords = torch.nonzero(C > 0).float()

    if coords.shape[0] < 2:
        zero = torch.tensor([0.0], device=C.device)
        return zero, zero

    distances = torch.cdist(coords, coords, p=2.0)
    distances.fill_diagonal_(float("inf"))
    min_distances = torch.min(distances, dim=1).values
    unique_min_distances = torch.unique(min_distances, sorted=True)
    return unique_min_distances, min_distances


def apply_plane_swap(
    C: torch.Tensor,
    i: int,
    ii: int,
    j: int,
) -> torch.Tensor:
    """Apply the original XZ/YZ plane-swap mutation."""
    NF = C.shape[2]
    candidate = C.clone()
    other_layer = int(torch.randint(0, NF, (1,), device=C.device).item())

    if other_layer == j:
        other_layer = (other_layer + 1) % NF

    if torch.rand(1, device=C.device).item() < 0.5:
        temp = candidate[i, :, j].clone()
        candidate[i, :, j] = candidate[i, :, other_layer]
        candidate[i, :, other_layer] = temp
    else:
        temp = candidate[:, ii, j].clone()
        candidate[:, ii, j] = candidate[:, ii, other_layer]
        candidate[:, ii, other_layer] = temp

    candidate.clamp_(0, 1)
    return candidate


def generate_candidate(C: torch.Tensor, cfg: Config) -> torch.Tensor:
    """Generate one candidate, usually targeting a point with minimum spacing."""
    coordinates = torch.nonzero(C > 0)
    if coordinates.shape[0] == 0:
        return C.clone()

    use_critical = torch.rand(1, device=C.device).item() < cfg.critical_move_probability

    if use_critical and coordinates.shape[0] >= 2:
        unique_distances, min_distances = unique_distance_irregular_sp(C)
        if unique_distances.numel() > 0:
            critical_positions = torch.nonzero(
                torch.isclose(min_distances, unique_distances[0])
            ).squeeze(1)
        else:
            critical_positions = torch.empty(0, dtype=torch.long, device=C.device)
    else:
        critical_positions = torch.empty(0, dtype=torch.long, device=C.device)

    if critical_positions.numel() > 0:
        random_position = torch.randint(
            0, critical_positions.numel(), (1,), device=C.device
        ).item()
        selected = int(critical_positions[random_position].item())
    else:
        selected = int(
            torch.randint(0, coordinates.shape[0], (1,), device=C.device).item()
        )

    # selected must be a Python integer. Keeping it as a one-element tensor
    # makes coordinates[selected] have shape (1, 3), which cannot be unpacked
    # into three scalar coordinates.
    i, ii, j = coordinates[selected].tolist()
    i, ii, j = int(i), int(ii), int(j)
    candidate = apply_plane_swap(C, i, ii, j)

    # Occasionally make one or two extra swaps to escape local optima.
    if torch.rand(1, device=C.device).item() < cfg.extra_random_moves_probability:
        n_moves = int(torch.randint(1, cfg.max_extra_random_moves + 1, (1,)).item())
        for _ in range(n_moves):
            coordinates = torch.nonzero(candidate > 0)
            selected = int(
                torch.randint(0, coordinates.shape[0], (1,), device=C.device).item()
            )
            i, ii, j = coordinates[selected].tolist()
            i, ii, j = int(i), int(ii), int(j)
            candidate = apply_plane_swap(candidate, i, ii, j)

    return candidate


# ============================================================
# Multi-objective comparison and Pareto archive
# ============================================================

def scheduled_value(start: float, end: float, progress: float) -> float:
    progress = min(max(progress, 0.0), 1.0)
    return start + (end - start) * progress


def constraint_violation(
    metrics: Dict[str, float],
    condition_limit: float,
    uniformity_limit: float,
) -> float:
    cond = metrics["condition"]
    uniformity = metrics["uniformity"]

    if not math.isfinite(cond):
        condition_violation = 1e6
    else:
        condition_violation = max(0.0, cond / max(condition_limit, 1e-12) - 1.0)

    uniformity_violation = max(0.0, uniformity - uniformity_limit)
    return condition_violation + uniformity_violation


def energy(
    metrics: Dict[str, float],
    refs: Dict[str, float],
    condition_limit: float,
    uniformity_limit: float,
    cfg: Config,
) -> float:
    density_norm = metrics["density"] / max(refs["density"], 1e-12)

    if math.isfinite(metrics["condition"]):
        condition_term = math.log(
            max(metrics["condition"], 1.0) / max(refs["condition"], 1.0)
        )
    else:
        condition_term = 1e6

    uniformity_norm = metrics["uniformity"] / max(refs["uniformity"], 1e-6)
    violation = constraint_violation(metrics, condition_limit, uniformity_limit)

    return (
        -cfg.weight_density * density_norm
        + cfg.weight_condition * condition_term
        + cfg.weight_uniformity * uniformity_norm
        + cfg.constraint_penalty * violation
    )


def is_feasible(
    metrics: Dict[str, float],
    condition_limit: float,
    uniformity_limit: float,
) -> bool:
    return (
        math.isfinite(metrics["condition"])
        and metrics["condition"] <= condition_limit
        and metrics["uniformity"] <= uniformity_limit
    )


def better_solution(
    a: Dict[str, float],
    b: Dict[str, float],
    condition_limit: float,
    uniformity_limit: float,
    density_tolerance: float = 1e-10,
) -> bool:
    """Feasibility-first, then density, then condition, then uniformity."""
    a_feasible = is_feasible(a, condition_limit, uniformity_limit)
    b_feasible = is_feasible(b, condition_limit, uniformity_limit)

    if a_feasible != b_feasible:
        return a_feasible

    if not a_feasible:
        a_violation = constraint_violation(a, condition_limit, uniformity_limit)
        b_violation = constraint_violation(b, condition_limit, uniformity_limit)
        if not math.isclose(a_violation, b_violation, rel_tol=1e-10, abs_tol=1e-12):
            return a_violation < b_violation

    if a["density"] > b["density"] + density_tolerance:
        return True
    if b["density"] > a["density"] + density_tolerance:
        return False

    if a["condition"] < b["condition"]:
        return True
    if b["condition"] < a["condition"]:
        return False

    return a["uniformity"] < b["uniformity"]


def dominates(a: Dict[str, float], b: Dict[str, float]) -> bool:
    """Pareto dominance: maximize density; minimize condition and uniformity."""
    no_worse = (
        a["density"] >= b["density"]
        and a["condition"] <= b["condition"]
        and a["uniformity"] <= b["uniformity"]
    )
    strictly_better = (
        a["density"] > b["density"]
        or a["condition"] < b["condition"]
        or a["uniformity"] < b["uniformity"]
    )
    return no_worse and strictly_better


def update_pareto_archive(
    archive: List[Tuple[torch.Tensor, Dict[str, float]]],
    candidate_C: torch.Tensor,
    candidate_metrics: Dict[str, float],
    max_size: int = 100,
) -> List[Tuple[torch.Tensor, Dict[str, float]]]:
    if not math.isfinite(candidate_metrics["condition"]):
        return archive

    for _, metrics in archive:
        if dominates(metrics, candidate_metrics):
            return archive

    archive = [
        (matrix, metrics)
        for matrix, metrics in archive
        if not dominates(candidate_metrics, metrics)
    ]
    archive.append((candidate_C.clone(), candidate_metrics.copy()))

    # Keep a diverse, bounded archive.
    if len(archive) > max_size:
        archive.sort(key=lambda item: item[1]["density"])
        indices = np.linspace(0, len(archive) - 1, max_size).round().astype(int)
        archive = [archive[index] for index in indices]

    return archive


# ============================================================
# Uniformization post-process
# ============================================================

def uniformize_matrix(
    C_in: torch.Tensor,
    M: int,
    max_iter_factor: int = 200,
) -> torch.Tensor:
    C_out = C_in.clone().float()
    N = C_out.shape[0]
    max_iter = max_iter_factor * N**2

    current_error = compute_uniformity_error(C_out, M, normalized=False)

    for _ in range(max_iter):
        if current_error == 0.0:
            break

        sum_x = torch.sum(C_out, dim=(1, 2))
        sum_y = torch.sum(C_out, dim=(0, 2))
        sum_z = torch.sum(C_out, dim=(0, 1))

        diff_x = sum_x - M
        diff_y = sum_y - M
        diff_z = sum_z - M

        is_row_ok = diff_x == 0
        is_col_ok = diff_y == 0
        is_layer_ok = diff_z == 0

        lock_mask = (
            is_row_ok[:, None, None]
            & is_col_ok[None, :, None]
            & is_layer_ok[None, None, :]
        )

        C_bool = C_out > 0
        movable_ones = torch.nonzero(C_bool & ~lock_mask)
        fillable_zeros = torch.nonzero((~C_bool) & ~lock_mask)

        if movable_ones.shape[0] == 0 or fillable_zeros.shape[0] == 0:
            break

        one_index = movable_ones[
            torch.randint(0, movable_ones.shape[0], (1,), device=C_out.device)
        ]
        zero_index = fillable_zeros[
            torch.randint(0, fillable_zeros.shape[0], (1,), device=C_out.device)
        ]

        candidate = C_out.clone()
        candidate[tuple(one_index.tolist())] = 0
        candidate[tuple(zero_index.tolist())] = 1

        candidate_error = compute_uniformity_error(candidate, M, normalized=False)
        if candidate_error <= current_error:
            C_out = candidate
            current_error = candidate_error

    return C_out


# ============================================================
# Optimization and plotting
# ============================================================

def optimize_matrix(
    C_initial: torch.Tensor,
    M: int,
    cfg: Config,
) -> Tuple[
    torch.Tensor,
    Dict[str, float],
    Dict[str, List[float]],
    List[Tuple[torch.Tensor, Dict[str, float]]],
]:
    current_C = C_initial.clone()
    current_metrics = evaluate(
        current_C, M, cfg.R, cfg.use_gram_condition
    )

    initial_uniformity_limit = max(
        current_metrics["uniformity"], cfg.uniformity_limit_start
    )
    condition_limit = max(
        current_metrics["condition"] * cfg.condition_limit_factor,
        current_metrics["condition"],
    )

    refs = {
        "density": max(current_metrics["density"], 1e-6),
        "condition": max(current_metrics["condition"], 1.0),
        "uniformity": max(current_metrics["uniformity"], 1e-3),
    }

    best_C = current_C.clone()
    best_metrics = current_metrics.copy()
    archive: List[Tuple[torch.Tensor, Dict[str, float]]] = []
    archive = update_pareto_archive(archive, current_C, current_metrics)

    history: Dict[str, List[float]] = {
        "candidate_condition": [],
        "candidate_density": [],
        "candidate_uniformity": [],
        "current_condition": [],
        "current_density": [],
        "current_uniformity": [],
        "best_condition": [],
        "best_density": [],
        "best_uniformity": [],
        "uniformity_limit": [],
        "temperature": [],
        "accepted": [],
    }

    no_improvement_count = 0

    for k in range(cfg.max_iterations):
        progress = k / max(cfg.max_iterations - 1, 1)
        uniformity_limit = scheduled_value(
            initial_uniformity_limit,
            cfg.uniformity_limit_final,
            progress,
        )
        temperature = scheduled_value(
            cfg.initial_temperature,
            cfg.final_temperature,
            progress,
        )

        candidate_C = generate_candidate(current_C, cfg)
        candidate_metrics = evaluate(
            candidate_C, M, cfg.R, cfg.use_gram_condition
        )

        current_energy = energy(
            current_metrics,
            refs,
            condition_limit,
            uniformity_limit,
            cfg,
        )
        candidate_energy = energy(
            candidate_metrics,
            refs,
            condition_limit,
            uniformity_limit,
            cfg,
        )

        delta = candidate_energy - current_energy
        accept = delta <= 0
        if not accept and temperature > 0 and math.isfinite(delta):
            probability = math.exp(-min(delta / temperature, 700.0))
            accept = random.random() < probability

        if accept:
            current_C = candidate_C
            current_metrics = candidate_metrics

        archive = update_pareto_archive(archive, candidate_C, candidate_metrics)

        if better_solution(
            current_metrics,
            best_metrics,
            condition_limit,
            uniformity_limit,
        ):
            best_C = current_C.clone()
            best_metrics = current_metrics.copy()
            no_improvement_count = 0
        else:
            no_improvement_count += 1

        history["candidate_condition"].append(candidate_metrics["condition"])
        history["candidate_density"].append(candidate_metrics["density"])
        history["candidate_uniformity"].append(candidate_metrics["uniformity"])
        history["current_condition"].append(current_metrics["condition"])
        history["current_density"].append(current_metrics["density"])
        history["current_uniformity"].append(current_metrics["uniformity"])
        history["best_condition"].append(best_metrics["condition"])
        history["best_density"].append(best_metrics["density"])
        history["best_uniformity"].append(best_metrics["uniformity"])
        history["uniformity_limit"].append(uniformity_limit)
        history["temperature"].append(temperature)
        history["accepted"].append(float(accept))

        if (k + 1) % cfg.print_every == 0:
            print(
                f"Iter {k + 1:5d} | "
                f"Current D={current_metrics['density']:.6f}, "
                f"C={current_metrics['condition']:.4f}, "
                f"U={current_metrics['uniformity']:.5f} | "
                f"Best D={best_metrics['density']:.6f}, "
                f"C={best_metrics['condition']:.4f}, "
                f"U={best_metrics['uniformity']:.5f} | "
                f"No improvement {no_improvement_count}/{cfg.patience}"
            )

        if k >= cfg.burn_in and no_improvement_count >= cfg.patience:
            print(f">> Convergence reached at iteration {k + 1}.")
            break

    return best_C, best_metrics, history, archive


def plot_histories(
    history: Dict[str, List[float]],
    archive: List[Tuple[torch.Tensor, Dict[str, float]]],
    M: int,
    cfg: Config,
) -> None:
    if not history["current_density"]:
        return

    iterations = np.arange(1, len(history["current_density"]) + 1)

    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    fig.suptitle(f"Multi-objective Optimization Evolution for M={M}", fontsize=16)

    axes[0].plot(iterations, history["candidate_condition"], alpha=0.25, label="Candidate")
    axes[0].plot(iterations, history["current_condition"], label="Current accepted")
    axes[0].plot(iterations, history["best_condition"], label="Best stored")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Condition number")
    axes[0].grid(True)
    axes[0].legend()

    axes[1].plot(iterations, history["candidate_density"], alpha=0.25, label="Candidate")
    axes[1].plot(iterations, history["current_density"], label="Current accepted")
    axes[1].plot(iterations, history["best_density"], label="Best stored")
    axes[1].set_ylabel("Density")
    axes[1].grid(True)
    axes[1].legend()

    axes[2].plot(iterations, history["candidate_uniformity"], alpha=0.25, label="Candidate")
    axes[2].plot(iterations, history["current_uniformity"], label="Current accepted")
    axes[2].plot(iterations, history["best_uniformity"], label="Best stored")
    axes[2].plot(iterations, history["uniformity_limit"], linestyle="--", label="Uniformity limit")
    axes[2].set_xlabel("Iterations")
    axes[2].set_ylabel("Normalized uniformity RMSE")
    axes[2].grid(True)
    axes[2].legend()

    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(os.path.join(cfg.graph_export_path, f"Evolution_Mux_{M}.png"), dpi=160)
    plt.close(fig)

    if archive:
        conditions = np.array([item[1]["condition"] for item in archive])
        densities = np.array([item[1]["density"] for item in archive])
        uniformities = np.array([item[1]["uniformity"] for item in archive])

        fig = plt.figure(figsize=(8, 6))
        scatter = plt.scatter(conditions, densities, c=uniformities)
        plt.xscale("log")
        plt.xlabel("Condition number")
        plt.ylabel("Density")
        plt.title(f"Pareto Archive for M={M}")
        plt.grid(True)
        plt.colorbar(scatter, label="Normalized uniformity RMSE")
        fig.tight_layout()
        fig.savefig(os.path.join(cfg.graph_export_path, f"Pareto_Mux_{M}.png"), dpi=160)
        plt.close(fig)


def history_to_numpy(history: Dict[str, List[float]]) -> Dict[str, np.ndarray]:
    return {key: np.asarray(value, dtype=np.float64) for key, value in history.items()}


# ============================================================
# Main
# ============================================================

def main() -> None:
    cfg = Config()

    random.seed(cfg.random_seed)
    np.random.seed(cfg.random_seed)
    torch.manual_seed(cfg.random_seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(cfg.random_seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    os.makedirs(cfg.output_folder, exist_ok=True)
    if cfg.plot_graphs:
        os.makedirs(cfg.graph_export_path, exist_ok=True)

    with open(cfg.log_filename, "w", encoding="utf-8") as file_id:
        file_id.write(f"Execution Summary - {datetime.datetime.now()}\n\n")
        header = (
            f"{'M':<4} | {'Initial Cond':<13} | {'Optimized Cond':<15} | "
            f"{'Final Cond':<12} | {'Initial Density':<15} | "
            f"{'Optimized Density':<17} | {'Final Density':<13} | "
            f"{'Initial U':<10} | {'Optimized U':<11} | {'Final U':<10}\n"
        )
        file_id.write(header)
        file_id.write("-" * (len(header) - 1) + "\n")

        for M in cfg.mux_values:
            print(f"\n{'=' * 72}\nStarting M={M}\n{'=' * 72}")

            G = run_eq10mux5(cfg.NF, M, device)
            C_initial, _ = generate_coded_aperture(
                G, cfg.NF, M, cfg.NF, device
            )
            initial_metrics = evaluate(
                C_initial, M, cfg.R, cfg.use_gram_condition
            )

            if M == 1:
                C_optimized = C_initial.clone()
                optimized_metrics = initial_metrics.copy()
                history = {
                    key: []
                    for key in (
                        "candidate_condition", "candidate_density", "candidate_uniformity",
                        "current_condition", "current_density", "current_uniformity",
                        "best_condition", "best_density", "best_uniformity",
                        "uniformity_limit", "temperature", "accepted",
                    )
                }
                archive = [(C_initial.clone(), initial_metrics.copy())]
            else:
                C_optimized, optimized_metrics, history, archive = optimize_matrix(
                    C_initial, M, cfg
                )

            if cfg.apply_final_uniformization and M > 1:
                print(">> Applying final uniformization...")
                C_final = uniformize_matrix(
                    C_optimized,
                    M,
                    max_iter_factor=cfg.uniformizer_max_iter_factor,
                )
            else:
                C_final = C_optimized.clone()

            final_metrics = evaluate(
                C_final, M, cfg.R, cfg.use_gram_condition
            )

            if cfg.plot_graphs and M > 1:
                plot_histories(history, archive, M, cfg)

            pareto_conditions = np.asarray(
                [item[1]["condition"] for item in archive], dtype=np.float64
            )
            pareto_densities = np.asarray(
                [item[1]["density"] for item in archive], dtype=np.float64
            )
            pareto_uniformities = np.asarray(
                [item[1]["uniformity"] for item in archive], dtype=np.float64
            )

            filename = os.path.join(
                cfg.output_folder,
                f"Bands_{cfg.NF}_Mux_{M}.mat",
            )
            scipy.io.savemat(
                filename,
                {
                    "C_initial": C_initial.detach().cpu().numpy(),
                    "C_optimized": C_optimized.detach().cpu().numpy(),
                    "C_final": C_final.detach().cpu().numpy(),
                    "initial_metrics": initial_metrics,
                    "optimized_metrics": optimized_metrics,
                    "final_metrics": final_metrics,
                    "history": history_to_numpy(history),
                    "pareto_condition": pareto_conditions,
                    "pareto_density": pareto_densities,
                    "pareto_uniformity": pareto_uniformities,
                },
            )

            print(f"-------------------- Summary M={M} --------------------")
            print(
                f"Condition : {initial_metrics['condition']:.4f} -> "
                f"{optimized_metrics['condition']:.4f} -> {final_metrics['condition']:.4f}"
            )
            print(
                f"Density   : {initial_metrics['density']:.6f} -> "
                f"{optimized_metrics['density']:.6f} -> {final_metrics['density']:.6f}"
            )
            print(
                f"Uniformity: {initial_metrics['uniformity']:.6f} -> "
                f"{optimized_metrics['uniformity']:.6f} -> {final_metrics['uniformity']:.6f}"
            )
            print(f"Pareto archive size: {len(archive)}")
            print(f"Saved: {filename}")

            row = (
                f"{M:<4} | "
                f"{initial_metrics['condition']:<13.4f} | "
                f"{optimized_metrics['condition']:<15.4f} | "
                f"{final_metrics['condition']:<12.4f} | "
                f"{initial_metrics['density']:<15.6f} | "
                f"{optimized_metrics['density']:<17.6f} | "
                f"{final_metrics['density']:<13.6f} | "
                f"{initial_metrics['uniformity']:<10.6f} | "
                f"{optimized_metrics['uniformity']:<11.6f} | "
                f"{final_metrics['uniformity']:<10.6f}\n"
            )
            file_id.write(row)
            file_id.flush()


if __name__ == "__main__":
    main()