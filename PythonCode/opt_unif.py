import datetime
import math
import os

import matplotlib.pyplot as plt
import scipy.io
import torch


# -----------------------------------------------------------------------------
# Reemplaza esta función por tu implementación real de RunEQ10MUX5.
# -----------------------------------------------------------------------------
def run_eq10mux5(NF, M, device):
    """Implementación provisional: genera una matriz aleatoria G."""
    return torch.randint(1, NF + 1, (NF, NF, M), device=device)


def generate_coded_aperture(G, NF, B, N, device):
    C = torch.zeros((NF, NF, NF), device=device)
    T = torch.zeros((N, N, NF), device=device)

    for j in range(1, NF + 1):
        for i in range(B):
            C[:, :, j - 1] += (G[:, :, i] == j).float()

    C.clamp_(0, 1)
    return C, T


def get_cond(C):
    NF = C.shape[2]
    sigma = torch.stack([C[:, :, j].flatten() for j in range(NF)], dim=1).float()

    max_val = sigma.max()
    if max_val > 0:
        sigma = sigma / max_val

    cond = torch.linalg.cond(sigma.T @ sigma)
    value = float(cond.item())
    return value if math.isfinite(value) else float("inf")


def compute_density_irregular_sp(T):
    NF = T.shape[2]
    coords = torch.nonzero(T > 0)

    if coords.size(0) < 2:
        return 0.0, 0.0

    X = coords.T.float()
    U, S, Vh = torch.linalg.svd(X, full_matrices=False)
    Q = torch.round(U @ torch.diag(S) @ Vh)

    distances = torch.cdist(Q.T, Q.T, p=2.0)
    distances.fill_diagonal_(float("inf"))
    min_distances = torch.min(distances, dim=1).values

    radii = torch.zeros(X.shape[1], device=T.device)
    valid = torch.isfinite(min_distances) & (min_distances > 0)
    radii[valid] = min_distances[valid] / 2.0

    cube_volume = (NF + 1) ** 3
    sphere_volume = torch.sum((4.0 / 3.0) * math.pi * radii**3)
    density = sphere_volume / cube_volume
    diameter = 2.0 * torch.mean(radii)

    return float(density.item()), float(diameter.item())


def unique_distance_irregular_sp(C):
    coords = torch.nonzero(C > 0).float()

    if coords.size(0) < 2:
        zero = torch.tensor([0.0], device=C.device)
        return zero, zero

    distances = torch.cdist(coords, coords, p=2.0)
    distances.fill_diagonal_(float("inf"))
    min_distances = torch.min(distances, dim=1).values
    return torch.unique(min_distances), min_distances


def uniformity_rmse(C, M):
    """
    RMSE de uniformidad de las proyecciones sobre X, Y y Z.

    El valor ideal es cero:
      sum(C, ejes YZ) = M para todo X
      sum(C, ejes XZ) = M para todo Y
      sum(C, ejes XY) = M para todo Z
    """
    sum_x = torch.sum(C, dim=(1, 2))
    sum_y = torch.sum(C, dim=(0, 2))
    sum_z = torch.sum(C, dim=(0, 1))
    errors = torch.cat((sum_x - M, sum_y - M, sum_z - M)).float()
    return float(torch.sqrt(torch.mean(errors**2)).item())


def normalized_condition(cond_value, initial_cond):
    if not math.isfinite(cond_value):
        return float("inf")

    denominator = math.log1p(max(initial_cond, 1.0))
    return math.log1p(max(cond_value, 1.0)) / max(denominator, 1e-12)


def combined_loss(
    density,
    cond_value,
    uniformity_value,
    initial_cond,
    M,
    density_target=0.74,
    cond_weight=1.0,
    uniformity_weight=1.0,
):
    """
    Loss de tres objetivos:

      loss = max(0, 0.74-density)
             + cond_weight * cond_normalizado
             + uniformity_weight * uniformity_normalizada

    La uniformidad se normaliza por M para hacer comparables distintos M.
    """
    density_gap = max(0.0, density_target - density)
    cond_term = normalized_condition(cond_value, initial_cond)
    uniformity_term = uniformity_value / max(float(M), 1.0)

    return (
        density_gap
        + cond_weight * cond_term
        + uniformity_weight * uniformity_term
    )


def propose_point_move(C, critical_point, M, candidate_count=64):
    """
    Mueve UN solo punto, no un plano.

    1. Elimina el 1 correspondiente a un punto de distancia mínima.
    2. Evalúa hasta `candidate_count` posiciones vacías.
    3. Prefiere posiciones pertenecientes a planos X/Y/Z deficitarios.
    4. Inserta el punto en una única posición vacía.

    El número total de unos permanece constante.
    """
    C_next = C.clone()
    source = tuple(int(v.item()) for v in critical_point)
    C_next[source] = 0.0

    empty_points = torch.nonzero(C_next == 0, as_tuple=False)
    if empty_points.numel() == 0:
        return C.clone()

    sample_size = min(candidate_count, empty_points.shape[0])
    indices = torch.randperm(empty_points.shape[0], device=C.device)[:sample_size]
    candidates = empty_points[indices]

    # Déficit después de retirar el punto. Un valor positivo indica que el
    # plano necesita recibir puntos para acercarse al objetivo M.
    deficit_x = M - torch.sum(C_next, dim=(1, 2))
    deficit_y = M - torch.sum(C_next, dim=(0, 2))
    deficit_z = M - torch.sum(C_next, dim=(0, 1))

    scores = (
        deficit_x[candidates[:, 0]]
        + deficit_y[candidates[:, 1]]
        + deficit_z[candidates[:, 2]]
    )

    # Se escoge aleatoriamente entre los candidatos con mejor puntuación para
    # conservar diversidad durante la búsqueda.
    best_score = torch.max(scores)
    best_indices = torch.nonzero(torch.isclose(scores, best_score)).squeeze(1)
    selected_local = int(
        best_indices[
            torch.randint(0, best_indices.numel(), (1,), device=C.device)
        ].item()
    )
    destination = tuple(int(v.item()) for v in candidates[selected_local])
    C_next[destination] = 1.0

    return C_next


def save_history_plot(history, M, export_path):
    if not history["loss"]:
        return

    iterations = range(1, len(history["loss"]) + 1)
    fig, axes = plt.subplots(4, 1, figsize=(10, 12))
    fig.suptitle(f"Point-move optimization, M={M}", fontsize=15)

    axes[0].plot(iterations, history["loss"])
    axes[0].set_ylabel("Loss")
    axes[0].grid(True)

    axes[1].plot(iterations, history["density"])
    axes[1].set_ylabel("Density")
    axes[1].grid(True)

    axes[2].plot(iterations, history["condition"])
    axes[2].set_yscale("log")
    axes[2].set_ylabel("Condition number")
    axes[2].grid(True)

    axes[3].plot(iterations, history["uniformity"])
    axes[3].set_ylabel("Uniformity RMSE")
    axes[3].set_xlabel("Iteration")
    axes[3].grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(os.path.join(export_path, f"loss_point_M_{M}.png"), dpi=160)
    plt.close(fig)


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # ---------------------------- Configuración -----------------------------
    NF = 31
    max_iterations = 20000
    patience = 500
    burn_in = 10
    rel_tol = 1e-7

    density_target = 0.74
    cond_weight = 1.0
    uniformity_weight = 1.0
    point_candidate_count = 64

    plot_graphs = True
    graph_export_path = "./PythonCode/graphs_loss_point_uniformity"
    output_folder = "./PythonCode/Patrones/LossPointUniformity"
    log_filename = "resumen_loss_point_uniformity.txt"
    mux_values_to_run = range(1, 17)

    os.makedirs(output_folder, exist_ok=True)
    if plot_graphs:
        os.makedirs(graph_export_path, exist_ok=True)

    with open(log_filename, "w", encoding="utf-8") as log_file:
        log_file.write(f"Execution Summary - {datetime.datetime.now()}\n\n")
        log_file.write(
            f"{'M':<4} | {'Initial Cond':<14} | {'Best Cond':<14} | "
            f"{'Initial Density':<17} | {'Best Density':<14} | "
            f"{'Initial Uniform':<16} | {'Best Uniform':<14} | "
            f"{'Initial Loss':<14} | {'Best Loss':<14}\n"
        )
        log_file.write("-" * 155 + "\n")

        for M in mux_values_to_run:
            G = run_eq10mux5(NF, M, device)
            C, _ = generate_coded_aperture(G, NF, M, NF, device)
            C_initial = C.clone()

            initial_cond = 1.0 if M == 1 else get_cond(C)
            initial_density, _ = compute_density_irregular_sp(C)
            initial_uniformity = uniformity_rmse(C, M)
            initial_loss = combined_loss(
                initial_density,
                initial_cond,
                initial_uniformity,
                initial_cond,
                M,
                density_target,
                cond_weight,
                uniformity_weight,
            )

            best_C = C.clone()
            best_cond = initial_cond
            best_density = initial_density
            best_uniformity = initial_uniformity
            best_loss = initial_loss
            no_improvement_count = 0

            history = {
                "loss": [],
                "density": [],
                "condition": [],
                "uniformity": [],
            }

            print(
                f"\nM={M} | initial loss={initial_loss:.6f} | "
                f"density={initial_density:.6f} | cond={initial_cond:.6g} | "
                f"uniformity={initial_uniformity:.6f}"
            )

            if M > 1:
                for k in range(max_iterations):
                    unique_distances, point_min_distances = unique_distance_irregular_sp(C)
                    occupied_points = torch.nonzero(C > 0, as_tuple=False)

                    if occupied_points.size(0) < 2:
                        break

                    minimum_distance = unique_distances[0]
                    critical_positions = torch.nonzero(
                        torch.isclose(point_min_distances, minimum_distance),
                        as_tuple=False,
                    ).squeeze(1)

                    if critical_positions.numel() == 0:
                        continue

                    selected = int(
                        torch.randint(
                            0,
                            critical_positions.numel(),
                            (1,),
                            device=C.device,
                        ).item()
                    )
                    critical_point = occupied_points[critical_positions[selected]]

                    C_next = propose_point_move(
                        C,
                        critical_point,
                        M,
                        candidate_count=point_candidate_count,
                    )

                    new_cond = get_cond(C_next)
                    new_density, _ = compute_density_irregular_sp(C_next)
                    new_uniformity = uniformity_rmse(C_next, M)
                    new_loss = combined_loss(
                        new_density,
                        new_cond,
                        new_uniformity,
                        initial_cond,
                        M,
                        density_target,
                        cond_weight,
                        uniformity_weight,
                    )

                    history["loss"].append(new_loss)
                    history["density"].append(new_density)
                    history["condition"].append(new_cond)
                    history["uniformity"].append(new_uniformity)

                    improvement = best_loss - new_loss
                    required_improvement = rel_tol * max(abs(best_loss), 1.0)

                    if math.isfinite(new_loss) and improvement > required_improvement:
                        C = C_next
                        best_C = C_next.clone()
                        best_cond = new_cond
                        best_density = new_density
                        best_uniformity = new_uniformity
                        best_loss = new_loss
                        no_improvement_count = 0
                    else:
                        no_improvement_count += 1

                    if (k + 1) % 20 == 0:
                        print(
                            f"Iter {k + 1:4d} | loss={new_loss:.6f} | "
                            f"best={best_loss:.6f} | dens={new_density:.6f} | "
                            f"cond={new_cond:.5g} | unif={new_uniformity:.5f} | "
                            f"no-improved={no_improvement_count}/{patience}"
                        )

                    if k >= burn_in and no_improvement_count >= patience:
                        print(f">> Convergence reached at iteration {k + 1}.")
                        break

            C_optimized = best_C

            scipy.io.savemat(
                os.path.join(output_folder, f"Bands_{NF}_Mux_{M}.mat"),
                {
                    "C_initial": C_initial.cpu().numpy(),
                    "C_optimized": C_optimized.cpu().numpy(),
                    "initial_condition": initial_cond,
                    "optimized_condition": best_cond,
                    "initial_density": initial_density,
                    "optimized_density": best_density,
                    "initial_uniformity_rmse": initial_uniformity,
                    "optimized_uniformity_rmse": best_uniformity,
                    "initial_loss": initial_loss,
                    "optimized_loss": best_loss,
                },
            )

            if plot_graphs:
                save_history_plot(history, M, graph_export_path)

            print(
                f"Summary M={M}: loss {initial_loss:.6f} -> {best_loss:.6f}; "
                f"density {initial_density:.6f} -> {best_density:.6f}; "
                f"cond {initial_cond:.6g} -> {best_cond:.6g}; "
                f"uniformity {initial_uniformity:.6f} -> {best_uniformity:.6f}"
            )

            log_file.write(
                f"{M:<4} | {initial_cond:<14.6g} | {best_cond:<14.6g} | "
                f"{initial_density:<17.6f} | {best_density:<14.6f} | "
                f"{initial_uniformity:<16.6f} | {best_uniformity:<14.6f} | "
                f"{initial_loss:<14.6f} | {best_loss:<14.6f}\n"
            )


if __name__ == "__main__":
    main()
