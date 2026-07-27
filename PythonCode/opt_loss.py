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
    """Número de condición de sigma.T @ sigma."""
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


def normalized_condition(cond_value, initial_cond):
    """
    Normaliza el número de condición para que sea comparable con 0.74-density.

    Se usa log(1+cond), porque el número de condición puede variar por varios
    órdenes de magnitud y, sin normalización, dominaría completamente la loss.
    """
    if not math.isfinite(cond_value):
        return float("inf")

    denominator = math.log1p(max(initial_cond, 1.0))
    return math.log1p(max(cond_value, 1.0)) / max(denominator, 1e-12)


def combined_loss(density, cond_value, initial_cond, density_target=0.74, cond_weight=1.0):
    """
    Loss solicitada:

        loss = (0.74 - density) + condition_number

    En la implementación, condition_number se usa normalizado para balancear
    las escalas:

        loss = max(0, 0.74-density) + cond_weight * cond_normalizado

    Minimizar la loss aumenta la densidad y reduce el número de condición.
    """
    density_gap = max(0.0, density_target - density)
    cond_term = normalized_condition(cond_value, initial_cond)
    return density_gap + cond_weight * cond_term


def propose_plane_swap(C, critical_point):
    """Propuesta original: intercambio de un plano XZ o YZ."""
    NF = C.shape[2]
    i, ii, j = [int(v.item()) for v in critical_point]
    C_next = C.clone()

    other_layer = int(torch.randint(0, NF, (1,), device=C.device).item())
    if other_layer == j:
        other_layer = (other_layer + 1) % NF

    if torch.rand(1, device=C.device).item() < 0.5:
        # Intercambio de la línea C[i, :, :] entre dos capas.
        temp = C_next[i, :, j].clone()
        C_next[i, :, j] = C_next[i, :, other_layer]
        C_next[i, :, other_layer] = temp
    else:
        # Intercambio de la línea C[:, ii, :] entre dos capas.
        temp = C_next[:, ii, j].clone()
        C_next[:, ii, j] = C_next[:, ii, other_layer]
        C_next[:, ii, other_layer] = temp

    C_next.clamp_(0, 1)
    return C_next


def save_history_plot(history, M, export_path):
    if not history["loss"]:
        return

    iterations = range(1, len(history["loss"]) + 1)
    fig, axes = plt.subplots(3, 1, figsize=(10, 10))
    fig.suptitle(f"Plane-swap optimization, M={M}", fontsize=15)

    axes[0].plot(iterations, history["loss"])
    axes[0].set_ylabel("Loss")
    axes[0].grid(True)

    axes[1].plot(iterations, history["density"])
    axes[1].set_ylabel("Density")
    axes[1].grid(True)

    axes[2].plot(iterations, history["condition"])
    axes[2].set_yscale("log")
    axes[2].set_ylabel("Condition number")
    axes[2].set_xlabel("Iteration")
    axes[2].grid(True)

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(os.path.join(export_path, f"loss_plane_M_{M}.png"), dpi=160)
    plt.close(fig)


def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    # ---------------------------- Configuración -----------------------------
    NF = 31
    max_iterations = 10000
    patience = 500
    burn_in = 10
    rel_tol = 1e-7

    density_target = 0.74
    cond_weight = 1.0

    plot_graphs = True
    graph_export_path = "./PythonCode/graphs_loss_density_condition"
    output_folder = "./PythonCode/Patrones/LossDensityCondition"
    log_filename = "resumen_loss_density_condition.txt"
    mux_values_to_run = range(1, 17)

    os.makedirs(output_folder, exist_ok=True)
    if plot_graphs:
        os.makedirs(graph_export_path, exist_ok=True)

    with open(log_filename, "w", encoding="utf-8") as log_file:
        log_file.write(f"Execution Summary - {datetime.datetime.now()}\n\n")
        log_file.write(
            f"{'M':<4} | {'Initial Cond':<14} | {'Best Cond':<14} | "
            f"{'Initial Density':<17} | {'Best Density':<14} | "
            f"{'Initial Loss':<14} | {'Best Loss':<14}\n"
        )
        log_file.write("-" * 115 + "\n")

        for M in mux_values_to_run:
            G = run_eq10mux5(NF, M, device)
            C, _ = generate_coded_aperture(G, NF, M, NF, device)
            C_initial = C.clone()

            initial_cond = 1.0 if M == 1 else get_cond(C)
            initial_density, _ = compute_density_irregular_sp(C)
            initial_loss = combined_loss(
                initial_density,
                initial_cond,
                initial_cond,
                density_target,
                cond_weight,
            )

            best_C = C.clone()
            best_cond = initial_cond
            best_density = initial_density
            best_loss = initial_loss
            no_improvement_count = 0

            history = {"loss": [], "density": [], "condition": []}

            print(
                f"\nM={M} | initial loss={initial_loss:.6f} | "
                f"density={initial_density:.6f} | cond={initial_cond:.6g}"
            )

            if M > 1:
                for k in range(max_iterations):
                    unique_distances, point_min_distances = unique_distance_irregular_sp(C)
                    occupied_points = torch.nonzero(C > 0)

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
                    C_next = propose_plane_swap(C, critical_point)

                    new_cond = get_cond(C_next)
                    new_density, _ = compute_density_irregular_sp(C_next)
                    new_loss = combined_loss(
                        new_density,
                        new_cond,
                        initial_cond,
                        density_target,
                        cond_weight,
                    )

                    history["loss"].append(new_loss)
                    history["density"].append(new_density)
                    history["condition"].append(new_cond)

                    improvement = best_loss - new_loss
                    required_improvement = rel_tol * max(abs(best_loss), 1.0)

                    if math.isfinite(new_loss) and improvement > required_improvement:
                        C = C_next
                        best_C = C_next.clone()
                        best_cond = new_cond
                        best_density = new_density
                        best_loss = new_loss
                        no_improvement_count = 0
                    else:
                        no_improvement_count += 1

                    if (k + 1) % 20 == 0:
                        print(
                            f"Iter {k + 1:4d} | loss={new_loss:.6f} | "
                            f"best={best_loss:.6f} | dens={new_density:.6f} | "
                            f"cond={new_cond:.5g} | no-improved="
                            f"{no_improvement_count}/{patience}"
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
                    "initial_loss": initial_loss,
                    "optimized_loss": best_loss,
                },
            )

            if plot_graphs:
                save_history_plot(history, M, graph_export_path)

            print(
                f"Summary M={M}: loss {initial_loss:.6f} -> {best_loss:.6f}; "
                f"density {initial_density:.6f} -> {best_density:.6f}; "
                f"cond {initial_cond:.6g} -> {best_cond:.6g}"
            )

            log_file.write(
                f"{M:<4} | {initial_cond:<14.6g} | {best_cond:<14.6g} | "
                f"{initial_density:<17.6f} | {best_density:<14.6f} | "
                f"{initial_loss:<14.6f} | {best_loss:<14.6f}\n"
            )


if __name__ == "__main__":
    main()
