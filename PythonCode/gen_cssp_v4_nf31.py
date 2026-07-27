import torch
import numpy as np
import os
import datetime
import math
import scipy.io
import matplotlib.pyplot as plt

# --- Dummy function for RunEQ10MUX5 ---
def run_eq10mux5(NF, M, device):
    """
    Dummy implementation of RunEQ10MUX5.
    Returns a random generator matrix G.
    """
    return torch.randint(1, NF + 1, (NF, NF, M), device=device)

# --- Translated generateCodedAperture ---
def generate_coded_aperture(G, NF, B, N, device):
    N1 = int(math.ceil(N / NF))
    # II = torch.ones((N1, N1), device=device)
    C = torch.zeros((NF, NF, NF), device=device)
    T = torch.zeros((N, N, NF), device=device)

    for j in range(1, NF + 1):
        for i in range(B):
            temp = (G[:, :, i] == j).float()
            C[:, :, j - 1] += temp

    return C, T

# --- Translated GetCond ---
def get_cond(C):
    NF = C.shape[2]
    # Flatten de cada capa y concatenación
    cols = [C[:, :, j].flatten().unsqueeze(1) for j in range(NF)]
    sigma = torch.cat(cols, dim=1).float()
    
    max_val = sigma.max()
    if max_val > 0:
        sigma = sigma / max_val
        
    return torch.linalg.cond(sigma.T @ sigma).item()

# --- Translated ComputeDensityIrregularSP ---
def compute_density_irregular_sp(T):
    NF = T.shape[2]
    
    coords = torch.nonzero(T > 0)
    if coords.size(0) == 0:
        return 0.0, 0.0
        
    X = coords.T.float()

    # SVD completamente en la GPU
    U, S, Vh = torch.linalg.svd(X, full_matrices=False)
    S_mat = torch.diag(S)
    
    Q = torch.round(U @ S_mat @ Vh)
    
    # Distancias en GPU (Reemplaza a pdist de scipy)
    D1 = torch.cdist(Q.T, Q.T, p=2.0)
    D1.fill_diagonal_(float('inf'))
    
    min_dists, _ = torch.min(D1, dim=1)
    
    R = torch.zeros(X.shape[1], device=T.device)
    mask = (min_dists != float('inf')) & (min_dists > 0)
    if mask.any():
        R[mask] = min_dists[mask] / 2.0

    cube = (NF + 1) ** 3
    sphere_vol = torch.sum((4/3) * math.pi * (R ** 3))
    density = sphere_vol / cube
    diameter = torch.mean(R) * 2.0
    
    return density.item(), diameter.item()

# --- Translated ComputeDensityRegularSP ---
def compute_density_regular_sp(T):
    return compute_density_irregular_sp(T)

# --- Translated uniqueDistanceIrregularSP ---
def unique_distance_irregular_sp(C):
    coords = torch.nonzero(C > 0).float()
    
    if coords.size(0) < 2:
        return torch.tensor([0.0], device=C.device), torch.tensor([0.0], device=C.device)

    # Distancias en GPU
    dists_xp = torch.cdist(coords, coords, p=2.0)
    dists_xp.fill_diagonal_(float('inf'))
    
    min_dists, _ = torch.min(dists_xp, dim=1)
    unique_min_dists = torch.unique(min_dists)
    
    return unique_min_dists, min_dists

# --- Translated uniformizeMatrix ---
def uniformize_matrix(C_in, M):
    C_out = C_in.clone().float()
    N, _, _ = C_out.shape
    max_iter = 200 * N**2

    for _ in range(max_iter):
        sum_x = torch.sum(C_out, dim=(1, 2))
        sum_y = torch.sum(C_out, dim=(0, 2))
        sum_z = torch.sum(C_out, dim=(0, 1))

        diff_x = sum_x - M
        diff_y = sum_y - M
        diff_z = sum_z - M

        if not torch.any(diff_x) and not torch.any(diff_y) and not torch.any(diff_z):
            return C_out

        is_row_ok = (diff_x == 0)
        is_col_ok = (diff_y == 0)
        is_layer_ok = (diff_z == 0)

        LockMask = (is_row_ok[:, None, None] & 
                    is_col_ok[None, :, None] & 
                    is_layer_ok[None, None, :])

        C_bool = C_out > 0
        movable_ones_indices = torch.nonzero(C_bool & ~LockMask)
        fillable_zeros_indices = torch.nonzero((~C_bool) & ~LockMask)

        if len(movable_ones_indices) == 0 or len(fillable_zeros_indices) == 0:
            print("Uniformizer: No valid moves found.")
            return C_out

        rand_one_idx = torch.randint(0, len(movable_ones_indices), (1,)).item()
        rand_zero_idx = torch.randint(0, len(fillable_zeros_indices), (1,)).item()

        one_to_move = tuple(movable_ones_indices[rand_one_idx].tolist())
        zero_to_fill = tuple(fillable_zeros_indices[rand_zero_idx].tolist())

        C_temp = C_out.clone()
        C_temp[one_to_move] = 0
        C_temp[zero_to_fill] = 1
        
        sum_x_new = torch.sum(C_temp, dim=(1, 2))
        sum_y_new = torch.sum(C_temp, dim=(0, 2))
        sum_z_new = torch.sum(C_temp, dim=(0, 1))

        errors_new = torch.cat([sum_x_new - M, sum_y_new - M, sum_z_new - M])
        new_rmse = torch.sqrt(torch.mean(errors_new.float()**2))
        
        errors_old = torch.cat([diff_x, diff_y, diff_z])
        current_rmse = torch.sqrt(torch.mean(errors_old.float()**2))

        if new_rmse <= current_rmse:
            C_out = C_temp

    print("Uniformizer: Max iterations reached.")
    return C_out

def main():
    # --- Check for GPU ---
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")
    
    # --- Configuration Parameters ---
    NF = 31
    N = NF
    R = 1
    max_iterations = 100000
    rel_tol = 1e-5
    patience = 2000
    optimization_metric = 'density'
    burn_in = 10
    plot_graphs = True
    graph_export_path = './PythonCode/graphs_pt'

    if plot_graphs:
        os.makedirs(graph_export_path, exist_ok=True)
    
    # --- Logging setup ---
    log_filename = 'resumen_metricas.txt'
    with open(log_filename, 'w') as fileID:
        line = 'Execution Summary - ' + str(datetime.datetime.now()) + '\n\n'
        fileID.write(line)

        header = (f"{'M':<4} | {'Initial Cond':<14} | {'Optimized Cond':<16} | {'Uniform Cond':<14} | "
                  f"{'Initial Density':<17} | {'Optimized Density':<19} | {'Uniform Density':<17}\n")
        separator = '-' * (len(header) - 1) + '\n'
        fileID.write(header)
        fileID.write(separator)

        mux_values_to_run = range(1, 17)

        for M in mux_values_to_run:
            
            # --- Matrix Initialization ---
            G = run_eq10mux5(NF, M, device)
            C, _ = generate_coded_aperture(G, NF, M, NF, device)
            C[C > 1] = 1
            C_initial = C.clone()

            # --- Initial Metrics ---
            initial_cond = get_cond(C)
            if M == 1:
                initial_cond = 1.0
            
            if R == 1:
                initial_density, _ = compute_density_irregular_sp(C)
            else:
                initial_density, _ = compute_density_regular_sp(C)

            sum_x_init = torch.sum(C_initial, dim=(1, 2))
            sum_y_init = torch.sum(C_initial, dim=(0, 2))
            sum_z_init = torch.sum(C_initial, dim=(0, 1))
            errors_init = torch.cat([sum_x_init - M, sum_y_init - M, sum_z_init - M])
            initial_rmse = torch.sqrt(torch.mean(errors_init.float()**2)).item()
            
            best_cond = initial_cond
            best_density = initial_density
            best_C = C.clone()
            
            if M > 1:
                cond_history = []
                density_history = []
                no_improvement_count = 0
                print(f'Starting Optimization M={M} | Initial Cond: {initial_cond:.2f} | Initial Density: {initial_density:.4f}')
                
                for k in range(max_iterations):
                    # --- Optimization loop ---
                    list_dist, diameters = unique_distance_irregular_sp(C)
                    X_coords = torch.nonzero(C > 0)

                    if len(list_dist) == 0 or len(X_coords) == 0:
                        print(">> No critical points found. Skipping iteration.")
                        continue
                    
                    positions = torch.nonzero(diameters == list_dist[0]).squeeze(1)
                    if len(positions) == 0:
                        continue
                        
                    val = torch.randint(0, len(positions), (1,)).item()
                    i, ii, j = X_coords[positions[val]]
                    
                    C_next = C.clone()
                    otherLayer = torch.randint(0, NF, (1,)).item()
                    
                    if torch.rand(1).item() < 0.5:
                        # XZ Plane swap
                        temp = C_next[i, :, j].clone()
                        C_next[i, :, j] = C_next[i, :, otherLayer]
                        C_next[i, :, otherLayer] = temp
                    else:
                        # YZ Plane swap
                        temp = C_next[:, ii, j].clone()
                        C_next[:, ii, j] = C_next[:, ii, otherLayer]
                        C_next[:, ii, otherLayer] = temp
                    
                    C_next[C_next > 1] = 1
                    
                    new_cond = get_cond(C_next)
                    if R == 1:
                        new_density, _ = compute_density_irregular_sp(C_next)
                    else:
                        new_density, _ = compute_density_regular_sp(C_next)

                    improvement_found = False
                    cond_history.append(new_cond)
                    density_history.append(new_density)

                    if optimization_metric == 'density':
                        if new_density > best_density:
                            improvement_ratio = (new_density - best_density) / best_density if best_density > 0 else float('inf')
                            if improvement_ratio > rel_tol:
                                best_density = new_density
                                best_cond = new_cond
                                improvement_found = True
                    else:
                        if new_cond < best_cond:
                             if (best_cond - new_cond) / best_cond > rel_tol:
                                best_cond = new_cond
                                best_density = new_density
                                improvement_found = True

                    if improvement_found:
                        best_C = C_next
                        C = C_next
                        no_improvement_count = 0
                    else:
                        no_improvement_count += 1
                    
                    if (k + 1) % 20 == 0:
                        print(f'Iter: {k+1} | Dens: {new_density:.4f} | Best Dens: {best_density:.4f} | Status: {no_improvement_count}/{patience}')

                    if k > burn_in and no_improvement_count >= patience:
                        print(f'>> Convergence reached at iteration {k+1}.')
                        break
                
                C_optimized = best_C

                # --- Plotting evolution graphs ---
                if plot_graphs and cond_history:
                    num_iters = len(cond_history)
                    iterations = range(1, num_iters + 1)

                    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
                    fig.suptitle(f'Optimization Evolution for M={M}', fontsize=16)

                    # Plot Condition Number
                    ax1.plot(iterations, cond_history, color='tab:blue', label='Condition Number')
                    min_cond_val = min(cond_history)
                    ax1.axhline(y=min_cond_val, color='r', linestyle='--', label=f'Min Cond: {min_cond_val:.2f}')
                    ax1.set_yscale('log')
                    ax1.set_title('Condition Number Evolution')
                    ax1.set_xlabel('Iterations')
                    ax1.set_ylabel('Condition Number (log scale)')
                    ax1.grid(True)
                    ax1.legend()

                    # Plot Density
                    ax2.plot(iterations, density_history, color='tab:orange', label='Density')
                    max_density_val = max(density_history)
                    ax2.axhline(y=max_density_val, color='g', linestyle='--', label=f'Max Density: {max_density_val:.4f}')
                    ax2.set_title('Density Evolution')
                    ax2.set_xlabel('Iterations')
                    ax2.set_ylabel('Density')
                    ax2.grid(True)
                    ax2.legend()

                    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
                    graph_filename = os.path.join(graph_export_path, f'Opt_Mux_{M}.png')
                    plt.savefig(graph_filename)
                    plt.close(fig)
                    print(f'>> Graph for M={M} saved to: {graph_filename}')

                # --- Uniformization ---
                print(">> Applying uniformization process...")
                C_uniformized = uniformize_matrix(C_optimized, M)
                final_cond = get_cond(C_uniformized)
                if R == 1:
                    density_after, _ = compute_density_irregular_sp(C_uniformized)
                else:
                    density_after, _ = compute_density_regular_sp(C_uniformized)
                
            else:
                C_optimized = C
                C_uniformized = C
                final_cond = best_cond
                density_after = best_density
                print("M=1 detected: Skipping optimization and uniformization.")

            # --- Saving results ---
            # Bajamos los tensores de GPU a CPU y los convertimos a NumPy para guardarlos en .mat
            folderName = "./PythonCode/Patrones/SSPv3"
            os.makedirs(folderName, exist_ok=True)
            filename = os.path.join(folderName, f'Bands_{NF}_Mux_{M}.mat')
            
            scipy.io.savemat(filename,
                             {'C_initial': C_initial.cpu().numpy(),
                              'C_optimized': C_optimized.cpu().numpy(),
                              'C_uniformized': C_uniformized.cpu().numpy()})

            print(f"-------------------- Summary M={M} --------------------")
            print(f"  Condition: Initial: {initial_cond:.2f} | Optimized: {best_cond:.2f} | Uniform: {final_cond:.4f}")
            print(f"  Density:   Initial: {initial_density:.4f} | Optimized: {best_density:.4f} | Uniform: {density_after:.4f}\n")

            # --- Log metrics to summary file ---
            data_row = (f"{M:<4} | {initial_cond:<14.2f} | {best_cond:<16.2f} | {final_cond:<14.4f} | "
                        f"{initial_density:<17.4f} | {best_density:<19.4f} | {density_after:<17.4f}\n")
            fileID.write(data_row)

if __name__ == '__main__':
    main()