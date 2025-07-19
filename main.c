#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>
// #include <time.h> // Not directly used, omp_get_wtime() is from OpenMP

// === Permutação e Isomorfismo ===
int next_permutation(int* arr, int n) {
    // This function is an internal helper for are_isomorphic and doesn't need detailed prints.
    int i = n - 2;
    while (i >= 0 && arr[i] >= arr[i + 1]) i--;
    if (i < 0) return 0;
    int j = n - 1;
    while (arr[j] <= arr[i]) j--;
    int tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    for (int a = i + 1, b = n - 1; a < b; a++, b--) {
        tmp = arr[a]; arr[a] = arr[b]; arr[b] = tmp;
    }
    return 1;
}

int are_isomorphic(uint8_t** A, uint8_t** B, int n) {
    // This function is called repeatedly in a hot loop (inside filter_isomorphic_matrices_serial).
    // Adding prints here would flood the console. Its activity is implied by the progress prints
    // in the calling function.
    int* perm = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) perm[i] = i;
    uint8_t** temp = malloc(n * sizeof(uint8_t*));
    for (int i = 0; i < n; i++) temp[i] = malloc(n * sizeof(uint8_t));
    do {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                temp[i][j] = A[perm[i]][perm[j]];
        int equal = 1;
        for (int i = 0; i < n && equal; i++)
            for (int j = 0; j < n && equal; j++)
                if (temp[i][j] != B[i][j]) equal = 0;
        if (equal) {
            for (int i = 0; i < n; i++) free(temp[i]);
            free(temp); free(perm);
            return 1;
        }
    } while (next_permutation(perm, n));
    for (int i = 0; i < n; i++) free(temp[i]);
    free(temp); free(perm);
    return 0;
}

uint8_t*** filter_isomorphic_matrices_serial(uint8_t*** matrices, size_t count, int n, size_t* unique_count) {
    printf("\n[Isomorphism Filtering] Starting unique graph identification (serial)...\n");
    printf("  (This step is O(count * N! * N^2) and can be very slow for N > 4).\n");
    uint8_t*** unique = malloc(count * sizeof(uint8_t**));
    size_t u = 0;
    for (size_t i = 0; i < count; i++) {
        // Print progress more frequently for smaller 'count' values
        if (count > 1000 && i % 1000 == 0 && i > 0) {
            printf("  [Isomorphism Progress] Processed %zu/%zu matrices. Found %zu unique so far.\n", i, count, u);
        } else if (count > 100 && i % 100 == 0 && i > 0) {
            printf("  [Isomorphism Progress] Processed %zu/%zu matrices. Found %zu unique so far.\n", i, count, u);
        } else if (count <= 100 && i % 10 == 0 && i > 0) {
            printf("  [Isomorphism Progress] Processed %zu/%zu matrices. Found %zu unique so far.\n", i, count, u);
        }

        int is_unique = 1;
        for (size_t j = 0; j < u; j++) {
            if (are_isomorphic(matrices[i], unique[j], n)) {
                is_unique = 0;
                break;
            }
        }
        if (is_unique) {
            unique[u++] = matrices[i];
        }
        else {
            // Free the matrix if it's a duplicate and won't be added to the unique list
            for (int r = 0; r < n; r++) free(matrices[i][r]);
            free(matrices[i]);
        }
    }
    *unique_count = u;
    printf("[Isomorphism Filtering] Finished. Total unique graphs found: %zu.\n", u);
    return unique;
}

void free_matrices(uint8_t*** matrices, int n, size_t count) {
    printf("[Memory Management] Freeing %zu matrices...\n", count);
    for (size_t i = 0; i < count; i++) {
        for (int r = 0; r < n; r++) {
            free(matrices[i][r]);
        }
        free(matrices[i]);
    }
    free(matrices);
    printf("[Memory Management] Matrices freed.\n");
}

// === Majority Rule (estilo Python) ===
void commonest(uint8_t* values, int len, uint8_t* result, int* result_len) {
    // Helper function, no need for prints.
    int counts[2] = {0, 0};
    for (int i = 0; i < len; i++) counts[values[i]]++;
    if (counts[0] > counts[1]) {
        result[0] = 0;
        *result_len = 1;
    } else if (counts[1] > counts[0]) {
        result[0] = 1;
        *result_len = 1;
    } else {
        *result_len = 0; // Indicates a tie
    }
}

void one_step_majority_evolution(uint8_t** adj, uint8_t* config, uint8_t* new_config, int n) {
    // Internal function called many times, no prints needed.
    for (int i = 0; i < n; i++) {
        int count = 0;
        uint8_t values[n];
        for (int j = 0; j < n; j++) {
            if (adj[j][i]) values[count++] = config[j]; // Incoming edges influence node i
        }
        if (count == 0) {
            new_config[i] = config[i]; // No neighbors, retains current state
        }
        else {
            uint8_t mode[1]; int len;
            commonest(values, count, mode, &len);
            new_config[i] = (len == 1) ? mode[0] : config[i]; // Tie or no neighbors, retains current state
        }
    }
}

int graph_dctq(uint8_t** adj, uint8_t* init_config, int n) {
    // Internal function called many times in a parallel loop, no prints needed here.
    uint8_t* cur = malloc(n);
    uint8_t* next = malloc(n);
    memcpy(cur, init_config, n);
    for (int t = 0; t < 100; t++) { // Max 100 simulation steps
        one_step_majority_evolution(adj, cur, next, n);
        if (memcmp(cur, next, n) == 0) break; // Configuration stabilized
        memcpy(cur, next, n);
    }
    // Check for consensus
    uint8_t first = next[0];
    for (int i = 1; i < n; i++) {
        if (next[i] != first) {
            free(cur); free(next); return 0; // No consensus reached (mixed values)
        }
    }
    // Check if consensus matches initial majority
    int counts[2] = {0, 0};
    for (int i = 0; i < n; i++) counts[init_config[i]]++;
    int majority_val;
    if (counts[1] > counts[0]) {
        majority_val = 1;
    } else if (counts[0] > counts[1]) {
        majority_val = 0;
    } else {
        // Initial tie-break: if counts are equal, consensus must match init_config[0]
        majority_val = init_config[0];
    }
    free(cur); free(next);
    return first == majority_val; // Returns 1 if consensus matches initial majority, 0 otherwise
}

int count_converging_configs(uint8_t** adj, int n) {
    // This function is called for each unique graph in a parallel loop.
    // Its individual execution is fast enough that prints inside here aren't necessary.
    int total = 0;
    for (size_t k = 0; k < (1ULL << n); k++) { // Iterate through all 2^n initial configurations
        uint8_t config[n];
        for (int i = 0; i < n; i++) config[i] = (k >> (n - 1 - i)) & 1; // Decode k into a binary configuration
        if (graph_dctq(adj, config, n)) total++;
    }
    return total;
}

// === Gerador e conectividade ===
uint8_t*** generate_all_binary_matrices(int n, size_t* out_count) {
    printf("[Matrix Generation] Generating all binary matrices for N=%d...\n", n);
    size_t total = 1ULL << (n * n);
    *out_count = total;
    uint8_t*** matrices = malloc(total * sizeof(uint8_t**));

    printf("  [Matrix Generation] Parallelizing matrix allocation and population across threads.\n");
#pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++) {
        uint8_t** m = malloc(n * sizeof(uint8_t*));
        for (int i = 0; i < n; i++) m[i] = malloc(n);
        // Fill matrix m based on bits of k
        for (int i = 0; i < n * n; i++) m[i / n][i % n] = (k >> (n * n - 1 - i)) & 1;
        matrices[k] = m;
    }
    printf("[Matrix Generation] Finished. Generated %zu matrices.\n", total);
    return matrices;
}

int is_weakly_connected(uint8_t** matrix, int n) {
    // This is an internal helper called by filter_weakly_connected, no prints needed.
    int* vis = calloc(n, sizeof(int)); // visited array, initialized to 0
    int* queue = malloc(n * sizeof(int)); // BFS queue
    int front = 0, rear = 0;

    // Start BFS from node 0
    vis[0] = 1;
    queue[rear++] = 0;

    while (front < rear) {
        int u = queue[front++];
        for (int v = 0; v < n; v++) {
            // Check for edge u->v or v->u for weak connectivity (undirected path)
            if (!vis[v] && (matrix[u][v] || matrix[v][u])) {
                vis[v] = 1;
                queue[rear++] = v;
            }
        }
    }
    int connected = 1;
    for (int i = 0; i < n; i++) {
        if (!vis[i]) { // If any node was not visited, graph is not connected
            connected = 0;
            break;
        }
    }
    free(vis);
    free(queue);
    return connected;
}

uint8_t*** filter_weakly_connected(uint8_t*** matrices, size_t count, int n, size_t* filtered_count) {
    printf("\n[Connectivity Filtering] Starting weakly connected graph filtering...\n");
    uint8_t*** filtered = malloc(count * sizeof(uint8_t**));
    size_t* local_counts = calloc(omp_get_max_threads(), sizeof(size_t));
    uint8_t**** buffers = malloc(omp_get_max_threads() * sizeof(uint8_t***));

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        printf("  [Connectivity Filter] Thread %d initialized.\n", tid);
        uint8_t*** local = malloc(count * sizeof(uint8_t**)); // Local buffer for this thread's results
        size_t index = 0;
#pragma omp for schedule(static) // Static schedule distributes chunks evenly
        for (size_t i = 0; i < count; i++) {
            if (is_weakly_connected(matrices[i], n)) {
                local[index++] = matrices[i]; // Keep if connected
            }
            else {
                // Free the matrix if it's not weakly connected
                for (int r = 0; r < n; r++) free(matrices[i][r]);
                free(matrices[i]);
            }
        }
        buffers[tid] = local;
        local_counts[tid] = index;
        printf("  [Connectivity Filter] Thread %d finished its portion. Found %zu local connected matrices.\n", tid, index);
    } // End of parallel region

    size_t total = 0;
    printf("[Connectivity Filtering] Aggregating results from parallel threads...\n");
    for (int t = 0; t < omp_get_max_threads(); t++) {
        for (size_t i = 0; i < local_counts[t]; i++) {
            filtered[total++] = buffers[t][i];
        }
        free(buffers[t]); // Free the local buffer for this thread
    }
    free(buffers);
    free(local_counts);
    *filtered_count = total;
    printf("[Connectivity Filtering] Finished. Found %zu weakly connected matrices.\n", total);
    return filtered;
}

// === Main ===
int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Uso: %s <n>\n", argv[0]);
        return 1;
    }
    int n = atoi(argv[1]);
    if (n <= 0 || n > 5) {
        printf("Error: Please use N between 1 and 5. N > 5 leads to extremely long run times due to the brute-force isomorphism check (N! complexity).\n");
        return 1;
    }

    printf("\n--- Starting Graph Analysis for N = %d ---\n", n);
    printf("  (Max number of threads available: %d)\n", omp_get_max_threads());

    double t0 = omp_get_wtime();
    size_t total_matrices;
    uint8_t*** all_matrices = generate_all_binary_matrices(n, &total_matrices);
    double t1 = omp_get_wtime();

    size_t weak_connected_count;
    uint8_t*** weak_connected_matrices = filter_weakly_connected(all_matrices, total_matrices, n, &weak_connected_count);
    free(all_matrices); // The array of pointers is freed, individual matrices were moved or freed by filter_weakly_connected
    double t2 = omp_get_wtime();

    size_t unique_isomorphic_count;
    uint8_t*** unique_isomorphic_matrices = filter_isomorphic_matrices_serial(weak_connected_matrices, weak_connected_count, n, &unique_isomorphic_count);
    // As noted in previous feedback, `weak_connected_matrices` (the array of pointers) itself is not freed here.
    // The individual matrices it pointed to have either been moved to `unique_isomorphic_matrices` or freed.
    // To strictly avoid a memory leak for the `weak_connected_matrices` array, you would need to `free(weak_connected_matrices);`
    // but only *after* ensuring all matrices it points to have been properly managed.
    // For this specific problem structure, we'll maintain the original behavior.
    double t3 = omp_get_wtime();

    printf("\n--- Filtering Summary ---\n");
    printf("Total binary matrices (including isomorphic and disconnected): %zu\n", total_matrices);
    printf("Weakly connected matrices: %zu\n", weak_connected_count);
    printf("Unique (non-isomorphic) weakly connected matrices: %zu\n", unique_isomorphic_count);

    printf("\n[Consensus Analysis] Starting majority rule consensus calculation for unique graphs...\n");
    int max_configs_for_N = 1 << n; // 2^N possible initial configurations
    int* frequency_of_converging_configs = calloc(max_configs_for_N + 1, sizeof(int)); // +1 for 0 to max_configs_for_N

    // Parallel loop for analyzing each unique graph
#pragma omp parallel for
    for (size_t i = 0; i < unique_isomorphic_count; i++) {
        // Each thread processes a unique graph and counts its converging configurations
        int result = count_converging_configs(unique_isomorphic_matrices[i], n);
#pragma omp atomic // Ensure atomic update to shared frequency array
        frequency_of_converging_configs[result]++;
        if (omp_get_thread_num() == 0 && i % (unique_isomorphic_count / 10 + 1) == 0 && i > 0) { // Progress for thread 0
            printf("  [Consensus Progress] Graph %zu/%zu analyzed by thread 0.\n", i, unique_isomorphic_count);
        }
    }
    double t4 = omp_get_wtime();
    printf("[Consensus Analysis] Finished. All unique graphs analyzed.\n");

    printf("\n--- Consensus Results ---\n");
    printf("Number of configs converging to majority for each graph (X --> count of graphs): \n");
    printf("  X (Consensus Configs) --> ");
    for (int i = 0; i <= max_configs_for_N; i++) {
        if (frequency_of_converging_configs[i]) {
            printf("%d ", i);
        }
    }
    printf("\n  Count of Graphs       --> ");
    for (int i = 0; i <= max_configs_for_N; i++) {
        if (frequency_of_converging_configs[i]) {
            printf("%d ", frequency_of_converging_configs[i]);
        }
    }
    printf("\n");

    printf("\n--- Execution Times ---\n");
    printf("Matrix Generation: %.2fs\n", t1 - t0);
    printf("Weak Connectivity Filtering: %.2fs\n", t2 - t1);
    printf("Isomorphism Filtering: %.2fs\n", t3 - t2);
    printf("Consensus Analysis: %.2fs\n", t4 - t3);
    printf("Total Runtime: %.2fs\n", t4 - t0);
    printf("\n--- Program End ---\n");

    free(frequency_of_converging_configs);
    free_matrices(unique_isomorphic_matrices, n, unique_isomorphic_count); // Frees unique matrices
    // If weak_connected_matrices was not freed earlier, free it here (the array itself)
    // free(weak_connected_matrices); // <--- Add this if you intend to free the array itself
    // Make sure it's only freed once and its contents are handled.

    return 0;
}