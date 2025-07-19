#include <iostream>
#include <cstdlib>
#include <string>
#include <cstdint>
#include <cstddef>
#include <omp.h>


//************************************ GENERATE AND FILTER MATRICES ************************************//

uint8_t ***generate_matrices(int n, size_t *out_count) {
    size_t total = 1ULL << (n * n);  // total = 2^(n^2)
    *out_count = total;


    uint8_t ***matrices = new uint8_t**[total];

#pragma omp parallel for schedule(static)
    for (size_t k = 0; k < total; k++) {

        uint8_t **matrix = new uint8_t*[n];
        for (int i = 0; i < n; i++) {
            matrix[i] = new uint8_t[n];
        }


        for (int i = 0; i < n * n; i++) {
            int row = i / n;
            int col = i % n;
            matrix[row][col] = (k >> (n * n - 1 - i)) & 1;
        }

        matrices[k] = matrix;
    }

    return matrices;
}


bool is_weakly_connected(uint8_t **matrix, int n) {
    int *visited = new int[n];
    int *queue = new int[n];

    for (int i = 0; i < n; i++) visited[i] = 0;

    int front = 0, rear = 0;
    queue[rear++] = 0;
    visited[0] = 1;

    while (front < rear) {
        int u = queue[front++];
        for (int v = 0; v < n; v++) {
            if (!visited[v] && (matrix[u][v] || matrix[v][u])) {
                visited[v] = 1;
                queue[rear++] = v;
            }
        }
    }

    bool connected = true;
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            connected = false;
            break;
        }
    }

    delete[] visited;
    delete[] queue;

    return connected;
}


uint8_t ***filter_weakly_connected(uint8_t ***matrices, size_t count, int n, size_t *filtered_count) {
    uint8_t ***filtered = new uint8_t[count];
    size_t *local_counts = new size_t[omp_get_max_threads()];
    uint8_t ***local_buffers = new uint8_t[omp_get_max_threads()];

#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        size_t local_cap = count / omp_get_num_threads() + 1;
        uint8_t ***local = new uint8_t[local_cap];
        size_t local_index = 0;

    }
}





//************************************    APPLY MAJORITY RULE     ************************************//


int main(int argc, char *argv[]) {

    size_t qtde_matrices;
    uint8_t ***matrices = generate_matrices(3, &qtde_matrices);

    std::cout << qtde_matrices << std::endl;

    return EXIT_SUCCESS;

}


