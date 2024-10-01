/*
   IGraph library.
   Copyright (C) 2024  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph.h>
#include <stdio.h>
#include <time.h>

#include "bench.h"

void run_bench(const igraph_t *graph, const igraph_vector_t *weights,
               const char *name, int rep) {
    igraph_vector_int_t membership;
    igraph_vector_t vertex_weight;
    igraph_integer_t vcount = igraph_vcount(graph);
    igraph_integer_t ecount = igraph_ecount(graph);
    char msg[256], msg2[128];

    igraph_vector_int_init(&membership, vcount);
    igraph_vector_init(&vertex_weight, vcount);

    igraph_strength(graph, &vertex_weight, igraph_vss_all(), IGRAPH_ALL, true, weights);

    snprintf(msg2, sizeof(msg2) / sizeof(msg2[0]),
             "%s, vcount=%" IGRAPH_PRId ", ecount=%" IGRAPH_PRId ", %s, %dx",
             name, vcount, ecount, weights == NULL ? "unweighted" : "weighted",
             rep);

    snprintf(msg, sizeof(msg) / sizeof(msg[0]), "1 Louvain, %s", msg2);
    BENCH(msg, REPEAT(igraph_community_multilevel(graph, weights, 1.0, &membership, NULL, NULL), rep));

    snprintf(msg, sizeof(msg) / sizeof(msg[0]), "2 Leiden , %s", msg2);
    BENCH(msg, REPEAT(igraph_community_leiden(graph, weights, &vertex_weight, 1.0 / igraph_vector_sum(weights), 0.01, false, 1, &membership, NULL, NULL), rep));

    printf("\n");

    igraph_vector_destroy(&vertex_weight);
    igraph_vector_int_destroy(&membership);
}

// Function to run benchmarks for creating induced subgraphs
void run_induced_subgraph_bench(const igraph_t *graph, const char *name, int rep) {
    igraph_vector_int_t vs_vector;
    igraph_vs_t vs;
    igraph_t subgraph;
    char msg[256];
    igraph_integer_t vcount = igraph_vcount(graph);

    // Initialize the vertex selector vector
    igraph_vector_int_init(&vs_vector, 0);

    for (double ratio = 0.1; ratio <= 1.0; ratio += 0.1) {
        int vertices_to_select = (int)(vcount * ratio);

        // Resize the vector to the required number of vertices
        igraph_vector_int_resize(&vs_vector, vertices_to_select);

        // Set the values manually instead of using the deprecated function
        for (int i = 0; i < vertices_to_select; i++) {
            VECTOR(vs_vector)[i] = i;
        }

        printf("Creating vertex selector for vertices 0 to %d\n", vertices_to_select - 1);
        igraph_vs_vector(&vs, &vs_vector);

        snprintf(msg, sizeof(msg) / sizeof(msg[0]), "Induced subgraph creation for %s, ratio %.2f, vcount=%d", name, ratio, vertices_to_select);

        // Measure time taken to create the induced subgraph
        BENCH(msg, {
            igraph_induced_subgraph(graph, &subgraph, vs, IGRAPH_SUBGRAPH_AUTO);
        });

        igraph_vs_destroy(&vs);
        igraph_destroy(&subgraph);
    }

    igraph_vector_int_destroy(&vs_vector);
}

// Function to generate random weights for the graph
void rand_weights(const igraph_t *graph, igraph_vector_t *weights) {
    igraph_integer_t ecount = igraph_ecount(graph);
    igraph_vector_resize(weights, ecount);
    for (igraph_integer_t i = 0; i < ecount; i++) {
        VECTOR(*weights)[i] = RNG_UNIF01();
    }
}

// Main function
int main(void) {
    igraph_t graph;
    igraph_vector_t weights;

    igraph_rng_seed(igraph_rng_default(), 137);
    BENCH_INIT();

    igraph_vector_init(&weights, 0);

    // Benchmark for Erdos-Renyi graphs
    igraph_erdos_renyi_game_gnm(&graph, 100, 500, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 1000);
    run_induced_subgraph_bench(&graph, "G(n,m)", 1000);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 1000, 5000, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 100);
    run_induced_subgraph_bench(&graph, "G(n,m)", 100);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 1000, 50000, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 10);
    run_induced_subgraph_bench(&graph, "G(n,m)", 10);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 10000, 50000, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 10);
    run_induced_subgraph_bench(&graph, "G(n,m)", 10);
    igraph_destroy(&graph);

    igraph_erdos_renyi_game_gnm(&graph, 100000, 500000, false, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "G(n,m)", 1);
    run_induced_subgraph_bench(&graph, "G(n,m)", 1);
    igraph_destroy(&graph);

    // Benchmark for Forest Fire graphs
    igraph_forest_fire_game(&graph, 1000, 0.2, 1, 2, false);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "forest fire", 100);
    run_induced_subgraph_bench(&graph, "forest fire", 100);
    igraph_destroy(&graph);

    // Benchmark for Barabasi-Albert graphs
    igraph_barabasi_game(&graph, 1000, 1, 5, NULL, true, 0, false, IGRAPH_BARABASI_PSUMTREE, NULL);
    rand_weights(&graph, &weights);
    run_bench(&graph, &weights, "PA", 100);
    run_induced_subgraph_bench(&graph, "PA", 100);
    igraph_destroy(&graph);

    igraph_vector_destroy(&weights);

    return 0;
}
