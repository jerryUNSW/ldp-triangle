# Triangle Counting Algorithm under Edge Local Differential Privacy

## Overview

This project implements a triangle counting algorithm under edge LDP.

## Usage

To run the triangle counting algorithm, use the following script format:

```bash
./counting <epsilon> <dataset> <num_rounds> <vertex_ratio> <algo_switch>
```

Where:
- `<epsilon>`: Privacy budget (controls the trade-off between data utility and privacy).
- `<dataset>`: Path to the graph dataset in edge list format.
- `<num_rounds>`: Number of rounds to run the algorithm.
- `<vertex_ratio>`: Ratio for vertex sampling (1.0 indicates no sampling).
- `<algo_switch>`: An integer to switch between different algorithms (1 indicates the PC algorithm. 2 indicates the VC algorithm). 

## Input Data Format

The input data for the triangle counting algorithm consists of:

1. **First line**: Two integers `n` and `m`, where:
   - `n` is the number of vertices.
   - `m` is the number of edges.

2. **Subsequent lines**: Each line contains two integers representing an edge between two vertices.

### Example Input:

```
829 7138
3 6
3 10
3 14
3 17
3 19
3 20
3 24
3 25
3 26
3 28
3 29
3 30
```

- **First line**: `829 7138` indicates 829 vertices and 7138 edges.
- **Subsequent lines**: Each pair (e.g., `3 6`) represents an edge between vertex 3 and vertex 6.