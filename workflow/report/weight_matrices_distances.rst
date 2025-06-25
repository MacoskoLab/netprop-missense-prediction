Weight Matrix Distance Comparison
===================================

This visualization shows the pairwise distances between weight matrices using three different distance metrics:

**Distance Metrics:**

- **Euclidean Distance**: Measures the straight-line distance between two matrices when flattened into vectors
- **Wasserstein Distance**: Also known as Earth Mover's Distance, measures the minimum cost to transform one distribution into another
- **Energy Distance**: A metric that measures the distance between two probability distributions based on energy statistics

**Matrix Comparisons:**

- **Control vs. Experimental**: Comparison between unperturbed and perturbed real weight matrices
- **Predicted vs. Experimental**: Comparison between predicted perturbed weights and actual perturbed weights  
- **Predicted vs. Control**: Comparison between predicted perturbed weights and unperturbed control weights

The plots show the distance values for each pairwise comparison, with larger values indicating greater differences between the weight matrices. This analysis helps evaluate:

1. How much the perturbation changed the network weights (Control vs. Experimental)
2. How well the model predicted the perturbation effects (Predicted vs. Experimental)  
3. Whether the predictions properly account for the baseline differences (Predicted vs. Control)

Generated from comparison data: ``{{ snakemake.input.distances }}``
