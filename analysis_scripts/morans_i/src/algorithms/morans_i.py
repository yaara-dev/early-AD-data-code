#!/usr/bin/env python3
"""
Moran's I Spatial Autocorrelation Algorithm

This module contains the core algorithm for calculating Moran's I spatial autocorrelation
statistic for gene expression data. The algorithm is independent of the spatial unit
(cells or grid) and focuses purely on the mathematical computation.
"""

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Any, Literal
from scipy.spatial.distance import cdist


def calculate_spatial_weights_matrix(coordinates: np.ndarray,
                                     weight_scheme: Literal["inverse_square", "inverse", "binary_4", "binary_6"] = "inverse_square",
                                     grid_ids: np.ndarray = None,
                                     max_x_cells: int = None,
                                     max_y_cells: int = None) -> np.ndarray:
    """
    Calculate spatial weights matrix using various schemes.

    Parameters:
        coordinates: Array of shape (n_units, 3) with X, Y, Z coordinates
        weight_scheme: Weight scheme to use:
            - 'inverse_square': 1/d^2 weights (default)
            - 'inverse': 1/d weights
            - 'binary_4': Binary adjacency in XY plane only
            - 'binary_6': Binary adjacency in full 3D
        grid_ids: Grid cell IDs (required for binary schemes)
        max_x_cells: Number of X cells (required for binary schemes)
        max_y_cells: Number of Y cells (required for binary schemes)

    Returns:
        Spatial weights matrix W of shape (n_units, n_units)

    Example:
        >>> coords = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
        >>> W = calculate_spatial_weights_matrix(coords, weight_scheme="inverse")
        >>> W.shape
        (3, 3)
    """
    # Handle binary adjacency schemes
    if weight_scheme in ["binary_4", "binary_6"]:
        if grid_ids is None or max_x_cells is None or max_y_cells is None:
            raise ValueError(f"binary schemes require grid_ids, max_x_cells, and max_y_cells parameters")
        
        from ..utils.spatial import calculate_binary_weights_from_grid_ids
        W = calculate_binary_weights_from_grid_ids(grid_ids, max_x_cells, max_y_cells, weight_scheme)
        
        # Debug: Check binary weights matrix
        if np.isnan(W).any():
            print(f"      🔍 DEBUG: NaN detected in binary weights matrix:")
            print(f"         grid_ids shape: {grid_ids.shape}")
            print(f"         max_x_cells: {max_x_cells}, max_y_cells: {max_y_cells}")
            print(f"         W shape: {W.shape}")
            print(f"         W has NaN: {np.isnan(W).any()}")
        
        return W
    
    # Handle distance-based schemes (original logic)
    # Calculate Euclidean distance matrix
    distance_matrix = cdist(coordinates, coordinates, metric='euclidean')

    # Replace zeros with a small value to avoid division by zero
    distance_matrix[distance_matrix == 0] = 1e-10

    # Calculate weights
    if weight_scheme == "inverse_square":
        W = 1.0 / (distance_matrix ** 2)
    elif weight_scheme == "inverse":
        W = 1.0 / distance_matrix
    else:
        raise ValueError(f"Unsupported weight_scheme: {weight_scheme}. Use 'inverse_square', 'inverse', 'binary_4', or 'binary_6'.")

    # Set diagonal to zero (no self-influence)
    np.fill_diagonal(W, 0)

    # Debug: Check for NaN values in weights matrix
    if np.isnan(W).any():
        print(f"      🔍 DEBUG: NaN detected in spatial weights matrix:")
        print(f"         coordinates shape: {coordinates.shape}")
        print(f"         distance_matrix has NaN: {np.isnan(distance_matrix).any()}")
        print(f"         W shape: {W.shape}")
        print(f"         Sample coordinates: {coordinates[:3] if len(coordinates) >= 3 else coordinates}")

    return W


def calculate_morans_i(expression_data: np.ndarray, weights_matrix: np.ndarray) -> float:
    """
    Calculate Moran's I statistic for a single gene.

    Parameters:
        expression_data: 1D array of expression values for spatial units
        weights_matrix: Spatial weights matrix W

    Returns:
        Moran's I statistic

    Raises:
        ValueError: If input arrays have incompatible shapes

    Example:
        >>> expression = np.array([1, 2, 3, 4, 5])
        >>> W = np.ones((5, 5)) - np.eye(5)  # Simple weights
        >>> morans_i = calculate_morans_i(expression, W)
        >>> isinstance(morans_i, float)
        True
    """
    if len(expression_data) != weights_matrix.shape[0]:
        raise ValueError("Expression data length must match weights matrix dimensions")

    n = len(expression_data)
    mean_expression = np.mean(expression_data)
    deviation = expression_data - mean_expression

    # Calculate numerator: sum of weighted cross-products
    numerator = np.sum(weights_matrix * np.outer(deviation, deviation))

    # Calculate denominator: sum of squared deviations
    denominator = np.sum(deviation ** 2)

    # Handle case where all values are the same (zero variance)
    if denominator == 0:
        return 0.0

    # Calculate Moran's I
    W_sum = np.sum(weights_matrix)
    if W_sum == 0:
        return 0.0

    # Debug: Check for potential NaN sources
    if np.isnan(numerator) or np.isnan(denominator) or np.isnan(W_sum):
        print(f"      🔍 DEBUG: NaN detected in calculation:")
        print(f"         expression_data: {expression_data}")
        print(f"         mean_expression: {mean_expression}")
        print(f"         deviation: {deviation}")
        print(f"         numerator: {numerator}")
        print(f"         denominator: {denominator}")
        print(f"         W_sum: {W_sum}")
        print(f"         weights_matrix has NaN: {np.isnan(weights_matrix).any()}")
        print(f"         weights_matrix shape: {weights_matrix.shape}")

    morans_i = (n / W_sum) * (numerator / denominator)

    # Check if result is NaN
    if np.isnan(morans_i):
        print(f"      🔍 DEBUG: Moran's I result is NaN: {morans_i}")
        print(f"         n: {n}, W_sum: {W_sum}, numerator: {numerator}, denominator: {denominator}")

    return morans_i


def build_expression_matrix(data: pd.DataFrame, spatial_units: List[str], 
                          genes: List[str], unit_col: str = 'Cell', 
                          gene_col: str = 'Gene') -> np.ndarray:
    """
    Build expression matrix from spatial data.

    Parameters:
        data: DataFrame with spatial data
        spatial_units: List of spatial unit IDs (cells or grid IDs)
        genes: List of gene names
        unit_col: Column name for spatial units
        gene_col: Column name for genes

    Returns:
        Expression matrix of shape (n_units, n_genes)

    Example:
        >>> data = pd.DataFrame({
        ...     'Cell': ['C1', 'C1', 'C2', 'C2'],
        ...     'Gene': ['ACTB', 'GAPDH', 'ACTB', 'GAPDH']
        ... })
        >>> units = ['C1', 'C2']
        >>> genes = ['ACTB', 'GAPDH']
        >>> matrix = build_expression_matrix(data, units, genes)
        >>> matrix.shape
        (2, 2)
    """
    expression_matrix = np.zeros((len(spatial_units), len(genes)))

    # Create mapping dictionaries for faster lookup
    unit_to_idx = {unit: idx for idx, unit in enumerate(spatial_units)}
    gene_to_idx = {gene: idx for idx, gene in enumerate(genes)}

    # Count gene expressions for each spatial unit
    for _, row in data.iterrows():
        unit_id = row[unit_col]
        gene_name = row[gene_col]

        if unit_id in unit_to_idx and gene_name in gene_to_idx:
            unit_idx = unit_to_idx[unit_id]
            gene_idx = gene_to_idx[gene_name]
            expression_matrix[unit_idx, gene_idx] += 1

    return expression_matrix


def calculate_morans_i_for_genes(expression_matrix: np.ndarray, 
                               weights_matrix: np.ndarray,
                               gene_names: List[str]) -> Dict[str, float]:
    """
    Calculate Moran's I for multiple genes.

    Parameters:
        expression_matrix: Matrix of shape (n_units, n_genes)
        weights_matrix: Spatial weights matrix
        gene_names: List of gene names corresponding to matrix columns

    Returns:
        Dictionary mapping gene names to Moran's I values

    Example:
        >>> expr_matrix = np.random.rand(10, 5)
        >>> W = np.random.rand(10, 10)
        >>> genes = ['ACTB', 'GAPDH', 'TUBB3', 'MAP2', 'SNAP25']
        >>> results = calculate_morans_i_for_genes(expr_matrix, W, genes)
        >>> len(results) == 5
        True
    """
    results = {}

    for gene_idx, gene_name in enumerate(gene_names):
        expression_data = expression_matrix[:, gene_idx]
        
        # Check if gene has any expression data
        if np.sum(expression_data) == 0:
            # Gene has no expression data - return NaN
            results[gene_name] = np.nan
        else:
            # Gene has expression data - calculate Moran's I
            morans_i = calculate_morans_i(expression_data, weights_matrix)
            results[gene_name] = morans_i

    return results


def process_morans_i_analysis(coordinates: np.ndarray, 
                            expression_matrix: np.ndarray,
                            gene_names: List[str],
                            weight_scheme: Literal["inverse_square", "inverse", "binary_4", "binary_6"] = "inverse_square",
                            grid_ids: np.ndarray = None,
                            max_x_cells: int = None,
                            max_y_cells: int = None) -> Dict[str, float]:
    """
    Complete Moran's I analysis pipeline.

    Parameters:
        coordinates: Array of shape (n_units, 3) with spatial coordinates
        expression_matrix: Matrix of shape (n_units, n_genes)
        gene_names: List of gene names
        weight_scheme: Weight scheme to use:
            - 'inverse_square': 1/d^2 weights (default)
            - 'inverse': 1/d weights
            - 'binary_4': Binary adjacency in XY plane only
            - 'binary_6': Binary adjacency in full 3D
        grid_ids: Grid cell IDs (required for binary schemes)
        max_x_cells: Number of X cells (required for binary schemes)
        max_y_cells: Number of Y cells (required for binary schemes)

    Returns:
        Dictionary mapping gene names to Moran's I values

    Example:
        >>> coords = np.random.rand(20, 3)
        >>> expr_matrix = np.random.randint(0, 10, (20, 3))
        >>> genes = ['ACTB', 'GAPDH', 'TUBB3']
        >>> results = process_morans_i_analysis(coords, expr_matrix, genes, weight_scheme="inverse")
        >>> len(results) == 3
        True
    """
    # Calculate spatial weights matrix
    weights_matrix = calculate_spatial_weights_matrix(
        coordinates, 
        weight_scheme=weight_scheme,
        grid_ids=grid_ids,
        max_x_cells=max_x_cells,
        max_y_cells=max_y_cells
    )

    # Calculate Moran's I for all genes
    results = calculate_morans_i_for_genes(expression_matrix, weights_matrix, gene_names)

    return results

