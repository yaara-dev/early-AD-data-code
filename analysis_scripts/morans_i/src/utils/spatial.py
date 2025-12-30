#!/usr/bin/env python3
"""
Spatial Utilities

This module contains utilities for spatial data processing including coordinate
conversion, grid creation, and data preprocessing for spatial analysis.
"""

import numpy as np
import pandas as pd
from typing import Tuple, List, Dict, Any


def convert_coordinates_to_microns(data: pd.DataFrame, 
                                 x_col: str = 'global_x', 
                                 y_col: str = 'global_y', 
                                 z_col: str = 'z',
                                 x_factor: float = 0.17,
                                 y_factor: float = 0.17,
                                 z_factor: float = 0.4) -> pd.DataFrame:
    """
    Convert pixel coordinates to microns.
    
    Parameters:
        data: DataFrame with coordinate columns
        x_col: Name of X coordinate column
        y_col: Name of Y coordinate column  
        z_col: Name of Z coordinate column
        x_factor: Conversion factor for X (microns per pixel)
        y_factor: Conversion factor for Y (microns per pixel)
        z_factor: Conversion factor for Z (microns per pixel)
        
    Returns:
        DataFrame with converted coordinates in microns
        
    Example:
        >>> data = pd.DataFrame({'global_x': [10, 20], 'global_y': [15, 25], 'z': [5, 10]})
        >>> converted = convert_coordinates_to_microns(data)
        >>> converted['X'].iloc[0] == 1.7  # 10 * 0.17
        True
    """
    result = data.copy()
    result['X'] = data[x_col] * x_factor
    result['Y'] = data[y_col] * y_factor  
    result['Z'] = data[z_col] * z_factor
    
    return result


def apply_gene_name_mapping(data: pd.DataFrame, 
                          gene_col: str = 'Gene',
                          mapping: Dict[str, str] = None) -> pd.DataFrame:
    """
    Apply gene name standardization mapping.
    
    Parameters:
        data: DataFrame with gene column
        gene_col: Name of gene column
        mapping: Dictionary mapping old names to new names
        
    Returns:
        DataFrame with standardized gene names
        
    Example:
        >>> data = pd.DataFrame({'Gene': ['Irg1', 'ACTB', 'HIF-1?']})
        >>> mapping = {'Irg1': 'ACOD1', 'HIF-1?': 'HIF-1α'}
        >>> result = apply_gene_name_mapping(data, mapping=mapping)
        >>> result['Gene'].iloc[0] == 'ACOD1'
        True
    """
    if mapping is None:
        mapping = {
            'Irg1': 'ACOD1',
            'CcrlG': 'Ccrl2',
            'HIF-1?': 'HIF-1α',
            'Idb3': 'Id3',
            'iNOS': 'NOS1'
        }
    
    result = data.copy()
    result[gene_col] = result[gene_col].replace(mapping)
    
    return result


def filter_unassigned_cells(data: pd.DataFrame, 
                          cell_col: str = 'Cell',
                          suffix: str = '_0') -> pd.DataFrame:
    """
    Filter out unassigned cells (those ending with specified suffix).
    
    Parameters:
        data: DataFrame with cell column
        cell_col: Name of cell column
        suffix: Suffix indicating unassigned cells
        
    Returns:
        DataFrame with unassigned cells removed
        
    Example:
        >>> data = pd.DataFrame({'Cell': ['C1', 'C2_0', 'C3', 'C4_0']})
        >>> filtered = filter_unassigned_cells(data)
        >>> len(filtered) == 2
        True
    """
    mask = ~data[cell_col].astype(str).str.endswith(suffix)
    return data[mask]


def create_grid(x_min: float, x_max: float, y_min: float, y_max: float, 
               grid_size: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Create grid edges and centers for spatial gridding.
    
    Parameters:
        x_min: Minimum X coordinate
        x_max: Maximum X coordinate
        y_min: Minimum Y coordinate
        y_max: Maximum Y coordinate
        grid_size: Size of grid cells
        
    Returns:
        Tuple of (x_edges, y_edges, x_centers, y_centers)
        
    Example:
        >>> x_edges, y_edges, x_centers, y_centers = create_grid(0, 10, 0, 10, 5)
        >>> len(x_centers) == 2
        True
    """
    # Create grid edges
    x_edges = np.arange(x_min, x_max + grid_size, grid_size)
    y_edges = np.arange(y_min, y_max + grid_size, grid_size)
    
    # Calculate grid centers
    x_centers = x_edges[:-1] + grid_size / 2
    y_centers = y_edges[:-1] + grid_size / 2
    
    return x_edges, y_edges, x_centers, y_centers



def assign_to_grid(x_coords: np.ndarray, y_coords: np.ndarray,
                  x_min: float, y_min: float, grid_size: float,
                  max_x_cells: int) -> np.ndarray:
    """
    Assign points to grid cells.
    
    Parameters:
        x_coords: X coordinates of points
        y_coords: Y coordinates of points
        x_min: Minimum X coordinate of grid
        y_min: Minimum Y coordinate of grid
        grid_size: Size of grid cells
        max_x_cells: Number of cells in X direction (for ID encoding)
        
    Returns:
        Array of grid cell IDs
        
    Example:
        >>> x = np.array([1.5, 6.5, 2.3])
        >>> y = np.array([1.2, 7.8, 3.4])
        >>> grid_ids = assign_to_grid(x, y, 0, 0, 5, 3)
        >>> len(grid_ids) == 3
        True
    """
    # Calculate grid indices
    x_indices = ((x_coords - x_min) / grid_size).astype(int)
    y_indices = ((y_coords - y_min) / grid_size).astype(int)
    
    # Create unique grid cell IDs using original encoding: y_idx * max_x_cells + x_idx
    grid_ids = y_indices * max_x_cells + x_indices
    
    return grid_ids


def create_grid_centroids(x_centers: np.ndarray, y_centers: np.ndarray,
                        z_center: float = 0.0) -> Tuple[Dict[int, Tuple[float, float, float]], Dict[int, int]]:
    """
    Create grid centroids mapping.
    
    Parameters:
        x_centers: X coordinates of grid centers
        y_centers: Y coordinates of grid centers
        z_center: Fixed Z coordinate for all grid centers
        
    Returns:
        Tuple of (grid_centroids dict, grid_id_to_index dict)
        
    Example:
        >>> x_centers = np.array([2.5, 7.5])
        >>> y_centers = np.array([2.5, 7.5])
        >>> centroids, id_map = create_grid_centroids(x_centers, y_centers)
        >>> len(centroids) == 4
        True
    """
    grid_centroids = {}
    grid_id_to_index = {}
    index = 0
    
    # Use original encoding: y_idx * max_x_cells + x_idx
    max_x_cells = len(x_centers)
    
    for y_idx, y_center in enumerate(y_centers):
        for x_idx, x_center in enumerate(x_centers):
            grid_id = y_idx * max_x_cells + x_idx  # Match original encoding
            grid_centroids[grid_id] = (x_center, y_center, z_center)
            grid_id_to_index[grid_id] = index
            index += 1
    
    return grid_centroids, grid_id_to_index


def calculate_cell_centroids(data: pd.DataFrame, 
                           cell_col: str = 'Cell',
                           x_col: str = 'X',
                           y_col: str = 'Y', 
                           z_col: str = 'Z') -> Dict[str, Tuple[float, float, float]]:
    """
    Calculate centroids for each cell.
    
    Parameters:
        data: DataFrame with spatial data
        cell_col: Name of cell column
        x_col: Name of X coordinate column
        y_col: Name of Y coordinate column
        z_col: Name of Z coordinate column
        
    Returns:
        Dictionary mapping cell IDs to (x, y, z) centroids
        
    Example:
        >>> data = pd.DataFrame({
        ...     'Cell': ['C1', 'C1', 'C2'], 
        ...     'X': [1, 2, 5], 'Y': [1, 2, 5], 'Z': [1, 2, 5]
        ... })
        >>> centroids = calculate_cell_centroids(data)
        >>> centroids['C1'][0] == 1.5  # Mean of 1 and 2
        True
    """
    centroids = {}
    
    for cell_id in data[cell_col].unique():
        cell_data = data[data[cell_col] == cell_id]
        centroid_x = cell_data[x_col].mean()
        centroid_y = cell_data[y_col].mean()
        centroid_z = cell_data[z_col].mean()
        centroids[cell_id] = (centroid_x, centroid_y, centroid_z)
    
    return centroids


def get_global_gene_list() -> List[str]:
    """
    Get the standardized list of 101 genes for analysis.
    
    Returns:
        List of gene names
        
    Example:
        >>> genes = get_global_gene_list()
        >>> len(genes) == 101
        True
    """
    return [
        'ACTB', 'APOD', 'Apoe', 'APP', 'AQP4', 'Baiap2l1', 'BSN', 'C1qa', 'C1qc', 'C1ql2',
        'C4a', 'CALB1', 'Casp8', 'Ccrl2', 'Cd68', 'Cd9', 'Clu', 'Creb1', 'Csf1r', 'Cst7',
        'Ctsa', 'Ctsd', 'Ctsh', 'Ctsl', 'Ctss', 'Cx3cr1', 'Cyba', 'DLG4', 'Fcer1g', 'Fcgr3',
        'Foxg1', 'Gfap', 'Gns', 'Gpx4', 'Grin3a', 'Grn', 'Gusb', 'Hexa', 'HIF-1α', 'HSPD1',
        'Id3', 'NOS1', 'ACOD1', 'Itga6', 'Itgam', 'Itgb5', 'Itm2b', 'KCNK2', 'Laptm5', 'Lgals3bp',
        'Ly86', 'MAPT', 'Mbp', 'Mfap4', 'Mpeg1', 'Neurod1', 'Neurod6', 'Nfix', 'NOV', 'Npas3',
        'Npc2', 'Nr3c2', 'Ntf3', 'Olfml3', 'P2ry1', 'Pax6', 'Plek', 'Plp1', 'Prox1', 'PSEN1',
        'PSEN2', 'PTGDS', 'PTN', 'PVALB', 'RAB27A', 'RAB3A', 'Rerg', 'Rest', 'RIMS1', 'RPH3A',
        'S100a6', 'Sdc4', 'Serpina3n', 'Slc1a3', 'SNAP25', 'SNCA', 'Sox2', 'Sox4', 'STX1A', 'STXBP1',
        'STXBP5', 'Syp', 'Tbr1', 'TIM50', 'Tmem119', 'TOMM20', 'TOMM40', 'Trem2', 'Tyrobp', 'UNC13A',
        'Vsir'
    ] 


def create_3d_grid(x_min: float, x_max: float, y_min: float, y_max: float, 
                   z_min: float, z_max: float, grid_size: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Create 3D grid edges and centers for spatial gridding.
    
    Parameters:
        x_min: Minimum X coordinate
        x_max: Maximum X coordinate
        y_min: Minimum Y coordinate
        y_max: Maximum Y coordinate
        z_min: Minimum Z coordinate
        z_max: Maximum Z coordinate
        grid_size: Size of grid cells
        
    Returns:
        Tuple of (x_edges, y_edges, z_edges, x_centers, y_centers, z_centers)
        
    Example:
        >>> x_edges, y_edges, z_edges, x_centers, y_centers, z_centers = create_3d_grid(0, 10, 0, 10, 0, 5, 5)
        >>> len(x_centers) == 2
        True
    """
    # Create 3D grid edges
    x_edges = np.arange(x_min, x_max + grid_size, grid_size)
    y_edges = np.arange(y_min, y_max + grid_size, grid_size)
    z_edges = np.arange(z_min, z_max + grid_size, grid_size)
    
    # Calculate 3D grid centers
    x_centers = x_edges[:-1] + grid_size / 2
    y_centers = y_edges[:-1] + grid_size / 2
    z_centers = z_edges[:-1] + grid_size / 2
    
    return x_edges, y_edges, z_edges, x_centers, y_centers, z_centers


def assign_to_3d_grid(x_coords: np.ndarray, y_coords: np.ndarray, z_coords: np.ndarray,
                      x_min: float, y_min: float, z_min: float, grid_size: float,
                      max_x_cells: int, max_y_cells: int) -> np.ndarray:
    """
    Assign points to 3D grid cells.
    
    Parameters:
        x_coords: X coordinates of points
        y_coords: Y coordinates of points
        z_coords: Z coordinates of points
        x_min: Minimum X coordinate of grid
        y_min: Minimum Y coordinate of grid
        z_min: Minimum Z coordinate of grid
        grid_size: Size of grid cells
        max_x_cells: Number of cells in X direction
        max_y_cells: Number of cells in Y direction
        
    Returns:
        Array of 3D grid cell IDs
        
    Example:
        >>> x = np.array([1.5, 6.5, 2.3])
        >>> y = np.array([1.2, 7.8, 3.4])
        >>> z = np.array([0.5, 2.1, 1.8])
        >>> grid_ids = assign_to_3d_grid(x, y, z, 0, 0, 0, 5, 3, 3)
        >>> len(grid_ids) == 3
        True
    """
    # Calculate 3D grid indices
    x_indices = ((x_coords - x_min) / grid_size).astype(int)
    y_indices = ((y_coords - y_min) / grid_size).astype(int)
    z_indices = ((z_coords - z_min) / grid_size).astype(int)
    
    # Create unique 3D grid cell IDs: z_idx * max_xy_cells + y_idx * max_x_cells + x_idx
    max_xy_cells = max_x_cells * max_y_cells
    grid_ids = z_indices * max_xy_cells + y_indices * max_x_cells + x_indices
    
    return grid_ids


def create_3d_grid_centroids(x_centers: np.ndarray, y_centers: np.ndarray, z_centers: np.ndarray) -> np.ndarray:
    """
    Create centroids for 3D grid cells.
    
    Parameters:
        x_centers: X coordinates of grid centers
        y_centers: Y coordinates of grid centers
        z_centers: Z coordinates of grid centers
        
    Returns:
        Array of shape (n_grid_cells, 3) with (x, y, z) centroids
        
    Example:
        >>> x_centers = np.array([2.5, 7.5])
        >>> y_centers = np.array([2.5, 7.5])
        >>> z_centers = np.array([1.0, 3.0])
        >>> centroids = create_3d_grid_centroids(x_centers, y_centers, z_centers)
        >>> centroids.shape[1] == 3
        True
    """
    centroids = []
    
    for z_center in z_centers:
        for y_center in y_centers:
            for x_center in x_centers:
                centroids.append([x_center, y_center, z_center])
    
    return np.array(centroids) 


def calculate_binary_weights_from_grid_ids(grid_ids: np.ndarray, 
                                         max_x_cells: int, 
                                         max_y_cells: int,
                                         kind: str = "binary_6") -> np.ndarray:
    """
    Build a binary adjacency weights matrix from 3D grid cell IDs.
    
    This function converts grid IDs back to (i,j,k) indices and creates
    a binary adjacency matrix where only immediate neighbors get weight=1.
    
    Parameters
    ----------
    grid_ids : (n,) int ndarray
        Grid cell IDs from assign_to_3d_grid function.
    max_x_cells : int
        Number of cells in X direction.
    max_y_cells : int
        Number of cells in Y direction.
    kind : {"binary_4", "binary_6"}
        - "binary_4": XY plane only (up, down, left, right).
        - "binary_6": full 3D face-adjacency (±x, ±y, ±z).
    
    Returns
    -------
    W : (n,n) ndarray
        Symmetric binary adjacency matrix (0/1), zero diagonal.
        
    Example
    -------
    >>> grid_ids = np.array([0, 1, 5, 6])  # 2x2x1 grid
    >>> W = calculate_binary_weights_from_grid_ids(grid_ids, 2, 2, "binary_4")
    >>> print(W)
    [[0 1 1 0]
     [1 0 0 1] 
     [1 0 0 1]
     [0 1 1 0]]
    """
    n_cells = len(grid_ids)
    max_xy_cells = max_x_cells * max_y_cells
    
    print(f"      🔧 Building binary adjacency matrix for {n_cells} cells...")
    
    # Convert grid IDs back to (i,j,k) indices
    indices_ijk = np.zeros((n_cells, 3), dtype=int)
    
    for idx, grid_id in enumerate(grid_ids):
        # Reverse the encoding: z_idx * max_xy_cells + y_idx * max_x_cells + x_idx
        z_idx = grid_id // max_xy_cells
        remainder = grid_id % max_xy_cells
        y_idx = remainder // max_x_cells
        x_idx = remainder % max_x_cells
        indices_ijk[idx] = [x_idx, y_idx, z_idx]
    
    # Build adjacency matrix using vectorized operations (MUCH faster!)
    W = np.zeros((n_cells, n_cells), dtype=int)
    
    # Define neighbor offsets based on kind
    if kind == "binary_4":
        # Only XY neighbors: ±x, ±y (ignore Z)
        neighbor_offsets = [
            [-1, 0, 0],  # left
            [1, 0, 0],   # right  
            [0, -1, 0],  # down
            [0, 1, 0],   # up
        ]
    elif kind == "binary_6":
        # Full 3D neighbors: ±x, ±y, ±z
        neighbor_offsets = [
            [-1, 0, 0],  # left
            [1, 0, 0],   # right
            [0, -1, 0],  # down
            [0, 1, 0],   # up
            [0, 0, -1],  # back
            [0, 0, 1],   # front
        ]
    else:
        raise ValueError(f"Invalid kind: {kind}. Must be 'binary_4' or 'binary_6'")
    
    # Vectorized approach: create all possible neighbor positions
    print(f"      📍 Checking {len(neighbor_offsets)} neighbor directions...")
    
    # For each cell, check if its neighbors exist in the active grid
    for i in range(n_cells):
        if i % 100 == 0:  # Progress update
            print(f"      📍 Processed {i}/{n_cells} cells...")
            
        # Get current cell position
        current_pos = indices_ijk[i]
        
        # Check each neighbor direction
        for offset in neighbor_offsets:
            neighbor_pos = current_pos + offset
            
            # Find if this neighbor exists in our active grid
            # Use vectorized comparison
            matches = np.all(indices_ijk == neighbor_pos, axis=1)
            neighbor_indices = np.where(matches)[0]
            
            # Set adjacency for all found neighbors
            for neighbor_idx in neighbor_indices:
                if neighbor_idx != i:  # Not self
                    W[i, neighbor_idx] = 1
                    W[neighbor_idx, i] = 1  # Make symmetric
    
    print(f"      ✅ Binary adjacency matrix built!")
    print(f"         Matrix shape: {W.shape}")
    print(f"         Non-zero elements: {np.sum(W > 0)}")
    print(f"         Average neighbors per cell: {np.sum(W > 0) / n_cells:.1f}")
    
    return W

