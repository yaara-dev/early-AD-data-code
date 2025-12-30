# Data Structure Documentation

This document describes the structure, columns, and relationships of the spatial transcriptomics data files used in the ExSeq Brain AD study.

## Overview

The data is stored in a **normalized CSV format** with **5 files per sample**. This normalized structure separates different types of information into distinct files. Each sample (e.g., `WT_1`, `5xFAD_1.1`) has its own set of 5 files located in `data/all_samples/`.

### Coordinate System

**Important**: All spatial coordinates (X, Y, Z) are in **pre-expansion units** (micrometers). The coordinates have been scaled by dividing post-expansion measurements by 3.5 to convert them back to pre-expansion physical dimensions.

---

## File Structure

Each sample consists of 5 CSV files:

1. **`*_points.csv`** - Main point data with spatial coordinates and gene information
2. **`*_regions.csv`** - Brain region metadata
3. **`*_points_regions.csv`** - Mapping between points and regions
4. **`*_cells.csv`** - Cell identifier information
5. **`*_points_cells.csv`** - Mapping between points and cells, including cell type annotations

---

## File 1: `*_points.csv`

**Purpose**: Contains the core spatial transcriptomics data - each row represents a single RNA molecule detection with its spatial location and gene identity.

**Columns**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `point_id` | Integer | Unique identifier for each point (primary key) | `1`, `2`, `3` |
| `gene` | String | Gene name/symbol detected at this point | `ACTB`, `SNAP25`, `APP` |
| `Z` | Float | Z-coordinate in micrometers (pre-expansion units, depth) | `1.02`, `3.56`, `7.21` |
| `X` | Float | X-coordinate in micrometers (pre-expansion units, horizontal) | `473.63`, `429.00`, `403.09` |
| `Y` | Float | Y-coordinate in micrometers (pre-expansion units, vertical) | `536.83`, `515.51`, `511.77` |
| `fov` | String | Field-of-View identifier where this point was detected | `4w_fem2_WT_F3_B_left_FOV_23F04` |

**Notes**:
- Each row represents one RNA molecule detection
- Multiple points can have the same gene (multiple detections of the same gene)
- Multiple points can share the same spatial coordinates (different genes at same location)
- `point_id` is unique within each sample
- Coordinates are in pre-expansion micrometers (scaled by 3.5 from post-expansion measurements)

**Example rows**:
```csv
point_id,gene,Z,X,Y,fov
1,PTN,1.0182857142857142,473.6346571428571,536.8265714285715,4w_fem2_WT_F3_B_left_FOV_23F04
2,ACTB,0.9283714285714285,429.0043428571429,515.5125714285715,4w_fem2_WT_F3_B_left_FOV_23F04
3,PSEN2,0.9211142857142857,403.09319999999997,511.77348571428575,4w_fem2_WT_F3_B_left_FOV_23F04
```

---

## File 2: `*_regions.csv`

**Purpose**: Contains metadata about hippocampal regions/anatomical areas in the sample.

**Columns**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `region_id` | Integer | Unique identifier for each region (primary key) | `1`, `2`, `3` |
| `region_name` | String | Name of the hippocampal region | `DG-CA1`, `CA1`, `SLM`, `unassigned` |
| `region_area` | Float | Area of the region in square micrometers | `2336905.43`, `1572711.11` |
| `region_proportion` | Float | Proportion of total sample area covered by this region (0.0-1.0) | `0.3015`, `0.2029` |

**Notes**:
- Each sample typically has 8-10 regions
- Some regions may have empty `region_area` and `region_proportion` (e.g., "unassigned" region)
- Common region names include: `DG-CA1`, `CA1`, `CA3`, `SLM`, `SO`, `DG`, `Hilus`, `unassigned`

**Example rows**:
```csv
region_id,region_name,region_area,region_proportion
1,DG-CA1,2336905.4302,0.3015
2,unassigned,,
3,SLM,1572711.1053,0.2029
4,CA3,682298.5665,0.0880
```

---

## File 3: `*_points_regions.csv`

**Purpose**: Maps each point to hippocampal regions.

**Columns**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `point_id` | Integer | Foreign key to `*_points.csv` | `1`, `2`, `3` |
| `region_id` | Integer | Foreign key to `*_regions.csv` | `1`, `3`, `6` |


**Example rows**:
```csv
point_id,region_id
1,1
2,1
3,1
4,1
```

**To join with points and regions**:
```sql
SELECT p.*, r.region_name 
FROM points p
JOIN points_regions pr ON p.point_id = pr.point_id
JOIN regions r ON pr.region_id = r.region_id
```

---

## File 4: `*_cells.csv`

**Purpose**: Contains cell identifier information. Each row represents a unique cell.

**Columns**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `cell_id` | Integer | Unique identifier for each cell (primary key) | `1`, `2`, `3` |
| `cell` | String | Original cell identifier/name from the analysis | `F04_2`, `F04_0`, `H05_1` |

**Notes**:
- Each row represents one unique cell
- `cell_id` is unique within each sample
- The `cell` column contains the original cell identifier (includes FOV prefix like `F04_`, `H05_`, etc.)

**Example rows**:
```csv
cell_id,cell
1,F04_2
2,F04_0
3,F04_4
4,F04_1
5,F04_7
```

---

## File 5: `*_points_cells.csv`

**Purpose**: Maps each point to a cell and includes cell type annotation.

**Columns**:

| Column | Type | Description | Example |
|--------|------|-------------|---------|
| `point_id` | Integer | Foreign key to `*_points.csv` | `1`, `2`, `3` |
| `cell_id` | Integer | Foreign key to `*_cells.csv` | `1`, `2`, `3` |
| `cell_type` | String | Annotated cell type for this point-cell assignment | `endothelial`, `activated microglia`, `astrocytes`, `unassigned` |

**Notes**:
- `cell_type` provides the cell type annotation (e.g., `endothelial`, `activated microglia`, `astrocytes`, `excitatory neurons`, `inhibitory neurons`, `microglia`, `oligodendrocytes`, `unassigned`)
- Every point should have exactly one cell assignment

**Example rows**:
```csv
point_id,cell_id,cell_type
1,1,endothelial
2,2,unassigned
3,2,unassigned
4,2,unassigned
```

**To join with points and cells**:
```sql
SELECT p.*, c.cell, pc.cell_type
FROM points p
JOIN points_cells pc ON p.point_id = pc.point_id
JOIN cells c ON pc.cell_id = c.cell_id
```

---

## File Relationships

The files are related through the following structure:

```
points (point_id) 
    ├── points_regions (point_id → region_id) → regions (region_id)
    └── points_cells (point_id → cell_id) → cells (cell_id)
```

### Relationship Diagram

```
┌─────────────┐
│   points    │
│ (point_id)  │
└──────┬──────┘
       │
       ├─────────────────┐
       │                 │
       ▼                 ▼
┌──────────────────┐  ┌──────────────────┐
│ points_regions   │  │  points_cells    │
│ point_id         │  │  point_id        │
│ region_id        │  │  cell_id         │
└────────┬─────────┘  │  cell_type       │
         │            └────────┬─────────┘
         │                     │
         ▼                     ▼
┌─────────────┐      ┌─────────────┐
│   regions   │      │    cells    │
│ (region_id) │      │  (cell_id)  │
└─────────────┘      └─────────────┘
```


## Complete Data Join Example

To get all information for a point (coordinates, gene, region, cell, cell type):

```sql
SELECT 
    p.point_id,
    p.gene,
    p.X, p.Y, p.Z,
    p.fov,
    r.region_name,
    c.cell,
    pc.cell_type
FROM points p
LEFT JOIN points_regions pr ON p.point_id = pr.point_id
LEFT JOIN regions r ON pr.region_id = r.region_id
LEFT JOIN points_cells pc ON p.point_id = pc.point_id
LEFT JOIN cells c ON pc.cell_id = c.cell_id
WHERE p.point_id = 1;
```

Or in Python/pandas:
```python
import pandas as pd

# Load all files
points = pd.read_csv('WT_1_points.csv')
regions = pd.read_csv('WT_1_regions.csv')
points_regions = pd.read_csv('WT_1_points_regions.csv')
cells = pd.read_csv('WT_1_cells.csv')
points_cells = pd.read_csv('WT_1_points_cells.csv')

# Join to get complete information
data = points.merge(
    points_regions, on='point_id', how='left'
).merge(
    regions, on='region_id', how='left'
).merge(
    points_cells, on='point_id', how='left'
).merge(
    cells, on='cell_id', how='left'
)
```

---

## Data Statistics (Example: WT_1)

- **Total points**: ~560,000 per sample
- **Unique genes**: ~100 genes per sample
- **Regions**: 8 regions per sample
- **Cells**: ~6,000-7,000 cells per sample
- **Coordinate ranges**:
  - X: ~100-1200 μm
  - Y: ~0-1000 μm  
  - Z: ~0.1-23 μm

---

## Common Cell Types

The `cell_type` column in `*_points_cells.csv` typically includes:
- `endothelial`
- `activated microglia`
- `astrocytes`
- `excitatory neurons`
- `inhibitory neurons`
- `microglia`
- `oligodendrocytes`
- `unassigned`

---

## Sample Files

All sample data files are located in `data/all_samples/`. Each sample has 5 files:
- `{SAMPLE_NAME}_points.csv`
- `{SAMPLE_NAME}_regions.csv`
- `{SAMPLE_NAME}_points_regions.csv`
- `{SAMPLE_NAME}_cells.csv`
- `{SAMPLE_NAME}_points_cells.csv`

Available samples:
- `WT_1`, `WT_2.1`, `WT_2.2`, `WT_3`
- `5xFAD_1.1`, `5xFAD_1.2`, `5xFAD_2`, `5xFAD_3`

