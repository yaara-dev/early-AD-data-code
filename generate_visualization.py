#!/usr/bin/env python3
"""
Generate interactive spatial visualizations for ExSeq brain AD data.

This script processes spatial transcriptomics data and creates static HTML
files with interactive Plotly visualizations for GitHub Pages deployment.
"""

import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from pathlib import Path


def load_and_process_data(csv_path):
    """Load the CSV data and create downsampled and full datasets."""
    print("Loading data from CSV...")
    df = pd.read_csv(csv_path)
    
    print(f"Loaded {len(df)} rows")
    print(f"Columns: {df.columns.tolist()}")
    
    # Get unique values for filters
    unique_genes = sorted(df['gene'].unique().tolist())
    unique_regions = sorted(df['region'].unique().tolist())
    unique_fovs = sorted(df['fov'].unique().tolist())
    
    # Create downsampled dataset (every 10th point for overview)
    downsample_factor = 10
    df_overview = df.iloc[::downsample_factor].copy()
    
    print(f"Downsampled to {len(df_overview)} points (every {downsample_factor}th)")
    
    return df, df_overview, unique_genes, unique_regions, unique_fovs


def save_overview_data(df_overview, output_dir):
    """Save downsampled overview data as JSON."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Convert DataFrame to JSON
    data_dict = df_overview.to_dict('records')
    
    output_file = output_dir / "overview.json"
    print(f"Saving overview data to {output_file}...")
    
    with open(output_file, 'w') as f:
        json.dump(data_dict, f)
    
    print(f"Saved {len(data_dict)} points to overview.json")
    
    return str(output_file)


def calculate_statistics(df, output_dir):
    """Calculate and save statistics per region."""
    stats = {}
    
    for region in df['region'].unique():
        region_df = df[df['region'] == region]
        stats[region] = {
            'cell_count': region_df['cell_id'].nunique(),
            'unique_genes': region_df['gene'].nunique(),
            'total_points': len(region_df),
            'area': region_df['region_area'].iloc[0],
            'proportion': region_df['region_proportion'].iloc[0]
        }
    
    output_file = Path(output_dir) / "stats.json"
    with open(output_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"Saved statistics for {len(stats)} regions")
    return stats


def create_2d_visualization(df_overview, df_full, genes, regions, fovs):
    """Create a 2D scatter plot visualization with gene and region filters."""
    print("Creating 2D visualization...")
    
    # Sample a subset for initial display (more points for better coverage)
    df_sample = df_overview.sample(min(10000, len(df_overview)))
    
    # Create figure
    fig = go.Figure()
    
    # Create traces for each gene (shown by default)
    gene_colors = px.colors.qualitative.Dark24 + px.colors.qualitative.Light24
    unique_genes_shown = sorted(df_sample['gene'].unique())
    
    for i, gene in enumerate(unique_genes_shown):
        gene_df = df_sample[df_sample['gene'] == gene]
        if len(gene_df) == 0:
            continue
        
        fig.add_trace(go.Scattergl(
            x=gene_df['x_coordinate'],
            y=gene_df['y_coordinate'],
            mode='markers',
            marker=dict(
                size=5,
                color=gene_colors[i % len(gene_colors)],
                opacity=0.6,
                line=dict(width=0)
            ),
            text=[f"Gene: {g}<br>Region: {r}<br>Cell: {c}<br>FOV: {f}<br>Z: {z:.2f}" 
                  for g, r, c, f, z in zip(gene_df['gene'], 
                                            gene_df['region'],
                                            gene_df['cell_id'],
                                            gene_df['fov'],
                                            gene_df['z_coordinate'])],
            hovertemplate='<b>%{text}</b><extra></extra>',
            name=gene,
            visible=True,
            legendrank=i,
            showlegend=True
        ))
    
    # Add dropdown menu for region filtering
    buttons = [dict(
        label="All Regions",
        method="restyle",
        args=[{"x": [gene_df['x_coordinate'].tolist() for gene in unique_genes_shown for gene_df in [df_sample[df_sample['gene'] == gene]] if len(df_sample[df_sample['gene'] == gene]) > 0], 
              "y": [gene_df['y_coordinate'].tolist() for gene in unique_genes_shown for gene_df in [df_sample[df_sample['gene'] == gene]] if len(df_sample[df_sample['gene'] == gene]) > 0]}]
    )]
    
    # Add button for each region
    for region in regions:
        region_df = df_sample[df_sample['region'] == region]
        
        # Get data for each gene in this region
        x_data = []
        y_data = []
        text_data = []
        
        for gene in unique_genes_shown:
            gene_region_df = region_df[region_df['gene'] == gene]
            if len(gene_region_df) > 0:
                x_data.append(gene_region_df['x_coordinate'].tolist())
                y_data.append(gene_region_df['y_coordinate'].tolist())
                text_data.append([f"Gene: {g}<br>Region: {r}<br>Cell: {c}<br>FOV: {f}<br>Z: {z:.2f}" 
                                 for g, r, c, f, z in zip(gene_region_df['gene'], 
                                                          gene_region_df['region'],
                                                          gene_region_df['cell_id'], 
                                                          gene_region_df['fov'],
                                                          gene_region_df['z_coordinate'])])
            else:
                x_data.append([])
                y_data.append([])
                text_data.append([])
        
        buttons.append(dict(
            label=region,
            method="restyle",
            args=[{"x": x_data, "y": y_data, "text": text_data}]
        ))
    
    # Add Select All / Select None buttons for genes
    gene_buttons = [
        dict(
            label="Select All Genes",
            method="update",
            args=[{"visible": [True] * len(fig.data)}]
        ),
        dict(
            label="Select None",
            method="update",
            args=[{"visible": [False] * len(fig.data)}]
        )
    ]
    
    updatemenus = [
        dict(
            buttons=buttons,
            direction="down",
            showactive=True,
            x=0.0,
            xanchor="left",
            y=1.02,
            yanchor="top"
        ),
        dict(
            buttons=gene_buttons,
            direction="down",
            showactive=False,
            x=0.15,
            xanchor="left",
            y=1.02,
            yanchor="top"
        )
    ]
    
    # Update layout
    fig.update_layout(
        title='ExSeq Brain AD - Spatial Transcriptomics (2D View)',
        xaxis_title='X Coordinate',
        yaxis_title='Y Coordinate',
        template='plotly_white',
        height=800,
        hovermode='closest',
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.01,
            bgcolor="rgba(255,255,255,0.8)"
        ),
        updatemenus=updatemenus,
        annotations=[
            dict(
                text='<b>Navigation:</b> <a href="index.html">2D View</a> | <a href="view_3d.html">3D View</a> | <a href="dashboard.html">Dashboard</a>',
                xref="paper", yref="paper",
                x=0.5, y=1.08,
                xanchor='center',
                showarrow=False,
                font=dict(size=12)
            ),
            dict(
                text='<i>Use dropdown menus above to filter by region or genes. Click legend to toggle individual genes.</i>',
                xref="paper", yref="paper",
                x=0.5, y=0.96,
                xanchor='center',
                showarrow=False,
                font=dict(size=11, color="gray")
            )
        ]
    )
    
    return fig


def create_3d_visualization(df_overview, df_full, genes, regions, fovs):
    """Create a 3D scatter plot visualization with region colors."""
    print("Creating 3D visualization...")
    
    # Create figure
    fig = go.Figure()
    
    # Sample subset for initial display
    df_sample = df_overview.sample(min(8000, len(df_overview)))
    
    # Create traces for each region
    region_colors = px.colors.qualitative.Set3
    for i, region in enumerate(regions):
        region_df = df_sample[df_sample['region'] == region]
        if len(region_df) == 0:
            continue
        
        fig.add_trace(go.Scatter3d(
            x=region_df['x_coordinate'],
            y=region_df['y_coordinate'],
            z=region_df['z_coordinate'],
            mode='markers',
            marker=dict(
                size=1.5,
                color=region_colors[i % len(region_colors)],
                opacity=0.6
            ),
            text=[f"Gene: {g}<br>Region: {r}<br>Cell: {c}<br>FOV: {f}" 
                  for g, r, c, f in zip(region_df['gene'], 
                                          region_df['region'],
                                          region_df['cell_id'],
                                          region_df['fov'])],
            hovertemplate='<b>%{text}</b><extra></extra>',
            name=region,
            showlegend=True
        ))
    
    # Update layout
    fig.update_layout(
        title='ExSeq Brain AD - Spatial Transcriptomics (3D View)',
        scene=dict(
            xaxis_title='X Coordinate',
            yaxis_title='Y Coordinate',
            zaxis_title='Z Coordinate',
            camera=dict(
                eye=dict(x=1.5, y=1.5, z=1.2)
            )
        ),
        template='plotly_white',
        height=800,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.01
        ),
        annotations=[
            dict(
                text='<b>Navigation:</b> <a href="index.html">2D View</a> | <a href="view_3d.html">3D View</a> | <a href="dashboard.html">Dashboard</a>',
                xref="paper", yref="paper",
                x=0.5, y=1.02,
                xanchor='center',
                showarrow=False,
                font=dict(size=14)
            )
        ]
    )
    
    return fig


def create_dashboard(df_overview, df_full, stats, genes, regions, fovs):
    """Create a combined dashboard with 2D and 3D views."""
    print("Creating dashboard...")
    
    # Create subplots
    from plotly.subplots import make_subplots
    
    fig = make_subplots(
        rows=2, cols=2,
        specs=[[{"type": "scatter"}, {"type": "table"}],
               [{"type": "scatter3d"}, {"type": "scatter"}]],
        subplot_titles=('2D View', 'Region Statistics', '3D View', 'Cell Counts by Region'),
        vertical_spacing=0.12,
        horizontal_spacing=0.15
    )
    
    # Sample data
    df_sample = df_overview.sample(min(6000, len(df_overview)))
    
    # Color by region for consistent coloring
    region_colors = px.colors.qualitative.Set3
    
    # Add 2D plot - color by region
    for i, region in enumerate(regions):
        region_df = df_sample[df_sample['region'] == region]
        if len(region_df) == 0:
            continue
        
        fig.add_trace(
            go.Scattergl(
                x=region_df['x_coordinate'],
                y=region_df['y_coordinate'],
                mode='markers',
                marker=dict(
                    size=3,
                    color=region_colors[i % len(region_colors)],
                    opacity=0.7,
                    line=dict(width=0)
                ),
                text=[f"Gene: {g}<br>Region: {r}" 
                      for g, r in zip(region_df['gene'], region_df['region'])],
                hovertemplate='<b>%{text}</b><extra></extra>',
                name=region,
                showlegend=True
            ),
            row=1, col=1
        )
    
    # Add statistics table
    table_data = []
    for region in regions:
        if region in stats:
            table_data.append([
                region,
                f"{stats[region]['cell_count']:,}",
                f"{stats[region]['unique_genes']}",
                f"{stats[region]['total_points']:,}"
            ])
    
    fig.add_trace(
        go.Table(
            header=dict(
                values=['Region', 'Cells', 'Genes', 'Points'],
                fill_color='paleturquoise',
                align='left'
            ),
            cells=dict(
                values=list(zip(*table_data)),
                fill_color='lavender',
                align='left'
            )
        ),
        row=1, col=2
    )
    
    # Add 3D plot
    for i, region in enumerate(regions):
        region_df = df_sample[df_sample['region'] == region]
        if len(region_df) == 0:
            continue
        
        fig.add_trace(
            go.Scatter3d(
                x=region_df['x_coordinate'],
                y=region_df['y_coordinate'],
                z=region_df['z_coordinate'],
                mode='markers',
                marker=dict(
                    size=1.5,
                    color=region_colors[i % len(region_colors)],
                    opacity=0.6
                ),
                text=[f"Gene: {g}<br>Region: {r}" 
                      for g, r in zip(region_df['gene'], region_df['region'])],
                hovertemplate='<b>%{text}</b><extra></extra>',
                name=region,
                showlegend=True
            ),
            row=2, col=1
        )
    
    # Create cell counts bar chart instead of gene counts
    cell_counts = []
    for region in regions:
        cell_count = stats[region]['cell_count']
        cell_counts.append(cell_count)
    
    fig.add_trace(
        go.Bar(
            x=regions,
            y=cell_counts,
            marker_color=region_colors[:len(regions)],
            showlegend=False
        ),
        row=2, col=2
    )
    
    # Update axes
    fig.update_xaxes(title_text="X Coordinate", row=1, col=1)
    fig.update_yaxes(title_text="Y Coordinate", row=1, col=1)
    fig.update_xaxes(title_text="Region", row=2, col=2)
    fig.update_yaxes(title_text="Cell Count", row=2, col=2)
    
    fig.update_scenes(
        xaxis_title='X Coordinate',
        yaxis_title='Y Coordinate',
        zaxis_title='Z Coordinate',
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.2)),
        row=2, col=1
    )
    
    # Update layout
    fig.update_layout(
        title='ExSeq Brain AD - Spatial Transcriptomics Dashboard',
        template='plotly_white',
        height=1000,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="right",
            x=1.02,
            bgcolor="rgba(255,255,255,0.8)"
        ),
        annotations=[
            dict(
                text='<b>Navigation:</b> <a href="index.html">2D View</a> | <a href="view_3d.html">3D View</a> | <a href="dashboard.html">Dashboard</a>',
                xref="paper", yref="paper",
                x=0.5, y=1.01,
                xanchor='center',
                showarrow=False,
                font=dict(size=14)
            )
        ]
    )
    
    return fig


def main():
    """Main function to generate all visualizations."""
    # Paths
    csv_path = "fem3_WTE1_B_R_regions_genes.csv"
    output_dir = Path(".")
    data_dir = Path("data")
    data_dir.mkdir(exist_ok=True)
    
    # Load and process data
    df, df_overview, genes, regions, fovs = load_and_process_data(csv_path)
    
    # Save overview data
    save_overview_data(df_overview, data_dir)
    
    # Calculate statistics
    stats = calculate_statistics(df, data_dir)
    
    # Create visualizations
    print("\nGenerating HTML files...")
    
    # 2D View
    fig_2d = create_2d_visualization(df_overview, df, genes, regions, fovs)
    fig_2d.write_html("index.html", include_plotlyjs='cdn')
    print("✓ Created index.html")
    
    # 3D View
    fig_3d = create_3d_visualization(df_overview, df, genes, regions, fovs)
    fig_3d.write_html("view_3d.html", include_plotlyjs='cdn')
    print("✓ Created view_3d.html")
    
    # Dashboard
    fig_dash = create_dashboard(df_overview, df, stats, genes, regions, fovs)
    fig_dash.write_html("dashboard.html", include_plotlyjs='cdn')
    print("✓ Created dashboard.html")
    
    print("\n✓ All visualizations generated successfully!")
    print(f"\nData summary:")
    print(f"  - Total points: {len(df):,}")
    print(f"  - Overview points: {len(df_overview):,}")
    print(f"  - Genes: {len(genes)}")
    print(f"  - Regions: {len(regions)}")
    print(f"  - FOVs: {len(fovs)}")
    print(f"\nRegions: {', '.join(regions)}")


if __name__ == "__main__":
    main()

