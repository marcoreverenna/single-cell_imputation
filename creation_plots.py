#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module containing the functions for testing Beagle imputation algorithm on single cell SNP array data

__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""

# load the packages
import os
import logging
import pandas as pd
import plotly.subplots as sp
import plotly.graph_objs as go
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def violin_plot_chr(file_path, variable, variable_title, y_max=1, y_min=0, colors=('salmon', 'mediumseagreen')):
    """
    Create violin plots for Jaccard and Recall scores.

    Args:
        file_path (str): Path to the Excel file containing the data.
        variable (str): Column name to be used for the x-axis.
        variable_title (str): Title to be used in the plot.
        output_dir (str): Directory to save the plot.
        y_max (float, optional): Maximum value for y-axis. Defaults to 1.
        y_min (float, optional): Minimum value for y-axis. Defaults to 0.
        colors (tuple, optional): Colors for the plots. Defaults to ('salmon', 'mediumseagreen').
    """

    table_conc = pd.read_excel(file_path)

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True,
                        subplot_titles=['Jaccard Scores', 'Recall Scores'],
                        vertical_spacing=0.07)

    x_labels = [str(i) for i in range(1, 23)]  # Define chromosome numbers as x-axis labels

    # Jaccard Scores Plot
    for score_type, color in zip(['j_score_pre', 'j_score_post'], colors):
        fig.add_trace(go.Violin(x=table_conc[variable], y=table_conc[score_type],
                                box_visible=False, line_color=color,
                                meanline_visible=True, opacity=1,
                                name=score_type.upper().replace('_', ' ')),
                      row=1, col=1)

    # Recall Scores Plot
    for score_type, color in zip(['recall_pre', 'recall_post'], colors):
        fig.add_trace(go.Violin(x=table_conc[variable], y=table_conc[score_type],
                                box_visible=False, line_color=color,
                                meanline_visible=True, opacity=1,
                                name=score_type.upper().replace('_', ' '),
                                showlegend=False),  # Hide legend for these traces
                      row=2, col=1)

    # Update layout
    fig.update_layout(height=600, width=900,
                      title=f'Comparison of Jaccard scores and Recall scores pre and post imputation in {variable_title}',
                      legend=dict(orientation="h", yanchor="bottom", y=1.04, xanchor="right", x=1),
                      xaxis=dict(tickmode='array', tickvals=list(range(1, 23)), ticktext=x_labels),
                      xaxis2=dict(tickmode='array', tickvals=list(range(1, 23)), ticktext=x_labels))

    fig.update_yaxes(range=[y_min, y_max], row=1, col=1)
    fig.update_yaxes(range=[y_min, y_max], row=2, col=1)

    # Save the plot
    output_dir = f'plots/'

    output_file = os.path.join(output_dir, f'violinplot_{variable_title}.pdf')
    fig.write_image(output_file)
    print(f"Plot saved to {output_file}")



if __name__ == "__main__":

    try:
        logging.info("Running module: violin_plot_chr")
        violin_plot_chr('results/tables/similarity_recall_concatenated.xlsx', 'chromosome', 'Chromosome', y_max=1, y_min=0)
        
    except Exception as e:
        logging.error(f"Error in violin_plot_chr: {e}")
