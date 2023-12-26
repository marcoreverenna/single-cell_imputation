#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Functions to run the data analysis which includes descriptive analysis, statistical test and correlation analysis.
The descriptive analysis calculates the statistics for the metrics pre and post imputation.
The statistical test (t-test paired) to determine if the differences between pre and post imputation values are statistically significant.
A correlation analysis to check if there is a correlation between the improvement of the metrics and other variables (chromosomes)


__author__ = Marco Reverenna
__copyright__ = Copyright 2022-2023
__version__ = 1.0
__maintainer__ = Marco Reverenna
__email__ = marcoreverenna@gmail.com
__status__ = Dev
"""


import logging
import pandas as pd
import os
import numpy as np
from scipy import stats
from scipy.stats import pearsonr

def calculate_statistics(file_path, group, stat):
    """ Calculate basic statistics such as mean, median, max, min, and standard deviation for specified metrics.

    Args:
        file_path (str): The path to the Excel file containing the data.
        group (str): The column to group by, either 'sample_name' or 'chromosome'.
        stat (str): The statistic to calculate, either 'j_score' or 'recall'.
    """
    # Load the table from the specified Excel file
    table = pd.read_excel(file_path)

    # Group the data and calculate all the statistics at once
    grouped_stats = table.groupby(group)[[f'{stat}_pre', f'{stat}_post']].agg(['mean', 'median', 'max', 'min', lambda x: round(np.std(x, ddof=0), 3)])
    
    # Flatten the MultiIndex and rename columns to make them more meaningful
    grouped_stats.columns = ['_'.join(col).strip() for col in grouped_stats.columns.values]
    rename_dict = {
        f'{stat}_pre_mean': f'mean_{stat}_pre',
        f'{stat}_post_mean': f'mean_{stat}_post',
        f'{stat}_pre_median': f'median_{stat}_pre',
        f'{stat}_post_median': f'median_{stat}_post',
        f'{stat}_pre_max': f'max_{stat}_pre',
        f'{stat}_post_max': f'max_{stat}_post',
        f'{stat}_pre_min': f'min_{stat}_pre',
        f'{stat}_post_min': f'min_{stat}_post',
        f'{stat}_pre_<lambda_0>': f'std_{stat}_pre',
        f'{stat}_post_<lambda_0>': f'std_{stat}_post',
    }
    grouped_stats.rename(columns=rename_dict, inplace=True)

    # Define the directory to save the results
    results_dir = os.path.join(os.path.dirname(file_path), 'statistics')
    os.makedirs(results_dir, exist_ok=True)
    output_file = f'results/statistics/{group}_{stat}_statistics.xlsx'
    grouped_stats.to_excel(output_file, index=True)

    print(f"Statistics saved to {output_file}")
    return grouped_stats




def paired_ttest(file_path, pre_column, post_column, stat):
    """
    Perform a paired t-test on two related samples of scores, pre and post.

    Parameters:
    file_path (str): The path to the Excel file containing the data.
    pre_column (str): The name of the column containing the pre-imputation scores.
    post_column (str): The name of the column containing the post-imputation scores.
    stat (str): The name of the statistic, used for naming the output file.

    Returns:
    DataFrame: A DataFrame containing the t-statistic and p-value.
    """

    table = pd.read_excel(file_path)
    t_statistic, p_value = stats.ttest_rel(table[pre_column], table[post_column])
    results = pd.DataFrame({
        't_statistic': [t_statistic],
        'p_value': [p_value]
    })

    results_dir = os.path.join('results', 'statistics')
    os.makedirs(results_dir, exist_ok=True)

    # Define the output file path
    output_file_path = os.path.join(results_dir, f'paired_test_{stat}.xlsx')

    # Save the results to an Excel file
    results.to_excel(output_file_path, index=False)

    print(f"Paired t-test results saved to {output_file_path}")

    return results


def calculate_correlation(file_path, score_pre, score_post, other_variable):
    """
    Calculate the correlation between the improvement of a metric and another variable.

    Args:
        file_path (str): The path to the Excel file containing the data.
        score_pre (str): The column name of the pre-imputation scores.
        score_post (str): The column name of the post-imputation scores.
        other_variable (str): The column name of the other variable to check for correlation.
    
    Returns:
        float: The Pearson correlation coefficient.
    """
    # Load the data
    df = pd.read_excel(file_path)

    # Calculate the improvement
    df['improvement'] = df[score_post] - df[score_pre]

    # Calculate the Pearson correlation coefficient
    correlation, _ = pearsonr(df['improvement'], df[other_variable])
    
    output_dir = f'results/statistics/correlation_{other_variable}.xlsx'

    os.makedirs(output_dir, exist_ok=True)

    # Define the output file path
    output_file_path = os.path.join(output_dir, f'correlation_{other_variable}.xlsx')
    pd.DataFrame({'Correlation': [correlation]}).to_excel(output_file_path, index=False)

    print(f"Correlation result saved to {output_file_path}")


    return correlation




if __name__ == "__main__":
    try:
        logging.info("Run function calculate_statistics")
        calculate_statistics('results/tables/similarity_recall_concatenated.xlsx', 'sample_name', 'j_score')
        calculate_statistics('results/tables/similarity_recall_concatenated.xlsx', 'sample_name', 'recall')
        calculate_statistics('results/tables/similarity_recall_concatenated.xlsx', 'chromosome', 'j_score')
        calculate_statistics('results/tables/similarity_recall_concatenated.xlsx', 'chromosome', 'recall')
    except Exception as e:
        logging.error(f"Error in calculate_statistics:{e}")
    
    #try:
    #    logging.info("Run function paired_ttest")
    #    paired_ttest('results/tables/similarity_recall_concatenated.xlsx', 'j_score_pre', 'j_score_post', 'jscore')
    #    paired_ttest('results/tables/similarity_recall_concatenated.xlsx', 'recall_pre', 'recall_post', 'recall')

    #except Exception as e:
    #    logging.error(f"Error in paired_ttest:{e}")

    #try:
    #    logging.info("Run function calculate_correlation")
    #    calculate_correlation(file_path='results/tables/similarity_recall_concatenated.xlsx',
    #                         score_pre = 'j_score_pre',
    #                         score_post = 'j_score_post',
    #                         other_variable = 'chromosome')
    #except Exception as e:
    #    logging.error(f"Error in calculate_correlation:{e}")
