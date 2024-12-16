"""
    This code aims at calculating the gap filling of the predictions of CII for each endpoint.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 30/03/2023, 15:05 PM
"""

import numpy as np
import pandas as pd

from natsort import natsorted
from typing import Tuple, Union

class GapFillingCalculator(object):
    """
        Class that calculates the gap filling after being passed as argument the current substances with annotations in CII.
    """

    def __init__(self, substances_cii: pd.DataFrame, model_predictions: pd.DataFrame):
        """
            Class initializer
        """

        self.subs = self.process_substances_cii(substances_cii)
        self.model_preds = self.process_model_predictions_dataframe(substances_cii, model_predictions)

    def process_substances_cii(self, substances_cii: pd.DataFrame) -> pd.DataFrame:
        """
            Generates a prepared dataframe to be used for calculations.

            :param substances_cii: cii substances dataframe

            :return eps_clean: cii substances dataframe cleaned
        """

        eps = substances_cii[['cmr','pbt','vpvb','endocrine_disruptor','c','m','r','p','b','t','vp','vb',
                                'androgen_rc','estrogen_rc','glucocorticoid_rc','thyroid_rc']]
        
        eps_clean = eps.apply(lambda x: x.value_counts()).reset_index().rename(columns={'index':'value'})

        eps_clean.fillna(0,inplace=True)

        eps_clean.loc[3] += eps_clean.loc[2]

        eps_clean.loc[eps_clean['value'] == 'YESPending', 'value'] = 'YES + Pending'

        eps_clean.drop(2, inplace=True)

        return eps_clean

    def process_model_predictions_dataframe(self, substances_cii: pd.DataFrame, model_predictions: pd.DataFrame) -> pd.DataFrame:
        """
            Generateds a prepared dataframe to be used for calculations

            :param model_predictions: dataframe containing the predictions for each endpoint and model version.

            :return grouped_ans: dataframe with grouped annotations per model and endpoint
        """

        new_predicted_compounds = self.get_new_predicted_compounds(substances_cii, model_predictions)
        group = new_predicted_compounds.groupby(['endpoint'])['value'].value_counts()
        
        grouped_ans = pd.DataFrame(group.unstack(level=0).to_records())

        grouped_ans.loc[grouped_ans['value'] == 'Uncertain', 'value'] = 'No information'
        grouped_ans.loc[grouped_ans['value'] == 'YES', 'value'] = 'YES + Pending'

        return grouped_ans

    def get_new_predicted_compounds(self, substances_cii: pd.DataFrame, model_predictions: pd.DataFrame) -> pd.DataFrame:
        """
            Gets the compounds that were non informed and now have a predicted annotations.
            Removes compounds previously informed in CII.

            :param model_precitions: dataframe containing the predictions for each endpoint and model version

            :return new_predicted_compounds: dataframe containing only newly annotated compounds
        """

        cmr_non_informed_cas = substances_cii.loc[substances_cii['cmr'] == 'No information', 'CAS'].values
        pbt_non_informed_cas = substances_cii.loc[substances_cii['pbt'] == 'No information', 'CAS'].values
        vpvb_non_informed_cas = substances_cii.loc[substances_cii['vpvb'] == 'No information', 'CAS'].values
        ed_non_informed_cas = substances_cii.loc[substances_cii['endocrine_disruptor'] == 'No information', 'CAS'].values

        new_predicted_compounds = model_predictions.loc[((model_predictions['endpoint'].str.contains('CMR')) &
                                                        (model_predictions['name'].isin(cmr_non_informed_cas))) |
                                                        ((model_predictions['endpoint'].str.contains('PBT')) &
                                                        (model_predictions['name'].isin(pbt_non_informed_cas))) |
                                                        ((model_predictions['endpoint'].str.contains('vPvB')) &
                                                        (model_predictions['name'].isin(vpvb_non_informed_cas))) |
                                                        ((model_predictions['endpoint'].str.contains('ED')) &
                                                        (model_predictions['name'].isin(ed_non_informed_cas)))]

        return new_predicted_compounds

    def get_cols_endpoint(self, df: pd.DataFrame, endpoint: str) -> list:
        """
            Returns a list of alphabetically ordered columns.

            :param df: substances dataframe
            :param endpoint: endpoint to take into account

            :return sorted_cols: list of sorted columns
        """

        cols = []
        for col in df.columns:
            upper_col = col.upper()
            if endpoint in upper_col:
                cols.append(col)
        
        sorted_cols = natsorted(cols, key=lambda x: x.split('_')[-1])
        sorted_cols.insert(0,'value')
        
        return sorted_cols

    def endpoint_check(self, endpoint: str) -> Tuple[str]:
        """
            Checks the endpoint name. Specially for endocrine_disruptor, which is stored using different names (ed, endocrine_disruptor)

            :param endpoint: endpoint to check

            :return endpoint_cii, endpoint_preds: endpoints to refer to cii dataframe and prediction dataframe
        """

        if endpoint == 'ed':
            endpoint_cii = 'endocrine_disruptor'
            endpoint_preds = 'ed'
        else:
            endpoint_cii = endpoint_preds = endpoint
            
        return endpoint_cii, endpoint_preds
    
    def process_column_no_info(self, col: np.ndarray, cii_no_info: Union[int, float]) -> np.ndarray:
        """
            Process a single column of the DataFrame.
            
            Args:
            col (np.ndarray): The column to process.
            cii_no_info (Union[int, float]): The 'No information' value from the cii column.
            
            Returns:
            np.ndarray: The processed column.
            
            This function does the following:
            1. Sums the 'NO' and 'YES + Pending' values (index 0 and 2).
            2. Subtracts this sum from the cii 'No information' value.
            3. Adds the result to the column's 'No information' value (index 1).
        """

        no_yes_sum = col[0] + col[2]  # Sum of NO and YES + Pending
        col[1] = cii_no_info - no_yes_sum # Adjust No information

        return col

    def add_cii_values(self, col: np.ndarray, cii_no: Union[int, float], cii_yes_pending: Union[int, float]) -> np.ndarray:
        """
            Add cii NO and YES + Pending values to the corresponding cells in the column.
            
            Args:
            col (np.ndarray): The column to process.
            cii_no (Union[int, float]): The 'NO' value from the cii column.
            cii_yes_pending (Union[int, float]): The 'YES + Pending' value from the cii column.
            
            Returns:
            np.ndarray: The processed column.
        """

        col[0] += cii_no  # Add cmr_cii NO to NO
        col[2] += cii_yes_pending  # Add cmr_cii YES + Pending to YES + Pending
        
        return col

    def add_total_space (self, col: np.ndarray) -> np.ndarray:
        """
            Calculate the total space for a given column.

            This function sums up the values in the first three elements of the input array,
            which typically represent 'NO', 'No information', and 'YES + Pending' categories.

            Args:
            self: The instance of the class this method belongs to.
            col (np.ndarray): The column to process, expected to be a NumPy array.

            Returns:
            np.ndarray: A single-element array containing the total space (sum of the first three elements).

            Note:
            - This function assumes that the input array has at least three elements.
            - The returned value represents the total chemical space covered by all categories.
        """

        total_space = col[0] + col[1] + col[2]

        return total_space

    def merge_source_ans_and_preds(self, endpoint: str) -> pd.DataFrame:
        """
            Merges the source annotations with the predictions in one dataframe to ultimately calculate the gap filling
            of the predictions when compared with the original annotations.

            :param endpoint: endpoint to select.

            :return merged_df: dataframe with the gap filling calculated.
        """

        cii_endpoint, pred_endpoint = self.endpoint_check(endpoint)
        merged_df = self.model_preds[self.get_cols_endpoint(self.model_preds,pred_endpoint.upper())]
        cii_endpoint_format = '_'.join([cii_endpoint.lower(),'cii'])
        
        merged_df.insert(1, cii_endpoint_format, self.subs[cii_endpoint.lower()].values)
        
        # Remove NaNs for 0 to avoid sum problems
        merged_df = merged_df.fillna(0)

        # Calculate the reduction of the No information in cii for each column
        # It uses the current NO and YES + Pending to subtract to No information
        
        # Get the necessary values from cii
        cii_no_info: float = merged_df.loc[merged_df['value'] == 'No information', [col for col in merged_df.columns if 'cii' in col][0]].values[0]
        cii_no: float = merged_df.loc[merged_df['value'] == 'NO', [col for col in merged_df.columns if 'cii' in col][0]].values[0]
        cii_yes_pending: float = merged_df.loc[merged_df['value'] == 'YES + Pending', [col for col in merged_df.columns if 'cii' in col][0]].values[0]

        # Apply the functions to each column except 'value' and 'cmr_cii'
        for column in merged_df.columns[2:]:  # Start from the third column
            merged_df[column] = self.process_column_no_info(merged_df[column].values, cii_no_info)
            merged_df[column] = self.add_cii_values(merged_df[column].values, cii_no, cii_yes_pending)

        # Calculate informed space (YES + NO)
        merged_df.loc[3] = merged_df.iloc[[0,2]].sum()
        
        # Clean new value in row 3
        merged_df.iloc[3,0] = 'Informed space'

        # Total space is the sum of all annotations from CII.
        # Models can't process all compounds, so if we sum all the annotations obtained from them
        # we will not take into account all the chemical space we have collected
        merged_df.loc[4] = [self.add_total_space(merged_df[column].values) for column in merged_df.columns]
        merged_df.iloc[4,0] = 'Total space'

        # Gap filling calculation in percentage
        merged_df.loc[5] = 'Gap filling (%)'
        merged_df.iloc[5,1:] = (merged_df.iloc[3,1:] / merged_df.iloc[4,1:]) * 100
        
        return merged_df