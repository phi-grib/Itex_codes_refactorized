"""
    This code aims at calculating the gap filling of the predictions of CII for each endpoint.

    Created by: Eric March Vila (eric.march@upf.edu)
    On: 30/03/2023, 15:05 PM
"""

import pandas as pd

from natsort import natsorted
from typing import Tuple

class GapFillingCalculator(object):
    """
        Class that calculates the gap filling after being passed as argument the current substances with annotations in CII.
    """

    def __init__(self, substances_cii: pd.DataFrame, model_predictions: pd.DataFrame):
        """
            Class initializer
        """

        self.subs = self.process_substances_cii(substances_cii)
        self.model_preds = self.process_model_predictions_dataframe(model_predictions)
    
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

    def process_model_predictions_dataframe(self, model_predictions: pd.DataFrame) -> pd.DataFrame:
        """
            Generateds a prepared dataframe to be used for calculations

            :param model_predictions: dataframe containing the predictions for each endpoint and model version.

            :return grouped_ans: dataframe with grouped annotations per model and endpoint
        """

        group = model_predictions.groupby(['endpoint'])['value'].value_counts()
        
        grouped_ans = pd.DataFrame(group.unstack(level=0).to_records())

        grouped_ans.loc[grouped_ans['value'] == 'Uncertain', 'value'] = 'No information'
        grouped_ans.loc[grouped_ans['value'] == 'YES', 'value'] = 'YES + Pending'

        return grouped_ans

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
        
        # Calculate informed space (YES + NO)
        merged_df.loc[3] = merged_df.iloc[[0,2]].sum()
        
        # Clean new value in row 3
        merged_df.iloc[3,0] = 'Informed space'
        
        # Total space is the sum of all annotations from CII.
        # Models can't process all compounds, so if we sum all the annotations obtained from them
        # we will not take into account all the chemical space we have collected
        merged_df.loc[4] = merged_df[cii_endpoint_format].values[0:3].sum()
        merged_df.iloc[4,0] = 'Total space'
        
        # Gap filling calculation in percentage
        merged_df.loc[5] = 'Gap filling (%)'
        merged_df.iloc[5,1:] = (merged_df.iloc[3,1:] / merged_df.iloc[4,1:]) * 100
        
        return merged_df