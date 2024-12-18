"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 08/04/2020, 12:04 AM
"""

import numpy as np
import pandas as pd
import sys

from typing import Optional

from UpdateDB_ref.Update_CII import UpdateDB

class Endpoint(UpdateDB):
    """
        Child class of UpdateDB. Uses Conenctor functions to interact with CII and UpdateDB functions to
        insert into CII new annotations if needed.
        This class aims to check in CII the presence of hazard annotations for each substance and generate 
        a new annotation for a given endpoint related to a certain hazard.
        The theoretical background comes from the USC Workflow to generate the endpoint annotations.

        Example:
        Formaldehyde has Carc. 2 hazard annotation. This means is positive for CMR, then we will have a YES annotation
        for formaldehyde in CMR.
        If no annotation is found, then No information is returned.
    """

    def __init__(self, host: str = None, dbname: str = None, user: str = None, password: str = None):
        """
            Initializes class with main arguments for psycopg2 connection to the database

            :param host:
            :param dbname:
            :param user:
            :param password:
        """

        if 'cii' in dbname.lower() or 'inventory' in dbname.lower():
            super().__init__(host, dbname, user, password)
            super().open_connection()
            self.db_tag = 'cii'
        elif dbname.lower() in ['cr']:
            super().compoundDB_connection()
            self.db_tag = 'cr'
    
    def get_annotations_for_chem_id(self, chemical_ids: np.ndarray, endpoint_annotations: dict) -> pd.DataFrame:
        """
            For each chemical id in the list, checks regulations table if there are certain annotations,
            which are the values of the dict endpoint_annotations. 

            :param chemical_ids:
            :param endpoint_annotations: dictionary which keys are endpoints (CMR, PBT...) and values are
                                        the hazard annotations associated to that endpoints

            :return chemid_endpoint_annotations:
        """

        chemid_endpoint_annotations = pd.DataFrame(index=range(len(chemical_ids)))
       
        for i, id_ in enumerate(chemical_ids):
            chemid_endpoint_annotations.loc[chemid_endpoint_annotations.index == i, 'chem_id'] = id_
            for endpoint in endpoint_annotations.keys():
                annotations = endpoint_annotations[endpoint]
                final_annotation, source, hazard = self.get_annotation_per_endpoint(id_, endpoint, annotations)
                chemid_endpoint_annotations.loc[chemid_endpoint_annotations.index == i, '{}_ECHA'.format(endpoint)] = final_annotation
                chemid_endpoint_annotations.loc[chemid_endpoint_annotations.index == i, '{}_source_ECHA'.format(endpoint)] = str(source)
                chemid_endpoint_annotations.loc[chemid_endpoint_annotations.index == i, '{}_hazard_ECHA'.format(endpoint)] = str(hazard)
        return chemid_endpoint_annotations

    def get_annotation_per_endpoint(self, chem_id: int, endpoint: str, annotations: list) -> str:
        """
            Checks if there are annotations for that substance id

            :param chem_id:
            :param annotations:

            :return final_annotation:
        """

        chemid_annotations = self.check_presence_in_table(chem_id, annotations)
        
        if chemid_annotations.empty:
            final_annotation = source = hazard = 'No information'
        else:
            if self.db_tag == 'cii':
                final_annotation, source, hazard = self.check_source_of_annotation(endpoint, chemid_annotations)
            elif self.db_tag == 'cr':
                final_annotation = self.check_cr_source(endpoint, chemid_annotations)

        return final_annotation, source, hazard

    def get_total_annotations_per_endpoint(self, chemid_endpoint_annotations: pd.DataFrame) -> pd.DataFrame:
        """
            Calculates the total number of annotations per endpoint in the input dataframe

            :param chemid_endpoint_annotations:

            :return total_annotations_endpoint:
        """

        endpoints = chemid_endpoint_annotations.columns[1:]
        total_annotations_endpoint = pd.DataFrame(index=range(len(endpoints)))

        for i, endpoint in enumerate(endpoints):
            yes_count = len(chemid_endpoint_annotations.loc[chemid_endpoint_annotations[endpoint] == 'YES',endpoint])
            pen_count = len(chemid_endpoint_annotations.loc[chemid_endpoint_annotations[endpoint] == 'Pending',endpoint])
            no_count = len(chemid_endpoint_annotations.loc[chemid_endpoint_annotations[endpoint] == 'NO',endpoint])
            no_info_count = len(chemid_endpoint_annotations.loc[chemid_endpoint_annotations[endpoint] == 'No information',endpoint])

            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'Endpoints'] = endpoint
            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'YES'] = yes_count
            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'Pending'] = pen_count
            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'NO'] = no_count
            total_annotations_endpoint.loc[total_annotations_endpoint.index == i, 'No information'] = no_info_count

        total_annotations_endpoint.loc[:, 'Total annotations'] = total_annotations_endpoint.sum(axis=1)
        
        return total_annotations_endpoint

    def check_presence_in_table(self, chem_id: int, annotations: str) -> pd.DataFrame:
        """
            Ask CII/CR if there are annotations for the input substance

            :param chem_id:
            :param annotations:

            :return substance_annotations:
        """

        if self.db_tag == 'cii':
            query_ = """SELECT reg.id, ci.id as chem_id, rco.country, rt."type", rg.general_regulation_name, 
                        rspec.specific_regulation_name, rsub.subspecific_regulation_name, 
                        rsc.special_cases_name, addr.additional_information_name, regn.names
                        FROM regulations reg
                        left join chem_id ci on ci.id = reg.chem_id 
                        LEFT JOIN substance sub ON sub.chem_id = ci.id
                        left join regulation_country rco on rco.id = reg.reg_country_id
                        left join regulation_type rt on rt.id = reg.reg_type_id
                        left join general_regulation rg on rg.id = reg.gen_reg_id
                        left join specific_regulation rspec on rspec.id = reg.spec_reg_id
                        LEFT JOIN subspecific_regulation rsub ON rsub.id = reg.subspec_reg_id
                        left join special_cases_regulation rsc on rsc.id = reg.special_cases_id
                        left join additional_information_regulation addr on addr.id = reg.additional_information_id
                        LEFT JOIN chem_type ct ON ct.id = reg.chem_type_id
                        LEFT JOIN regulation_names regn ON regn.id = reg.regulation_id 
                        WHERE reg.chem_id = {} and regn.names in {}""".format(chem_id, tuple(annotations))
            
            conn = self.conn
        
        elif self.db_tag == 'cr':
            query_ = """SELECT synonym.name as reg_number, source.name as source_name, 
                        subs_ann.original_annotation, annotation.general, annotation.category
                        FROM substance sub
                        left join synonym on synonym.subsid = sub.id
                        left join source on source.id = sub.sourceid
                        left join subs_ann on subs_ann.subsid = sub.id
                        left join annotation on annotation.id = subs_ann.annid
                        where synonym.name = '{}'
                        and subs_ann.original_annotation in {}""".format(chem_id,tuple(annotations))
            
            conn = self.compounddb_conn
        
        chem_id_annotations = pd.read_sql_query(query_, conn)
        chem_id_annotations.drop_duplicates(inplace=True)
        
        return chem_id_annotations
    
    def check_cr_source(self, endpoint: str, chem_id_annotations: pd.DataFrame) -> str:
        """
            Checks which source the annotation comes from in CR database.

            :param chem_id_annotations:

            :return final_annotation:
        """

        sources = chem_id_annotations['source_name'].drop_duplicates()

        # List of sources that should be checked. 
        # Decision taken:
        # - REACH Annex VI is considered as yes_source since is the one including Harmonised entries for CLP (Harmonised C&L)
        # - EPA Genetox is considered yes_sources since it's a compendium of 3000 substances analyzed for genetic toxicology, currently outdated
        # - CosmosDB, EFSA OpenFoodTox and Inditex related sources aren't considered since they don't have GHS hazard annotations

        not_considered = ['CosmosDB', 'EFSA OpenFoodTox', 'Inditex',
                          'Inditex MRSL for Parts Suppliers',
                          'Inditex MRSL for Wet Processing Units']

        yes_sources = ['SVHC','REACH Annex VI','EPA Genetox']
        pending_sources = ['CLP Notification', 'REACH Registration','REACH Annex III']

        yes_ann = self.check_yes_or_no(sources_df=sources, cr_source=yes_sources)
        if yes_ann:
            final_annotation = yes_ann
        else:
            pending_ann = self.check_pending(sources_df=sources,cr_source=pending_sources)
            final_annotation = pending_ann

        return final_annotation

    def check_source_of_annotation(self, endpoint: str, chem_id_annotations: pd.DataFrame) -> str:
        """
            Checks which source the annotation comes from and assign either a YES or a Pending annotation
            for the endpoint of interest.

            :param chem_id_annotations:

            :return final_annotation
        """
        
        sources = chem_id_annotations[['general_regulation_name','specific_regulation_name','subspecific_regulation_name',
        'special_cases_name','additional_information_name','names']].drop_duplicates()
       
        # We use this lists to check the presence of annotations in these regulations,
        # which are the ones that are used in the USC Workflow. 
        
        gen_regs = ['clp', 'pbt_vpvb', 'endocrine_disruptors','GHS']
        spec_regs = ['svhc', 'harmonised_C&L','annex_vi']
        subspec_regs = ['candidate_list','hazard_class']

        # Registration dossiers and notification. To decide whether to add Pending or not
        reg_dos_not = ['registration_dossier','notification']
        drafts = ['Submitted SVHC intentions','Withdrawn SVHC intentions and submissions','Amendment 2016/1179']

        # These list include terms that indicates a negative endpoint (NO)
        Negative_endpoint = ['Not PBT', 'Not vPvB']
            
        yes_or_no_ann, source, hazard = self.check_yes_or_no(sources_df=sources, general_regs=gen_regs, specific_regs=spec_regs, subspec_regs=subspec_regs,
        spec_cases=reg_dos_not, drafts=drafts, neg_ans=Negative_endpoint)
        if yes_or_no_ann:
            final_annotation = yes_or_no_ann
        else:
            pending_ann, source, hazard = self.check_pending(sources_df=sources, spec_cases=reg_dos_not, drafts=drafts)
            final_annotation = pending_ann

        return final_annotation, source, hazard
    
    def check_yes_or_no(self, sources_df: pd.DataFrame, cr_source: list = None, general_regs: list = None, 
                    specific_regs: list = None, subspec_regs: list = None, spec_cases: list = None, 
                    drafts: list = None, neg_ans: list = None) -> Optional[str]:
        """
            Checks regulations to get a YES

            :param sources_df:
            :param general_regs: general regulations list
            :param specific_regs: specific regulations list
            :param no_presence: annotations that indicate no presence of hazard annotation

            :return annotation:
        """

        if self.db_tag == 'cii':
            yes_df = sources_df[(((sources_df['general_regulation_name'].isin(general_regs)) |
                    (sources_df['specific_regulation_name'].isin(specific_regs)) |
                    (sources_df['subspecific_regulation_name'].isin(subspec_regs))) &
                    ~(sources_df['special_cases_name'].isin(spec_cases)) &
                    ~(sources_df['additional_information_name'].isin(drafts)))]
            
            no_df = sources_df[sources_df['names'].isin(neg_ans)]
            
        elif self.db_tag == 'cr':
            yes_df = sources_df.isin(cr_source)
        
        if not no_df.empty:
            annotation = 'NO'
            source = sources_df.loc[sources_df['names'].isin(neg_ans), 'general_regulation_name'].unique()
            hazard = sources_df.loc[sources_df['names'].isin(neg_ans), 'names'].unique()
        elif not yes_df.empty:
            annotation = 'YES'
            source = sources_df.loc[(((sources_df['general_regulation_name'].isin(general_regs)) |
                    (sources_df['specific_regulation_name'].isin(specific_regs)) |
                    (sources_df['subspecific_regulation_name'].isin(subspec_regs))) &
                    ~(sources_df['special_cases_name'].isin(spec_cases)) &
                    ~(sources_df['additional_information_name'].isin(drafts))), 'general_regulation_name'].unique()
            hazard = sources_df.loc[(((sources_df['general_regulation_name'].isin(general_regs)) |
                    (sources_df['specific_regulation_name'].isin(specific_regs)) |
                    (sources_df['subspecific_regulation_name'].isin(subspec_regs))) &
                    ~(sources_df['special_cases_name'].isin(spec_cases)) &
                    ~(sources_df['additional_information_name'].isin(drafts))), 'names'].unique()
        else:
            annotation = source = hazard = None
        
        return annotation, source, hazard

    def check_pending(self, sources_df: pd.DataFrame, cr_source: list = None, 
                        spec_cases: list = None, drafts: list = None) -> str:
        """
            Checks regulations to get Pending if YES hasn't been annotated previously.

            :param sources_df:
            :param spec_cases: general regulations list
            :param drafts: specific regulations list

            :return annotation:
        """

        if self.db_tag == 'cii':
            pending_df = sources_df[(sources_df['special_cases_name'].isin(spec_cases)) |
                                (sources_df['additional_information_name'].isin(drafts))]
        elif self.db_tag == 'cr':
            pending_df = sources_df.isin(cr_source)
        
        if not pending_df.empty:
            annotation = 'Pending'
            source = sources_df.loc[(sources_df['special_cases_name'].isin(spec_cases)) |
                                (sources_df['additional_information_name'].isin(drafts)), 'special_cases_name'].unique()
            hazard = sources_df.loc[(sources_df['special_cases_name'].isin(spec_cases)) |
                                (sources_df['additional_information_name'].isin(drafts)), 'names'].unique()
        else:
            annotation = 'No information'
            source = hazard = None
        
        return annotation, source, hazard
    
    def add_endpoint_annotations_to_database(self, chemid_endpoint_annotations:pd.DataFrame, table_name: str):
        """
            Adds dataframe information to CII database.

            :param chemid_endpoint_annotations:
        """

        for i, row in chemid_endpoint_annotations.iterrows():
            chem_id = int(row['chem_id'])
            cmr = row['cmr']
            pbt = row['pbt']
            vpvb = row['vpvb']
            endoc = row['endocrine_disruptor']
            c = row['c']
            m = row ['m']
            r = row['r']
            p = row['p']
            b = row['b']
            t = row['t']
            vp = row['vp']
            vb = row['vb']
            ar = row['androgen_rc']
            er = row['estrogen_rc']
            gr = row['glucocorticoid_rc']
            tr = row['thyroid_rc']

            self.add_annotation(chem_id,cmr,pbt,vpvb,endoc,c,m,r,p,b,t,vp,vb,ar,er, gr,tr,table_name)
    
    def add_annotation(self, chem_id: int, cmr: str, pbt: str, vpvb: str, endoc: str,
                        c: str, m: str, r: str, p: str, b: str, t: str, vp: str, vb: str, 
                        ar: str, er: str, gr: str, tr:str, table_name: str):
        """
            Adds new annotation into database.
            Checks if there is already one and updates it if necessary

            :param chem_id:
            :param cmr:
            :param pbt:
            :parab vpvb:
            :param endoc:
        """
        
        ep_cmd = "SELECT * FROM {} WHERE chem_id = '{}';".format(table_name, chem_id)
        ep_list = self.check_presence_or_absence_endpoint(ep_cmd)
        
        if ep_list:
            cmr_db = ep_list[1]
            pbt_db = ep_list[2]
            vpvb_db = ep_list[3]
            endoc_db = ep_list[4]
            c_db = ep_list[7]
            m_db = ep_list[8]
            r_db = ep_list[9]
            p_db = ep_list[10]
            b_db = ep_list[11]
            t_db = ep_list[12]
            vp_db = ep_list[13]
            vb_db = ep_list[14]
            ar_db = ep_list[15]
            er_db = ep_list[16]
            gr_db = ep_list[17]
            tr_db = ep_list[18]
            
            self.update_annotations('cmr',cmr_db,cmr,chem_id,table_name)
            self.update_annotations('pbt',pbt_db,pbt,chem_id,table_name)
            self.update_annotations('vpvb',vpvb_db,vpvb,chem_id,table_name)
            self.update_annotations('endocrine_disruptor',endoc_db,endoc,chem_id,table_name)
            self.update_annotations('c',c_db,c,chem_id,table_name)
            self.update_annotations('m',m_db,m,chem_id,table_name)
            self.update_annotations('r',r_db,r,chem_id,table_name)
            self.update_annotations('p',p_db,p,chem_id,table_name)
            self.update_annotations('b',b_db,b,chem_id,table_name)
            self.update_annotations('t',t_db,t,chem_id,table_name)
            self.update_annotations('vp',vp_db,vp,chem_id,table_name)
            self.update_annotations('vb',vb_db,vb,chem_id,table_name)
            self.update_annotations('androgen_rc',ar_db,ar,chem_id,table_name)
            self.update_annotations('estrogen_rc',er_db,er,chem_id,table_name)
            self.update_annotations('glucocorticoid_rc',gr_db,gr,chem_id,table_name)
            self.update_annotations('thyroid_rc',tr_db,tr,chem_id,table_name)
        else:
            max_id_cmd = """SELECT max(id) from {}""".format(table_name)
            insert_query = """INSERT INTO """+table_name+""" (id, chem_id, cmr, pbt, vpvb, endocrine_disruptor, c,m,r,p,b,t,vp,vb,androgen_rc,estrogen_rc,glucocorticoid_rc,thyroid_rc)
                             VALUES ({},{},'{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}','{}')"""
            
            self.insert_in_database(max_id_cmd,insert_query,chem_id,cmr,pbt,vpvb,endoc,c,m,r,p,b,t,vp,vb,ar,er,gr,tr)

    def update_annotations(self, endpoint: str, old_endpoint_annotation: str, new_endpoint_annotation: str, chem_id: int, table_name: str):
        """
            Updates old annotations with new ones if those are different from each other.

            :param endpoint: endpoint to update (CMR, PBT, vPvB, Endocrine disruptor)
            :param old_endpoint_annotation: annotation in database
            :param new_endpoint_annotation: annotation assigned after algorithm
            :param chem_id: substance id of the compound in the database
        """

        if (old_endpoint_annotation == 'YES') or \
            (old_endpoint_annotation == new_endpoint_annotation) or \
            (old_endpoint_annotation != 'No information' and old_endpoint_annotation is not None and new_endpoint_annotation == 'No information'):
            pass
        else:
            update_query = """UPDATE {} SET {} = '{}' WHERE chem_id = '{}';""".format(table_name, endpoint, new_endpoint_annotation,chem_id)
            self.curs.execute(update_query)
            self.conn.commit()

    def check_presence_or_absence_endpoint(self, query: str) -> list:
        """
            Checks if the desired input is already in DB and returns the list of endpoints

            :param query: query to check

            :return endpoints_:
        """

        self.curs.execute(query)
        try:
            endpoints_ = self.curs.fetchall()[0]
            self.conn.commit()
        except (TypeError, IndexError):
            endpoints_ = None
        
        return endpoints_

    def calculate_gap_filling(self, endpoint: str) -> pd.DataFrame:
        """
            Function that calls gap filling object and returns a dataframe with the calculation.

            :param endpoint: endpoint to calculate.

            :return gap_filling_frame: gap filling dataframe
        """

        endpoint_list = ['cmr', 'pbt', 'vpvb', 'ed']

        if endpoint.lower() not in endpoint_list:
            sys.stderr.write('Wrong endpoint added. Please write one of the following: {}\n'.format(endpoint_list))
            sys.exit(1)
            
        from .gap_filling import GapFillingCalculator as gapfill

        substances_cii = self.get_substances_with_merged_sources_and_exp_annotations()

        model_preds = pd.read_sql_query("""select ci."name", p.value, md.endpoint, md.creation_date 
                                        from predictions p 
                                        join chem_id ci on ci.id = p.chem_id
                                        join model_detail md on md.id = p.modelid 
                                        order by p.modelid, p.chem_id """,self.conn)
        
        Gap_filling = gapfill(substances_cii=substances_cii, model_predictions=model_preds)
        gap_filling_frame = Gap_filling.merge_source_ans_and_preds(endpoint=endpoint)

        return gap_filling_frame