"""
    Created by: Eric March Vila (eric.march@upf.edu)
    On: 14/02/2020, 10:00 AM
"""

import numpy as np
import pandas as pd
import psycopg2

from compoundDB import inputtools as it
from psycopg2.extensions import register_adapter, AsIs, Float
from typing import Union

class Connector():
    """
       Main class to connect to CII and CR database.
    """

    def __init__(self, host: str = None, dbname: str = None, user: str = None, password: str = None):
        """
            Initializes class with main arguments for psycopg2 connection to the database.
            
            :param host:
            :param dbname:
            :param user:
            :param password:
        """

        self.host = host
        self.dbname = dbname
        self.user = user
        self.password = password

        register_adapter(float, self.nan_to_null)

    def open_connection(self) -> psycopg2.extensions.connection:
        """
            Opens a psycopg2 connection to the PostgreSQL DB and returns the connection object.

            :return self.conn:
        """

        self.conn = psycopg2.connect(host=self.host, dbname=self.dbname, user=self.user, password=self.password)
        self.curs = self.conn.cursor()

        return self.conn

    def compoundDB_connection(self) -> psycopg2.extensions.connection:
        """
            Connects to compoundDB (CR)

            :return self.compounddb_conn:
        """

        self.compounddb_conn = it.openconnection(host='172.20.16.74', password='DBAdmin')

    def nan_to_null(self, f, _NULL: str =AsIs('NULL'), _NaN: np.nan =np.NaN, _Float: float =Float) -> Union[float,str]:
        """
            Special function to handle NaN values from numpy when passed as part of a query in postgres.
            NaN should be interpreted as NULL.

            :param _NULL:
            :param _NaN:
            :param _Float:

            :returns _Float(f)
            :returns _NULL

            BUG: not sure if this is working
        """

        if f is not _NaN:
            return _Float(f)
        return _NULL

    #### Substances block
    #### Functions that interact with the substances tables and extract them as dataframes

    def get_substances(self) -> pd.DataFrame:
        """
            Get subsances from database

            :return substances_df
        """

        substances_df = pd.read_sql_query("""SELECT id, chem_id, class_name, preferred_name, mol_formula FROM public.substance;""",self.conn)

        return substances_df
    
    def get_curated_substances(self) -> pd.DataFrame:
        """
            Get curated substances from database

            :return curated_substances_df:
        """

        curated_substances_df = pd.read_sql_query("""SELECT id, chem_id, class_name_curated, preferred_name_curated, mol_formula_curated FROM public.substance;""",self.conn)

        return curated_substances_df
    
    def get_substances_with_cas_number(self) -> pd.DataFrame:
        """
            Get substances with CAS numbers

            :return substances_cas:
        """

        substances_cas = pd.read_sql_query("""SELECT sub.id as substance_id, sub.chem_id as chemical_id, sub.class_name_curated, 
                                            sub.preferred_name_curated, sub.mol_formula_curated, cid."name" as Chemical_identifier, 
                                            ct."type" as Type_of_identifier
                                            FROM substance sub
                                            left join chem_id cid on cid.id = sub.chem_id
                                            left join chem_type ct on ct.id = cid.chem_type_id
                                            where cid.chem_type_id = 1""", self.conn)

        return substances_cas

    def get_substances_with_chemical_identifiers(self) -> pd.DataFrame:
        """
            Get substances with chemical identifiers (CAS, EC, Index) from database

            :return substances_chem_id:
        """
        
        substances_chem_id = pd.read_sql_query("""SELECT sub.id as substance_id, sub.chem_id as chemical_id, sub.class_name_curated, sub.preferred_name_curated, 
                                                sub.mol_formula_curated, cid."name" as Chemical_identifier, ct."type" as Type_of_identifier
                                                FROM substance sub
                                                left join chem_id cid on cid.id = sub.chem_id
                                                left join chem_type ct on ct.id = cid.chem_type_id""",self.conn)

        return substances_chem_id
    
    def get_substances_with_structure(self) -> pd.DataFrame:
        """
            Get substances with SMILES from database

            :return substance_structures:
        """

        substance_structures = pd.read_sql_query("""SELECT s.class_name_curated, s.preferred_name_curated, cid."name" , struc.chem_id, struc."structure",
                                                struc.structure_curated, struc.substance_type_id, st.type
                                                FROM substance_structure struc
                                                left join chem_id cid on cid.id = struc.chem_id 
                                                left join substance s on s.chem_id = struc.chem_id 
                                                left join substance_type st on st.id = struc.substance_type_id
                                                order by struc.subs_id ASC""", self.conn)

        return substance_structures
    
    def get_endpoints(self) -> pd.DataFrame:
        """
            Gets all the content from endpoint_annotation table

            :return endpoint_annotation:
        """

        endpoint_annotation = pd.read_sql_query("""SELECT *
                                                FROM endpoint_annotation  
                                                order by id ASC""", self.conn)
        
        return endpoint_annotation
    
    def get_endpoints_with_chemical_id(self) -> pd.DataFrame:
        """
            Get all the substances with endpoint annotations and the chemical id
            corresponding to each.

            :return endpoint_chem_id_annotation:
        """

        endpoint_chem_id_annotation = pd.read_sql_query("""SELECT s.class_name_curated, s.preferred_name_curated, ci."name",
                                                        ea.cmr, ea.pbt, ea.vpvb, ea.endocrine_disruptor, ea.c, ea.m, ea.r,
                                                        ea.p, ea.b, ea.t, ea.vp, ea.vb, ep.androgen_rc, ep.estrogen_rc, ep.glucocorticoid_rc, ep.thyroid_rc
                                                        FROM endpoint_annotation ea
                                                        left join chem_id ci on ci.id = ea.chem_id 
                                                        left join substance s on s.chem_id = ea.chem_id
                                                        order by s.chem_id ASC """, self.conn)
        
        return endpoint_chem_id_annotation

    def get_substances_with_endpoint_annotations_and_structure(self) -> pd.DataFrame:
        """
            Get substances with SMILES and endpoint annotations. The aim is to generate an sdf from 
            the resulting dataframe

            :return sub_ann_struc:
        """
        
        sub_ann_struc = pd.read_sql_query("""SELECT sub.class_name_curated, sub.preferred_name_curated, cid.id as chem_id, cid."name", sub.mol_formula_curated,
                                                str.structure, str.structure_curated, st.type, ep.cmr, ep.pbt, ep.vpvb, ep.endocrine_disruptor, ep.c, ep.m, 
                                                ep.r, ep.p, ep.b, ep.t, ep.vp, ep.vb, ep.androgen_rc, ep.estrogen_rc, ep.glucocorticoid_rc, ep.thyroid_rc,
                                                ep.ppar,ep.ahr,ep.rxr,ep.rora,ep.nr3c2,ep.nr1h3,ep.pxr
                                                FROM substance sub
                                                left join chem_id cid on sub.chem_id = cid.id 
                                                left join substance_structure str on str.chem_id = cid.id
                                                left join endpoint_annotation ep on ep.chem_id = cid.id
                                                left join substance_type st on st.id = str.substance_type_id
                                                where str.structure_curated is not null 
                                                and cid.chem_type_id = 1
                                                order by sub.chem_id asc;""", self.conn)
        
        return sub_ann_struc

    def get_substances_with_experimental_endpoint_annotation_and_structure(self) -> pd.DataFrame:
        """
            Returns a dataframe with the experimental endpoint annotations and the SMILES of each compound.

            :return exp_endpoint_ann:
        """

        exp_endpoint_ann = pd.read_sql_query("""SELECT sub.class_name_curated, sub.preferred_name_curated, cid.id as chem_id, cid."name", sub.mol_formula_curated,
                                                str.structure, str.structure_curated, st.type, ep.cmr, ep.pbt, ep.vpvb, ep.endocrine_disruptor, ep.c, ep.m, 
                                                ep.r, ep.p, ep.b, ep.t, ep.vp, ep.vb, ep.androgen_rc, ep.estrogen_rc, ep.glucocorticoid_rc, ep.thyroid_rc,
                                                ep.ppar,ep.ahr,ep.rxr,ep.rora,ep.nr3c2,ep.nr1h3,ep.pxr
                                                FROM experimental_endpoint_annotation ep
                                                left join substance sub on sub.chem_id = ep.chem_id 
                                                left join chem_id cid on ep.chem_id = cid.id 
                                                left join substance_structure str on str.chem_id = cid.id
                                                left join substance_type st on st.id = str.substance_type_id
                                                where str.structure_curated is not null
                                                and cid.chem_type_id = 1  
                                                order by ep.chem_id asc;""", self.conn)

        return exp_endpoint_ann

    def get_substances_with_predicted_endpoint_annotations_and_structure(self) -> pd.DataFrame:
        """
            Returns a dataframe with the predicted endpoint annotations and the SMILES of each compound

            :return predicted_endpoint_ann:
        """

        pass
    
    #### Chemical Identifier block (CAS/EC/Index)
    #### Functions that interact with chem id tables and extract them as dataframes

    def get_chem_id_dataframe(self) -> pd.DataFrame:
        """
            Get chem id dataframe from database

            :return chem_id_df
        """
        
        chem_id_df = pd.read_sql_query("""SELECT id, "name", chem_type_id, subs_id FROM public.chem_id;""",self.conn)

        return chem_id_df

    def get_chem_id_type_dataframe(self) -> pd.DataFrame:
        """
            Get chemical identifier type from database

            :return chem_type_df:
        """
        
        chem_type_df = pd.read_sql_query("""SELECT id, "type" FROM public.chem_type;""",self.conn)

        return chem_type_df

    #### Regulations block
    #### Functions that interact with the regulations tables and extract them as dataframes

    def get_regulation_country(self) -> pd.DataFrame:
        """
            Get the regulation country from the database

            :return regulation_country:
        """

        regulation_country = pd.read_sql_query("""SELECT id, country
                                                FROM public.regulation_country;""", self.conn)
        
        return regulation_country

    def get_regulation_type(self) -> pd.DataFrame:
        """
            Get the regulation type from the database

            :return regulation_type:
        """

        regulation_type = pd.read_sql_query("""SELECT id, "type"
                                            FROM public.regulation_type;""", self.conn)

        return regulation_type

    def get_general_regulations(self) -> pd.DataFrame:
        """
            Get the general regulations from the database

            :return general_regulations:
        """

        general_regulations = pd.read_sql_query("""SELECT id, general_regulation_name
                                                FROM public.general_regulation;""", self.conn)
        
        return general_regulations
    
    def get_specific_regulations(self) -> pd.DataFrame:
        """
            Get the specific regulations from the database

            :return specific_regulations:
        """

        specific_regulations = pd.read_sql_query("""SELECT id, specific_regulation_name
                                                FROM public.specific_regulation;""", self.conn)
        
        return specific_regulations

    def get_subspecific_regulations(self) -> pd.DataFrame:
        """
            Get subspecific regulations from the database

            :return subspecific_regulations:
        """

        subspecific_regulations = pd.read_sql_query("""SELECT id, subspecific_regulation_name
                                                    FROM public.subspecific_regulation;""",self.conn)
        
        return subspecific_regulations

    def get_special_cases(self) -> pd.DataFrame:
        """
            Get special cases names from the database

            :return special_cases:
        """

        special_cases = pd.read_sql_query("""SELECT id, special_cases_name
                                            FROM public.special_cases_regulation;""",self.conn)
        
        return special_cases

    def get_additional_information(self) -> pd.DataFrame:
        """
            Get additional information names from the database

            :return additional_information:
        """

        additional_information = pd.read_sql_query("""SELECT id, additional_information_name
                                                    FROM public.additional_information_regulation;""",self.conn)
        
        return additional_information
    
    def get_regulation_names(self) -> pd.DataFrame:
        """
            Get regulation names from the database

            :return regulation_names:
        """

        regulation_names = pd.read_sql_query("""SELECT id, names
                                            FROM public.regulation_names;""",self.conn)

        return regulation_names
    
    def get_big_regulations_table(self) -> pd.DataFrame:
        """
            Gets the big regulations table wiht only the id's

            :return big_regulations_table:
        """
        
        big_regulations_table = pd.read_sql_query("""SELECT id, chem_id, reg_country_id, reg_type_id, gen_reg_id, 
                                                    spec_reg_id, subspec_reg_id, special_cases_id, additional_information_id, 
                                                    chem_id_name, chem_type_id, regulation_id FROM public.regulations;""",self.conn)
        
        return big_regulations_table
        
    def get_regulations_per_substance(self) -> pd.DataFrame:
        """
            Get regulations for each substance from the database

            :return regulations_per_substance:
        """
        
        regulations_per_substance = pd.read_sql_query("""SELECT reg.id, cid."name" as chemical_identifier, sub.class_name_curated , 
                                                        sub.preferred_name_curated , rco.country, rt."type", rg.general_regulation_name, 
                                                        rspec.specific_regulation_name, rsub.subspecific_regulation_name, 
                                                        rsc.special_cases_name, addr.additional_information_name, regn.names
                                                        FROM regulations reg
                                                        LEFT JOIN chem_id cid ON cid.id = reg.chem_id
                                                        LEFT JOIN substance sub ON sub.chem_id = cid.id
                                                        left join regulation_country rco on rco.id = reg.reg_country_id
                                                        left join regulation_type rt on rt.id = reg.reg_type_id
                                                        left join general_regulation rg on rg.id = reg.gen_reg_id
                                                        left join specific_regulation rspec on rspec.id = reg.spec_reg_id
                                                        LEFT JOIN subspecific_regulation rsub ON rsub.id = reg.subspec_reg_id
                                                        left join special_cases_regulation rsc on rsc.id = reg.special_cases_id
                                                        left join additional_information_regulation addr on addr.id = reg.additional_information_id
                                                        LEFT JOIN chem_type ct ON ct.id = cid.chem_type_id
                                                        LEFT JOIN regulation_names regn ON regn.id = reg.regulation_id
                                                        order by reg.chem_id asc;""",self.conn)

        return regulations_per_substance
    
    def get_regulations_with_id(self) -> pd.DataFrame:
        """
            Get regulations with id from the database

            :return regulations_with_id:
        """
        
        regulations_with_id = pd.read_sql_query("""SELECT DISTINCT rco.country, rco.id as "country_id", rt."type", rt.id as "reg_type_id",
                                                    rg.general_regulation_name, rg.id as "general_reg_id", rspec.specific_regulation_name,
                                                    rspec.id as "spec_reg_id", rsub.subspecific_regulation_name, rsub.id as "subspec_reg_id",
                                                    rsc.special_cases_name, rsc.id as "special_cases_id"
                                                    FROM regulations reg
                                                    left join regulation_country rco on rco.id = reg.reg_country_id
                                                    left join regulation_type rt on rt.id = reg.reg_type_id
                                                    left join general_regulation rg on rg.id = reg.gen_reg_id
                                                    left join specific_regulation rspec on rspec.id = reg.spec_reg_id
                                                    LEFT JOIN subspecific_regulation rsub ON rsub.id = reg.subspec_reg_id
                                                    left join special_cases_regulation rsc on rsc.id = reg.special_cases_id""",self.conn)

        return regulations_with_id
    
    def get_substances_with_merged_sources_and_exp_annotations(self, remove_class: bool = False) -> pd.DataFrame:
        """
            Merges the experimental annotations with the source's ones to have a proper dataframe containing all the information
            regarding endpoint annotations
            
            :param remove_class: boolean. If True, removes class name and only takes into account preferred name curated, which is changed into name.
            :return merged_dataframe:
        """

        sources = self.get_substances_with_endpoint_annotations_and_structure()
        experimental = self.get_substances_with_experimental_endpoint_annotation_and_structure()
        
        sources.rename(columns={'name':'CAS'}, inplace=True)
        experimental.rename(columns={'name':'CAS'}, inplace=True)


        # endpoint_list = ['cmr', 'pbt', 'vpvb', 'endocrine_disruptor', 
        # 'c', 'm', 'r', 'p', 'b', 't', 'vp', 'vb', 'androgen_rc', 'estrogen_rc', 'glucocorticoid_rc', 'thyroid_rc']

        endpoint_list = ['cmr', 'pbt', 'vpvb', 'endocrine_disruptor', 
        'c', 'm', 'r', 'p', 'b', 't', 'vp', 'vb', 'androgen_rc', 'estrogen_rc', 'glucocorticoid_rc', 'thyroid_rc',
        'ppar','ahr','rxr','rora','nr3c2','nr1h3','pxr']

        for i, row in sources.iterrows():
            source_cas = row['CAS']
            if source_cas in experimental.CAS.values:
                exp_ans = experimental.loc[experimental['CAS'] == source_cas, endpoint_list]
                exp_ans_dict = self.proper_dict(exp_ans)
                for key,value in exp_ans_dict.items():
                    if sources.loc[sources.index == i, key].values != 'YES':
                        sources.loc[sources.index==i, key] = value

        if remove_class:
            sources.loc[sources['preferred_name_curated'].isna(),'preferred_name_curated'] = sources.loc[sources['preferred_name_curated'].isna(),'class_name_curated']
            sources.drop('class_name_curated', axis=1, inplace=True)
            sources.rename(columns={'preferred_name_curated':'name'},inplace = True)

        return sources
    
    def proper_dict(self, experimental_dataframe: pd.DataFrame) -> dict:
        """
            Converts the experimental dataframe into a dictionary of annotations and values

            :param experimental_dataframe: dataframe with experimental endpoint annotations

            :return proper_dict: dictionary with the number of experimental annotations per endpoint
        """
        
        dict_ = experimental_dataframe.to_dict('records')
        proper_dict = {}
        for key, value in dict_[0].items():
            if value == 'No information' or value is None:
                pass
            else:
                proper_dict.update({key:value})
        
        return proper_dict

    # List of regulation dataframes to iterate over

    def get_dict_of_regulation_dataframes(self) -> dict:
        """
            Puts in a list the main regulation dataframe so one can iterate over them

            :return self.df_dict:
        """

        gen_reg = self.get_general_regulations()
        spec_reg = self.get_specific_regulations()
        subspec_reg = self.get_subspecific_regulations()
        special_reg = self.get_special_cases()

        self.df_dict = {'general_regulation':gen_reg, 'specific_regulation':spec_reg, 
                        'subspecific_regulation':subspec_reg, 'special_cases_names':special_reg}
        
        return self.df_dict