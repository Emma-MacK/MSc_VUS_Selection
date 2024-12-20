import sys
import numpy as np
import functools
import labkey
import pandas as pd
import matplotlib as mpl
from collections import Counter

print("searching for gene: ", sys.argv[1])

# set up labkey query

def labkey_to_df(sql_query, database, maxrows):
    ver = labkey.__version__

    if ver == '1.2.0' or ver == '1.4.0' or ver == '1.4.1':
        server_context = labkey.utils.create_server_context(
            domain = "labkey-embassy.gel.zone",
            container_path = database,
            context_path = "labkey",
            use_ssl = True
        )

        results =  labkey.query.execute_sql(
            server_context,
            schema_name = "lists",
            sql = sql_query,
            max_rows = maxrows
        )

    if ver == '2.4.0':
        from labkey.api_wrapper import APIWrapper

        labkey_server = "labkey-embassy.gel.zone"
        project_name = database
        contextPath = "labkey"
        schema = 'lists'
        api = APIWrapper(
            labkey_server,
            project_name,
            contextPath,
            use_ssl=True
        )
        results = api.query.execute_sql(sql=sql_query, schema_name=schema, max_rows = maxrows)
    return(pd.DataFrame(results['rows']))

version = "/main-programme/main-programme_v18_2023-12-21"

# get gene of intests
gene_name = sys.argv[1]

# get information for variants occuring in  gene of interest
exomiser_hgvs_sql = (f'''
    SELECT participant_id, genomic_feature_hgnc, sift, score, hgvs, poly_phen, genotype, mode_of_inheritance, chromosome, position, reference, alternate
    FROM exomiser
    WHERE (genomic_feature_hgnc like '{gene_name}' ) AND (assembly = 'GRCh38') AND (participant_type = 'Proband') AND participant_id NOT IN (SELECT participant_id FROM submitted_diagnostic_discovery) AND participant_id NOT IN (SELECT participant_id FROM gmc_exit_questionnaire WHERE case_solved_family LIKE 'yes') AND participant_id IN (SELECT participant_id FROM participant WHERE programme_consent_status = 'Consenting') 
''')

# filter to only get probands
# only get build 38 variants
# participant_id NOT IN (SELECT participant_id FROM submitted_diagnostic_discovery) - removes cases that have been solved or possibly solved 
# AND participant_id IN (SELECT participant_id FROM programme_consent_status WHERE programme_consent_status = 'Consenting')  - only get info from patients where there is consent 

exomiser_hgvs_query = labkey_to_df(exomiser_hgvs_sql, version, 100000)

# check dataframe is not empty
if exomiser_hgvs_query.empty:
    print('DataFrame is empty!')


# get participant ids from table to find panels applied
participants = list(exomiser_hgvs_query['participant_id'].astype(str))
participants_string = ','.join(participants)

exomiser_panels = (f'''
    SELECT participant_id, panel_name, panel_version
    FROM panels_applied
    WHERE participant_id IN ({participants_string})
''')

panel_applied_to_participants = labkey_to_df(exomiser_panels, version, 100000)

# combine panel name and version
panel_applied_to_participants['panel_name_version'] = panel_applied_to_participants['panel_name'] + "-" + panel_applied_to_participants['panel_version'].astype(str)
panel_applied_to_participants = panel_applied_to_participants.drop(['panel_name', 'panel_version'], axis=1)
panel_applied_to_participants['participant_id'] = panel_applied_to_participants['participant_id'].astype(str)

# create vcf file
vcf_file_name = gene_name + "_vep_input.txt"
vcf_builder_df = exomiser_hgvs_query[['chromosome', 'position', 'reference', 'alternate', 'hgvs']]

# remove duplicate lines
vcf_no_duplicates = vcf_builder_df.drop_duplicates()

def print_row(row):
    vcf_line = str(row['chromosome']) + "\t" + str(row['position']) + "\t.\t" + str(row['reference']) + "\t" + str(row['alternate']) + "\t.\t.\t" + str(row['hgvs']) + "\n"
    with open(vcf_file_name, "a") as vcf_out:
       vcf_out.write(vcf_line)

# create cvs file
vcf_no_duplicates.apply(print_row, axis=1)

# sort by exomiser score
# this is unecessary as sort later
exomiser_hgvs_query_sorted = exomiser_hgvs_query.sort_values(by=['score'], ascending=False)

# drop rows with NaN exomiser

exomiser_hgvs_query_sorted = exomiser_hgvs_query_sorted.dropna(subset=['score'])
exomiser_hgvs_query_sorted['participant_id'] = exomiser_hgvs_query_sorted['participant_id'].astype(str)

# add panel name_name, pane version

exomiser_hgvs_query_sorted = pd.merge(exomiser_hgvs_query_sorted, panel_applied_to_participants, on='participant_id', how='left')
exomiser_hgvs_query_sorted[['panel_name_version']] = exomiser_hgvs_query_sorted[['panel_name_version']].fillna('No Panel Found')

# print(exomiser_hgvs_query_sorted)
# sort by exomiser_score, so highest score at the top
# dropping duplicate hgvs removes all after the first occurance (which is the occurrence with the highest exomiser score
exomiser_table = exomiser_hgvs_query_sorted.sort_values('score', ascending=False).drop_duplicates(['hgvs'])
exomiser_table = exomiser_table.reindex(columns=['hgvs', 'participant_id', 'genomic_feature_hgnc', 'sift', 'score', 'poly_phen',
                 'genotype', 'mode_of_inheritance',
                 'chromosome', 'position', 'reference', 'alternate', 'panel_name_version'])


# save table to csv
output_table_name = gene_name + "_output_table.csv"
exomiser_table.to_csv(output_table_name)
