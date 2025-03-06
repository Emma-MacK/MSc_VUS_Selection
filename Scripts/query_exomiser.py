import sys
import numpy as np
import functools
import labkey
import pandas as pd

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

version = "/main-programme/main-programme_v19_2024-10-31"

# get gene of intests
gene_name = sys.argv[1]

# get information for variants occuring in  gene of interest
exomiser_hgvs_sql = (f'''
    SELECT participant_id, genomic_feature_hgnc, sift, score, hgvs, poly_phen, genotype, mode_of_inheritance, chromosome, position, reference, alternate
    FROM exomiser
    WHERE (genomic_feature_hgnc like '{gene_name}' ) AND (assembly = 'GRCh38') AND (participant_type = 'Proband') AND participant_id IN (SELECT participant_id FROM participant WHERE programme_consent_status = 'Consenting') 
''')

# filter to only get probands
# only get build 38 variants
# participant_id NOT IN (SELECT participant_id FROM submitted_diagnostic_discovery) - removes cases that have been solved or possibly solved 
# AND participant_id IN (SELECT participant_id FROM programme_consent_status WHERE programme_consent_status = 'Consenting')  - only get info from patients where there is consent 

### TRYING WITH all cases

exomiser_hgvs_query = labkey_to_df(exomiser_hgvs_sql, version, 100000)

# check dataframe is not empty
if exomiser_hgvs_query.empty:
    print('DataFrame is empty!')



hgvcs_part_IDs = exomiser_hgvs_query[['hgvs', 'participant_id']]
hgvcs_part_IDs['participant_id'] = hgvcs_part_IDs['participant_id'].astype(str)

# get participant_panels with exomiser score > 75

row_high_exomiser = exomiser_hgvs_query[exomiser_hgvs_query['score'] > 0.70]
participants_high_exomiser = list(row_high_exomiser['participant_id'].astype(str))
high_participants_string = ','.join(participants_high_exomiser)

exomiser_panels = (f'''
    SELECT participant_id, panel_name, panel_version
    FROM panels_applied
    WHERE participant_id IN ({high_participants_string})
''')

panel_applied_to_high_participants = labkey_to_df(exomiser_panels, version, 100000)

# combine panel name and version
panel_applied_to_high_participants['panels_associated_with_exomiser_scores_over_0.70'] = panel_applied_to_high_participants['panel_name'] + "-" + panel_applied_to_high_participants['panel_version'].astype(str)
panel_applied_to_high_participants['participant_id'] = panel_applied_to_high_participants['participant_id'].astype(str)
high_exomiser_socre_partID_hgvs = pd.merge(panel_applied_to_high_participants, hgvcs_part_IDs, on='participant_id', how='left')

# find top exomiser panel seen in >60% of participants with variant, only done when > 2 participants with variant

# get participant ids from table to find panels applied
participants = list(exomiser_hgvs_query['participant_id'].astype(str))
participants_string = ','.join(participants)

exomiser_panels = (f'''
    SELECT participant_id, panel_name
    FROM panels_applied
    WHERE participant_id IN ({participants_string})
''')

panel_applied_to_participants = labkey_to_df(exomiser_panels, version, 100000)

# combine panel name and version
panel_applied_to_participants['participant_id'] = panel_applied_to_participants['participant_id'].astype(str)
panels_hgvs = pd.merge(panel_applied_to_participants, hgvcs_part_IDs, on='participant_id', how='left')

common_panel_per_hgvs = pd.DataFrame(columns=['hgvs'])

# number of unique participants with that variant
for variant in list(set(panels_hgvs['hgvs'])):
   #print(variant)
   variant_subset = panels_hgvs[panels_hgvs['hgvs'] == variant]
   #print(variant_subset)
   n_participants = len(pd.unique(variant_subset['participant_id']))
   if n_participants <= 1:
      continue
   panel_counts = variant_subset['panel_name'].value_counts().rename_axis("panel_name").to_frame('counts')
   common_panel = panel_counts[panel_counts['counts'] / n_participants >= 0.6]
   if common_panel.empty:
      continue
   # print(common_panel.index.values.tolist())
   freq_panel = common_panel.index.values.tolist()
   high_freq_panel = '; '.join(freq_panel)
   new_row = pd.DataFrame({"hgvs": [variant], "panels_associated_with_at_least_60_percent_of_occurrences": [high_freq_panel]})
   common_panel_per_hgvs = pd.concat([common_panel_per_hgvs, new_row], ignore_index=True)

print("common panels for hgvs")
print(common_panel_per_hgvs)

# get panels where  n_particpance/count >0.6

if not common_panel_per_hgvs.empty:
   agg_functions = {'panels_associated_with_at_least_60_percent_of_occurrences' :  lambda x: ';  '.join(x.unique())}
   common_panel_per_hgvs = common_panel_per_hgvs.groupby('hgvs').agg(agg_functions)


high_exomiser_socre_partID_hgvs = high_exomiser_socre_partID_hgvs.drop(['panel_name', 'panel_version', 'participant_id'], axis=1)
print(high_exomiser_socre_partID_hgvs)

# fing hgvs corresponding to each participant ID
agg_functions = {'panels_associated_with_exomiser_scores_over_0.70' :  lambda x: ';  '.join(x.unique())}
all_panels_high_exomiser_partipants = high_exomiser_socre_partID_hgvs.groupby('hgvs').agg(agg_functions)
 
#print(all_panels_high_exomiser_partipants)

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

#exomiser_hgvs_query_sorted = pd.merge(exomiser_hgvs_query_sorted, panel_applied_to_participants, on='participant_id', how='left')
#exomiser_hgvs_query_sorted[['top_exomiser_scoring_panel']] = exomiser_hgvs_query_sorted[['top_exomiser_scoring_panel']].fillna('No Panel Found')

partic_ID = exomiser_hgvs_query_sorted['participant_id']
participants_file_name = gene_name + "participants_prior_to_filtering.csv"
partic_ID.to_csv(participants_file_name)

# print(exomiser_hgvs_query_sorted)
# sort by exomiser_score, so highest score at the top
# dropping duplicate hgvs removes all after the first occurance (which is the occurrence with the highest exomiser score
exomiser_table = exomiser_hgvs_query_sorted.sort_values('score', ascending=False).drop_duplicates(['hgvs'])


exomiser_table = exomiser_table.reindex(columns=['hgvs', 'participant_id', 'genomic_feature_hgnc', 'sift', 'score', 'poly_phen',
                 'genotype', 'mode_of_inheritance',
                 'chromosome', 'position', 'reference', 'alternate'])


exomiser_table = pd.merge(exomiser_table, all_panels_high_exomiser_partipants, on='hgvs', how='left')
exomiser_table = pd.merge(exomiser_table, common_panel_per_hgvs, on='hgvs', how='left')
# save table to csv
output_table_name = gene_name + "_output_table.csv"
exomiser_table.to_csv(output_table_name)
