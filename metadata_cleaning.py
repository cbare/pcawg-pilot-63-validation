

# df = query2df(syn.chunkedQuery('select * from file where benefactorId=="%s"' % DATA_SOURCE_PROJECT_ID))

annotation_keys = ['center', 'center_name', 'workflow_name', 'variant_workflow',
                   'call_type', 'dataType', 'dataSubType', 'variant_type', 'fileType']

for df in dfs:
    series = df.apply(lambda x: Counter(x).most_common())
    print
    print series['source']

    for key in annotation_keys:
        if key in series:
            print key, '=', series[key]


pd.crosstab(df_all['center'].replace(float('NaN'), 'NA'), df_all['center_name'].replace(float('NaN'), 'NA'))

pd.crosstab(df_all['center'].replace(float('NaN'), 'NA'), df_all['workflow_name'].replace(float('NaN'), 'NA'))

pd.crosstab(df_all['center'].replace(float('NaN'), 'NA'), df_all['call_type'].replace(float('NaN'), 'NA'))
pd.crosstab(df_all['center'].replace(float('NaN'), 'NA'), df_all['dataSubType'].replace(float('NaN'), 'NA'))




def find_sample_id_source(sample_id):
    results = []
    for col in sample_df.columns:
        if sample_id in sample_df[col].values:
            results.append(col)
    return tuple(results)


results = syn.tableQuery('select * from syn2887105')
sample_df = results.asDataFrame()

## figure out which centers named their samples with which IDs
sample_id_sources = df_all.sample_id.apply(find_sample_id_source)
pd.crosstab(sample_id_sources, df_all.source)

# In [254]: pd.crosstab(sample_id_sources, df_all.source)
# Out[254]: 
# source                                   BSC  Broad  DKFZ  EMBL  MDA_HGSC  \
# sample_id                                                                   
# ()                                         0      0     0     0         0   
# (Donor ID,)                                0      0     0     0         0   
# (Normal Analysis ID,)                      0      0     0     0         0   
# (Tumour Analysis ID,)                    567    998  1068   252       108   
# (Tumour Analyzed Sample/Aliquot GUUID,)    0      0     0     0         0   

# source                                   MDA_KChen  McGill  OICR  SANGER  SFU  \
# sample_id                                                                       
# ()                                               0       0     0    9024    0   
# (Donor ID,)                                      0       0     0       0    0   
# (Normal Analysis ID,)                            0       0     0       0    0   
# (Tumour Analysis ID,)                          216      78   252       0  126   
# (Tumour Analyzed Sample/Aliquot GUUID,)          0       0     0     896    0   

# source                                   UCSC  WUSTL  Yale  
# sample_id                                                   
# ()                                          0      6     0  
# (Donor ID,)                               165      0     0  
# (Normal Analysis ID,)                       0      0    59  
# (Tumour Analysis ID,)                       0    358     0  
# (Tumour Analyzed Sample/Aliquot GUUID,)     0      0     0  


## fix UCSC sample_id -> donor_id and Tumour Analysis ID -> analysis_id_tumor
results = clean_query_results(syn.chunkedQuery('select id, name from file where parentId=="syn3107237"'))
for result in results:
    print result['id'], result['name']
    e = syn.get(result['id'], downloadFile=False)
    a = e.annotations
    a['donor_id'] = a['sample_id']
    a['analysis_id_tumor'] = sample_df['Tumour Analysis ID'][ sample_df['Donor ID']==e.sample_id[0] ][0]
    print syn.setAnnotations(e, a)


# source = Source("UCSC", "syn3107237")
# results = syn.chunkedQuery('select * from file where parentId=="%s"' % source.folder_id)
# df = query2df(results)

results = clean_query_results(syn.chunkedQuery('select id, name from file where parentId=="syn3107237"'))
for result in results:
    m = re.match(r"([a-zA-Z0-9_\-]*)\.(gatk_muse_0.9.9.5|gatk_mutect|muse_0.9.9.5)\.(snv_mnv)\.(vcf.gz)", result['name'])
    if m:
        e = syn.get(result['id'], downloadFile=False)
        a = e.annotations

        a['donor_id'] = m.group(1)
        a['workflow_name'] = m.group(2)
        a['call_type'] = 'somatic'
        a['analysis_id_tumor'] = sample_df['Tumour Analysis ID'][ sample_df['Donor ID']==a['donor_id'] ][0]
        a['dataSubType'] = m.group(3)
        a['dataType'] = 'DNA'
        a['disease']  = 'Cancer'
        a['fileType'] = 'vcf'
        a['center'] = 'ucsc'

        print syn.setAnnotations(e, a)
    else:
        raise Exception("Couldn't parse filename: "+result['name'])


## fix sample_id -> analysis_id_tumor and Donor ID -> donor_id
for source in sources:
    if source.name in ["BSC", "Broad", "DKFZ", "EMBL", "MDA_HGSC", "MDA_KChen", "McGill", "OICR", "SFU", "WUSTL"]:
        results = clean_query_results(syn.chunkedQuery('select id, name from file where parentId=="%s"' % source.folder_id))
        for result in results:
            if (result['id'] not in already_done) and (result['id'] not in second_batch_done):
                print result['id'], result['name']
                e = syn.get(result['id'], downloadFile=False)
                a = e.annotations
                a['analysis_id_tumor'] = a['sample_id']
                donor_id_series = sample_df['Donor ID'][ sample_df['Tumour Analysis ID']==e.sample_id[0] ]
                if len(donor_id_series) == 0:
                    print "no donor ID associated with sample_id: ", a['sample_id']
                else:
                    a['donor_id'] = donor_id_series[0]
                print syn.setAnnotations(e, a)


for source in sources:
    if source.name in ["Yale"]:
        results = clean_query_results(syn.chunkedQuery('select id, name from file where parentId=="%s"' % source.folder_id))
        for result in results:
            if (result['id'] not in already_done) and (result['id'] not in second_batch_done):
                print result['id'], result['name']
                e = syn.get(result['id'], downloadFile=False)
                a = e.annotations
                a['analysis_id_tumor'] = a['sample_id']
                if 'donor_id' not in a:
                    donor_id_series = sample_df['Donor ID'][ sample_df['Normal Analysis ID']==e.sample_id[0] ]
                    if len(donor_id_series) == 0:
                        print "no donor ID associated with sample_id: ", a['sample_id']
                    else:
                        a['donor_id'] = donor_id_series[0]
                print syn.setAnnotations(e, a)


