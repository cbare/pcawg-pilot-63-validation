from itertools import izip
import argparse
import os
import re
import sys
import synapseclient
import synapseclient.utils as utils
from synapseclient import Project, File, Folder, Activity
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm



## used to remove properties from Synapse query results
property_keys= ['id', 'name', 'description',
                'modifiedOn', 'modifiedByPrincipalId',
                'createdOn', 'createdByPrincipalId',
                'versionLabel', 'versionComment', 'versionNumber',
                'nodeType', 'concreteType',
                'benefactorId', 'parentId', 'eTag']

## remove the annoying prefixes from synapse query results
def unlist_singletons(value):
    if isinstance(value,list) and len(value)==1:
        return value[0]
    else:
        return value

def clean_key(key):
    """remove the prefix that Synapse queries add to field names"""
    prefix, new_key = key.split('.', 1)
    return new_key

def clean_result(result):
    return {clean_key(key):unlist_singletons(value) for key,value in result.iteritems()}

def filter_keys(dictionary, keys_to_remove):
    return {key:value for key,value in dictionary.iteritems() if key not in keys_to_remove}

def clean_query_results(results, keys_to_remove=[]):
    return [filter_keys(clean_result(result), keys_to_remove) for result in results]

def query2df(queryContent, keepSynapseFields=set(['id', 'name'])):
    """Converts the returned query object from Synapse into a Pandas DataFrame
    
    Arguments:
    - `queryContent`: content returned from query or chunkedQuery
    - `keepSynapseFields`: Synapse properties are removed from query results except those named here
    """
    return pd.DataFrame(clean_query_results(queryContent,
        keys_to_remove=set(k for k in property_keys if k not in set(keepSynapseFields))))

def extract_uuid(name):
    ## example: 0e90fb64-00b2-4b53-bbc7-df8182b84060
    m = re.match("([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})", name, re.IGNORECASE)
    if m:
        return m.group(1)
    else:
        return None


DATA_SOURCE_PROJECT_ID = "syn2351328"
PROJECT_SYNAPSE_ID     = "syn2875157"
THIS_SCRIPT_SYNAPSE_ID = "syn3159349"
BAR_CHART_SYNAPSE_ID   = "syn3158825"
TABLE_SYNAPSE_ID       = "syn3159075"


class Source(object):
    def __init__(self, name, folder_id):
        self.name = name
        self.folder_id = folder_id
    def __repr__(self):
        return "Source(name=\"{name}\", folder_id=\"{folder_id}\")".format(**self.__dict__)

sources = [
#          Center        Folder
    Source("Broad",      "syn3165121"),
    Source("BSC",        "syn3165143"),
    Source("DKFZ",       "syn3104289"),
    Source("EMBL",       "syn3153526"),
    Source("McGill",     "syn3165151"),
    Source("MDA_HGSC",   "syn3167886"),
    Source("MDA_KChen",  "syn3165149"),
    Source("OICR",       "syn3167076"),
#    Source("SANGER",     "syn3155834"), ## ???
    Source("SFU",        "syn3165152"),
    Source("UCSC",       "syn3107237"),
    Source("WUSTL",      "syn3165146"),
#    Source("Yale",       "syn3165120")  ## can't tell if these are germline or tumor
    ]


## Once the table is save in Synapse, summary stats can be computed with queries:
## SELECT sample_id, variant_type, count(*) FROM syn3159075 GROUP BY sample_id, variant_type
## SELECT source, workflow, variant_type, count(*) FROM syn3159075 GROUP BY source, workflow, variant_type

def create_metadata_df(sources):
    dfs = []
    for source in sources:
        results = syn.chunkedQuery('select * from file where parentId=="%s"' % source.folder_id)
        df = query2df(results)
        print source.name, df.shape
        df['source'] = source.name
        dfs.append(df)

    df_all = pd.concat(dfs, ignore_index=True)

    ## find duplicates where an equivalent .vcf and .vcf.gz file both exist
    def is_dup(name):
        return name.endswith('.vcf') and (name+".gz" in df_all.name.values)
    dups = df_all.name.apply(is_dup)

    print "%d duplicate .vcf and .vcf.gz files found" % sum(dups)

    ## assume wustle VCFs are somatic variants only ??
    df_all.call_type[ df_all.center=='wustl' ] = 'somatic'
    df_all.workflow_name[ df_all.center=='wustl' ] = 'wustl'

    ## fix DKFZ dataSubType
    df_all.dataSubType[ df_all.center=='DKFZ' ] = df_all.dataType[ df_all.center=='DKFZ' ]

    ## fix EMBL dataSubType and workflow
    df_all.dataSubType[ df_all.center=='EMBL' ] = 'sv'
    df_all.workflow_name[ df_all.center=='EMBL' ] = 'pcawg-delly'

    df_all.dataSubType[ df_all.dataSubType=='indels' ] = 'indel'

    df = df_all.ix[ (~df_all.center.isnull()) & (df_all.fileType=='vcf') & (df_all.call_type=='somatic') & (~dups), :]

    df_progress = df.ix[:,['analysis_id_tumor', 'id', 'dataSubType', 'source', 'workflow_name']]
    df_progress.columns = ['sample_id', 'synapse_id', 'variant_type', 'source', 'workflow']

    return df_all, df, df_progress


def add_new_rows_to_table(df, replace_table=False):
    """Add rows for synapse IDs not already represented in the table or replace the whole table"""
    schema = syn.get(TABLE_SYNAPSE_ID)
    if replace_table:
        ## delete previous entries in pilot-63-progress table
        results = syn.tableQuery('select * from %s' % utils.id_of(schema), resultsAs='rowset')
        syn.delete(results)
    else:
        results = syn.tableQuery('select synapse_id from %s' % utils.id_of(schema), includeRowIdAndRowVersion=False)
        synapse_ids = [row[0] for row in results]
        df = df[ [synapse_id not in synapse_ids for synapse_id in df['synapse_id']] ]

    if df.shape[0] > 0:
        print "Adding %d rows to pilot-63-progress table" % df.shape[0]
        return syn.store(Table(schema, df))
    else:
        print "No new rows for pilot-63-progress table"
        return None


def plot_progress(df, sources, image_filename):
    #plt.close()

    width = 0.8
    padding = 0.2
    variant_types = ['snv_mnv', 'sv', 'indel']
    palette = {'snv_mnv':'#336699', 'sv':'#993333', 'indel':'#669933'}

    counts = []
    colors = []
    labels = []
    variant_type_indexes = {}
    for source in sources:
        workflows = df.workflow[df.source==source.name].unique()
        for i, variant_type in enumerate(variant_types):
            for workflow in workflows:
                count = len(df[ (df['source']==source.name) &
                                (df['workflow']==workflow) &
                                (df['variant_type']==variant_type) ])
                if count > 0:
                    variant_type_indexes[variant_type] = len(counts)
                    counts.append(count)
                    colors.append(palette[variant_type])
                    if workflow == source.name:
                        labels.append(source.name)
                    else:
                        labels.append(source.name + " - " + workflow)

    ind = np.arange(len(counts))    # the x locations for the bars
    bar = plt.bar(ind+padding, counts, width=width, color=colors, alpha=0.6)

    plt.title('Count of VCF files by source and variant type')
    #plt.xticks(ind+width/2., labels)
    for i, label in enumerate(labels):
        plt.text(i+width/2+padding, 3, str(label), rotation='vertical',
            horizontalalignment='center',
            verticalalignment='bottom')

    ## build a legend with a color patch for each variant type
    plt.legend([bar[variant_type_indexes[v]] for v in variant_types], variant_types,
        loc='upper right', bbox_to_anchor=(1.05, 0.6))

    ## label each bar with its count
    for i,count in enumerate(counts):
        plt.text(i+width/2+padding, count+1,str(count), horizontalalignment="center")

    plt.xlim([0,23.2])
    #plt.show()

    fig = plt.gcf()
    fig.set_size_inches(14,9)
    fig.savefig(image_filename, bbox_inches='tight')



def update_figure_and_table(sources, script_commit_url=None, replace_table=False, force_update=False):
    df_all, df, df_progress = create_metadata_df(sources)
    print df_progress.groupby(["source", "variant_type"])["synapse_id"].count()

    table = add_new_rows_to_table(df_progress, replace_table)

    if table or force_update:
        script_entity = syn.get(THIS_SCRIPT_SYNAPSE_ID, downloadFile=False)
        if script_commit_url:
            script_entity.externalURL = script_commit_url
            fileHandle = syn._addURLtoFileHandleService(script_commit_url, mimetype="text/x-python")
            script_entity.dataFileHandleId = fileHandle['id']
            script_entity = syn.store(script_entity)

        activity = Activity(
            name='Pilot-63-progress',
            description='Track VCF files uploaded for the PCAWG Pilot-63 project',
            used=list(set(source.folder_id for source in sources)),
            executed=[script_entity])

        activity = syn.setProvenance(TABLE_SYNAPSE_ID, activity)

        image_filename="pilot-63-progress.png"
        plot_progress(df_progress, sources, image_filename)

        bar_chart = syn.get(BAR_CHART_SYNAPSE_ID, downloadFile=False)
        bar_chart.path = "pilot-63-progress.png"
        bar_chart.synapseStore=True
        bar_chart = syn.store(bar_chart, activity=activity)



def main():
    global syn

    parser = argparse.ArgumentParser()

    parser.add_argument("-u", "--user", help="UserName", default=None)
    parser.add_argument("-p", "--password", help="Password", default=None)
    parser.add_argument("--debug", help="Show verbose error output from Synapse API calls", action="store_true", default=False)
    parser.add_argument("--replace-table", help="Replace whole pilot-63-progress table", action="store_true", default=False)
    parser.add_argument("--force-update", help="Update graph even if there are no new table entries", action="store_true", default=False)
    parser.add_argument("script_commit_url", metavar="SCRIPT_COMMIT_URL", help="Github URL of this script")

    args = parser.parse_args()

    syn = synapseclient.Synapse(debug=args.debug)
    if not args.user:
        args.user = os.environ.get('SYNAPSE_USER', None)
    if not args.password:
        args.password = os.environ.get('SYNAPSE_PASSWORD', None)
    syn.login(email=args.user, password=args.password)

    update_figure_and_table(sources, args.script_commit_url, 
        replace_table=args.replace_table, force_update=args.force_update)



if __name__ == "__main__":
    main()


