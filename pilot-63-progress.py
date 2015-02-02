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


syn = None

## Project: ICGC-TCGA Whole Genome Pan-Cancer Analysis
## /Files/VCF files/  syn2351328/syn2882202
class Source(object):
    def __init__(self, name, workflow, folder_id, prefix):
        self.name = name
        self.workflow = workflow
        self.folder_id = folder_id
        self.prefix = prefix
    def __repr__(self):
        return "Source(name=\"{name}\", workflow=\"{workflow}\", folder_id=\"{folder_id}\", prefix=\"{prefix}\")".format(**self.__dict__)

sources = [
    Source("DKFZ",   "DKFZ",        "syn3104289", ".somatic"),
    Source("EMBL",   "delly",       "syn3153529", ".somatic"),
    Source("SANGER", "svcp",        "syn3155834", ".somatic"),
    Source("UCSC",   "gatk-muse",   "syn3107237", ".gatk_muse_0.9.9.5"),
    Source("UCSC",   "gatk-mutect", "syn3107237", ".gatk_mutect"),
    Source("UCSC",   "muse",        "syn3107237", ".muse_0.9.9.5")]

PROJECT_SYNAPSE_ID     = "syn2875157"
THIS_SCRIPT_SYNAPSE_ID = "syn3159349"
BAR_CHART_SYNAPSE_ID   = "syn3158825"
TABLE_SYNAPSE_ID       = "syn3159075"

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
    return re.sub(r'(entity\.|file\.|folder\.|project\.|data\.)','', key)

def clean_result(result):
    return {clean_key(key):unlist_singletons(value) for key,value in result.iteritems()}

def clean_query_results(results):
    return [clean_result(result) for result in results]

def remove_properties(result):
    return {key:value for key,value in result.iteritems() if key not in property_keys}


def extract_uuid(name):
    ## example: 0e90fb64-00b2-4b53-bbc7-df8182b84060
    m = re.match("([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})", name, re.IGNORECASE)
    if m:
        return m.group(1)
    else:
        return None


def get_annotations(folder_id, prefix=".somatic"):
    synapse_ids   = []
    sample_ids    = []
    variant_types = []
    results = syn.chunkedQuery('select * from entity where parentId=="%s" order by name' % folder_id)
    for result in clean_query_results(results):

        ## assign variant type or give up
        if result['name'].endswith(prefix+".indel.vcf.gz"):
            variant_type = "indel"
        elif result['name'].endswith(prefix+".snv_mnv.vcf.gz"):
            variant_type = "snv"
        elif result['name'].endswith(prefix+".sv.vcf.gz"):
            variant_type = "sv"
        else:
            continue

        ## assign sample id or give up
        sample_id = extract_uuid(result['name'])
        if 'sample_id' in result:
            if sample_id and sample_id != result['sample_id']:
                raise Exception("Sample ID {sample_id} doesn't match name {name}".format(**result))
            if not sample_id:
                sample_id = result['sample_id']
        if not sample_id:
            continue
        synapse_id = result['id']

        if 'variant_type' in result and variant_type != result['variant_type'].lower():
            raise Exception("Variant type {variant_type} doesn't match name {name}!".format(**result))

        synapse_ids.append(synapse_id)
        sample_ids.append(sample_id)
        variant_types.append(variant_type)
        print synapse_id, sample_id, variant_type

    return pd.DataFrame({
            "synapse_id": synapse_ids,
            "sample_id":  sample_ids,
            "variant_type": variant_types})


def create_metadata_df(sources):

    dfs = []
    for source in sources:
        print "\n\n", "-"*60, "\n", source.name
        df = get_annotations(source.folder_id, source.prefix)
        df['source'] = source.name
        df['workflow'] = source.workflow
        dfs.append(df)

    return pd.concat(dfs)


def plot_progress(df, image_filename):
    width = 0.8
    variant_types = ['snv', 'sv', 'indel']
    palette = {'snv':'#336699', 'sv':'#993333', 'indel':'#669933'}

    counts = []
    colors = []
    labels = []
    variant_type_indexes = {}
    for source in sources:
        for i, variant_type in enumerate(variant_types):
            count = len(df[ (df['source']==source.name) &
                            (df['workflow']==source.workflow) &
                            (df['variant_type']==variant_type) ])
            if count > 0:
                variant_type_indexes[variant_type] = len(counts)
                counts.append(count)
                colors.append(palette[variant_type])
                if source.workflow == source.name:
                    labels.append(source.name)
                else:
                    labels.append(source.workflow + "\n" + source.name)

    ind = np.arange(len(counts))    # the x locations for the bars
    bar = plt.bar(ind, counts, width=width, color=colors, alpha=0.6)

    plt.title('Count of VCF files by source and variant type')
    plt.xticks(ind+width/2., labels)

    ## build a legend with a color patch for each variant type
    plt.legend([bar[variant_type_indexes[v]] for v in variant_types], variant_types)

    ## label each bar with its count
    for i,count in enumerate(counts):
        plt.text(i+width/2, count+1,str(count), horizontalalignment="center")

    fig = plt.gcf()
    fig.set_size_inches(14,9)
    fig.savefig(image_filename, bbox_inches='tight')


def create_table_schema(project, activity):
    cols = [
        Column(name='sample_id',  columnType='STRING', maximumSize=36),
        Column(name='synapse_id', columnType='STRING', maximumSize=13),
        Column(name='variant_type', columnType='STRING', maximumSize=10),
        Column(name='source', columnType='STRING', maximumSize=30),
        Column(name='workflow', columnType='STRING', maximumSize=100)
    ]

    schema = syn.store(
        Schema(name='pilot-63-progress', 
               columns=cols,
               parent=project,),
        activity=activity)
    print "Schema created:", schema.id

    return schema, activity


def add_new_rows_to_table(schema, df):
    results = syn.tableQuery('select synapse_id from %s' % utils.id_of(schema), includeRowIdAndRowVersion=False)
    synapse_ids = [row[0] for row in results]
    df_new_rows = df[ [synapse_id not in synapse_ids for synapse_id in df['synapse_id']] ]
    if df_new_rows.shape[0] > 0:
        return syn.store(Table(schema, df_new_rows))
    else:
        return None

def store_this_script():
    return syn.store(
        File("https://github.com/cbare/pcawg-pilot-63-validation/",
             name="pilot-63-progress.py",
             parentId="syn2875157",
             synapseStore=False),
        used="https://github.com/cbare/pcawg-pilot-63-validation/")

def update_figure_and_table(sources, script_commit_url=None):
    df = create_metadata_df(sources)
    print df.groupby(["source", "variant_type"])["synapse_id"].count()

    schema = syn.get(TABLE_SYNAPSE_ID)
    table = add_new_rows_to_table(schema, df)

    if table:
        script_entity = syn.get(THIS_SCRIPT_SYNAPSE_ID, downloadFile=False)
        if script_commit_url:
            script_entity.externalURL = script_commit_url
            script_entity = syn.store(script_entity)

        activity = Activity(
            name='Pilot-63-progress',
            description='Track VCF files uploaded for the PCAWG Pilot-63 project',
            used=list(set(source.folder_id for source in sources)),
            executed=[script_entity])

        activity = syn.setProvenance(schema, activity)

        image_filename="pilot-63-progress.png"
        plot_progress(df, image_filename)

        bar_chart = syn.get(BAR_CHART_SYNAPSE_ID, downloadFile=False)
        bar_chart.path = "pilot-63-progress.png"
        bar_chart.synapseStore=True
        bar_chart = syn.store(bar_chart, activity=activity)

    return bar_chart


def main():
    global syn

    parser = argparse.ArgumentParser()

    parser.add_argument("-u", "--user", help="UserName", default=None)
    parser.add_argument("-p", "--password", help="Password", default=None)
    parser.add_argument("--debug", help="Show verbose error output from Synapse API calls", action="store_true", default=False)
    parser.add_argument("script_commit_url", metavar="SCRIPT_COMMIT_URL", help="Github URL of this script")

    args = parser.parse_args()

    syn = synapseclient.Synapse(debug=args.debug)
    if not args.user:
        args.user = os.environ.get('SYNAPSE_USER', None)
    if not args.password:
        args.password = os.environ.get('SYNAPSE_PASSWORD', None)
    syn.login(email=args.user, password=args.password)

    update_figure_and_table(sources, args.script_commit_url)



if __name__ == "__main__":
    main()


