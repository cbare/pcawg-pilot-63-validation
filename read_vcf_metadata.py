from itertools import izip
import re
import sv
import synapseclient
import pandas as pd


syn = synapseclient.Synapse()
syn.login()


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

property_keys= ['id', 'name', 'description',
                'modifiedOn', 'modifiedByPrincipalId',
                'createdOn', 'createdByPrincipalId',
                'versionLabel', 'versionComment', 'versionNumber',
                'nodeType', 'concreteType',
                'benefactorId', 'parentId', 'eTag']

def remove_properties(result):
    return {key:value for key,value in result.iteritems() if key not in property_keys}


def extract_uuid(name):
    ## example: 0e90fb64-00b2-4b53-bbc7-df8182b84060
    m = re.match("([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})", name, re.IGNORECASE)
    if m:
        return m.group(1)
    else:
        return None


## Project: ICGC-TCGA Whole Genome Pan-Cancer Analysis 
## /Files/VCF files/  syn2351328/syn2882202
class Source(object):
    def __init__(self, name, folder_id, prefix):
        self.name = name
        self.folder_id = folder_id
        self.prefix = prefix
    def __repr__(self):
        return "Source(name=\"{name}\", folder_id=\"{folder_id}\", prefix=\"{prefix}\")".format(**self.__dict__)

sources = [
    Source("DKFZ", "syn3104289", ".somatic"),
    Source("EMBL", "syn3153529", ".somatic"),
    Source("SANGER", "syn3155834", ".somatic"),
    Source("UCSC-gatk_muse", "syn3107237", ".gatk_muse_0.9.9.5"),
    Source("UCSC-gatk_mutect", "syn3107237", ".gatk_mutect"),
    Source("UCSC-muse", "syn3107237", ".muse_0.9.9.5")]

def get_annotations(folder_id, prefix=".somatic"):
    synapse_ids   = []
    sample_ids    = []
    variant_types = []
    results = syn.chunkedQuery('select * from entity where parentId=="%s" order by name' % folder_id)
    for result in clean_query_results(results):
        if result['name'].endswith(prefix+".indel.vcf.gz"):
            variant_type = "indel"
        elif result['name'].endswith(prefix+".snv_mnv.vcf.gz"):
            variant_type = "snv"
        elif result['name'].endswith(prefix+".sv.vcf.gz"):
            variant_type = "sv"
        else:
            continue
        sample_id = extract_uuid(result['name'])
        synapse_id = result['id']

        if 'variant_type' in result and variant_type != result['variant_type'].lower():
            sys.stderr.write("Variant type {variant_type} doesn't match name {name}!".format(**result))

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
        dfs.append(df)

    return pd.concat(dfs)


df = create_metadata_df(sources)
print df.groupby(["source", "variant_type"])["synapse_id"].count()


import numpy as np
import matplotlib
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt


def plot_progress(df):
    ind = np.arange(len(sources)) + 0.1    # the x locations for the groups
    width = 0.80
    variant_types = ['snv', 'sv', 'indel']

    colors = {'snv':'#336699', 'sv':'#0000AA', 'indel':'#6633AA'}
    bars = {}
    bottoms = np.zeros(len(sources))
    for i, variant_type in enumerate(variant_types):
        counts = np.array(
            [len(df[ (df['source']==source.name) & (df['variant_type']==variant_type) ]) for source in sources])
        print i, variant_type, counts
        bars[variant_type] = plt.bar(ind, counts, width=width, bottom=bottoms, color=colors[variant_type])
        bottoms += counts

    plt.title('Count of VCF files by source and variant type')
    plt.xticks(ind+width/2., [source.name for source in sources] )
    #plt.yticks(np.arange(0,81,10))
    plt.legend([bars[variant_type] for variant_type in variant_types], variant_types)

    plt.show()


def plot_progress_old(data):

    sample_ids = list(set(data['sample_ids']))
    sample_ids.sort()
    n = len(sample_ids)

    fig, ax = plt.subplots()

    ys = np.linspace(0.95,0.05,n)
    patches = []
    for y in ys:
        patches.append(Rectangle((0.33,y), 0.1, 0.9/(1.1*n), fill=True))

    p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

    colors = 100*np.random.rand(len(patches))
    p.set_array(np.array(colors))

    ax.add_collection(p)

    for sample_id, y in izip(sample_ids, ys):
        ax.text(0.1,y,sample_id[0:8])

    plt.axis('equal')
    plt.axis('off')

    plt.show()


