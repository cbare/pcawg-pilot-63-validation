import synapseclient
from synapseclient import Project, File, Folder
from synapseclient import Schema, Column, Table, Row, RowSet, as_table_columns

syn.login('chris.bare')
project = syn.get('syn2875157')
examples_folder = syn.store(Folder("Examples", parent=project))

evaluation = syn.getEvaluation(3060780)


## submit a VCF file annotated with the tumor and normal sample IDs
entity = syn.store(File("foo.vcf", description="Test submission to PCAWG validation", parent=examples_folder,
        normal_sample_id="59b99c9a-d1c0-4cc0-953b-4358bcf9f392",
        tumor_sample_id="e9889071-4c6e-4761-9d65-06c6b5989fb7",
        variant_type="snv" ## must be: snv, sv, or indel
    ))

submission = syn.submit(evaluation, entity, name="Testing 1", teamName="SMC Admin team")


## submit using command line client
# synapse -u chris.bare store foo.vcf --parentId syn0000  --used ... --executed ... --annotations '{"normal_sample_id": "59b99c9a-d1c0-4cc0-953b-4358bcf9f392", "tumor_sample_id": "e9889071-4c6e-4761-9d65-06c6b5989fb7", "variant_type": "snv"}'
# synapse -u chris.bare submit --evalID 3060780 --entity syn0001


## get annotations from submitted File entity
import json
entity_bundle = json.loads(submission.entityBundleJSON)
properties = entity_bundle['entity']
annotations = synapseclient.annotations.from_synapse_annotations(entity_bundle['annotations'])
sub_entity = synapseclient.entity.Entity.create(properties, annotations)


