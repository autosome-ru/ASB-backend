import sys
import pandas as pd

from ASB_app.releases import ReleaseZanthar as r

if __name__ == '__main__':
    uniprot_conv = pd.read_table('~/PARAMETERS/gene_names.tab')
    tfs = []
    for index, row in uniprot_conv.iterrows():
        tf = r.TranscriptionFactor.query.filter(
            r.TranscriptionFactor.name == row['Entry name']
        ).one()
        print(tf.name)
        if not pd.isna(row['Gene names']):
            if row['Gene names'] == 'T-Cell Receptor V-alpha region':
                tf.gene_name = row['Gene names']
            else:
                tf.gene_name = row['Gene names'].split(' ')[0]
    r.session.commit()
