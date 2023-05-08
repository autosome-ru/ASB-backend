import json
import os

import pathlib

ananastra_stats_file = os.path.join(pathlib.Path(__file__).parent.absolute(), 'stats.json')

chromosomes = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
               'chr20', 'chr21', 'chr22', 'chrX', 'chrY')

nucleotides = ('A', 'C', 'G', 'T')

bads = ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')
max_nrows = 10000
max_comments = 5000
db_name_property_dict = {
    x: 'has_{}_associations'.format(x.lower())
    for x in ['clinvar', 'ebi', 'finemapping', 'grasp', 'phewas', 'QTL']
}

default_fdr_tr = lambda v: '0.1' if v == 3 else '0.05'
default_es_tr = lambda v: 'all' if v >= 3 else 'all'
fdr_choices = ['0.01', '0.05', '0.1', '0.15', '0.25']
es_choices = ['2', '1', '0']
fdr_classes = fdr_choices + ['1']
es_classes = es_choices + ['all']
background_choices = ('WG', 'LOCAL', 'LD-ASN', 'LD-EUR', 'LD-AFR')


def read_ananastra_constants():
    if not os.path.isfile(ananastra_stats_file):
        print('No ananastra statistics file!')
        return {
            'stats_dict': {},
            'tf_stats_dict': {},
            'cl_stats_dict': {},
            'chr_stats_dict': {},
        }
    with open(ananastra_stats_file) as f:
        return json.load(f)


ananastra_constants = read_ananastra_constants()


stats_dict = ananastra_constants

