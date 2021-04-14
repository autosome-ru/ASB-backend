from ASB_app import models, releases

chromosomes = ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
               'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19',
               'chr20', 'chr21', 'chr22', 'chrX', 'chrY')

nucleotides = ('A', 'C', 'G', 'T')

bads = ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')

db_name_property_dict = {
    x: 'has_{}_associations'.format(x.lower())
    for x in ['clinvar', 'ebi', 'finemapping', 'grasp', 'phewas', 'QTL']
}

default_fdr_tr = lambda v: '0.1' if v >= 3 else '0.05'
default_es_tr = lambda v: 'all' if v >= 3 else 'all'
fdr_choices = ['0.01', '0.05', '0.1', '0.15', '0.25']
es_choices = ['1.2', '0.6', '0']
fdr_classes = fdr_choices + ['1']
es_classes = es_choices + ['all']

stats_dict = {'0.01': {'possible_tf_asbs': 100178,
                       'possible_cl_asbs': 156005,
                       'possible_all_asbs': 256183,
                       'possible_tf_asbs_rs': 79713,
                       'possible_cl_asbs_rs': 115757,
                       'possible_all_asbs_rs': 126978},
              '0.05': {'possible_tf_asbs': 183756,
                       'possible_cl_asbs': 262013,
                       'possible_all_asbs': 445769,
                       'possible_tf_asbs_rs': 136045,
                       'possible_cl_asbs_rs': 180299,
                       'possible_all_asbs_rs': 198235},
              '0.1': {'possible_tf_asbs': 255812,
                      'possible_cl_asbs': 349490,
                      'possible_all_asbs': 605302,
                      'possible_tf_asbs_rs': 180149,
                      'possible_cl_asbs_rs': 229260,
                      'possible_all_asbs_rs': 253182},
              '0.15': {'possible_tf_asbs': 323556,
                       'possible_cl_asbs': 427750,
                       'possible_all_asbs': 751306,
                       'possible_tf_asbs_rs': 219099,
                       'possible_cl_asbs_rs': 270903,
                       'possible_all_asbs_rs': 300951},
              '0.25': {'possible_tf_asbs': 468555,
                       'possible_cl_asbs': 585750,
                       'possible_all_asbs': 1054305,
                       'possible_tf_asbs_rs': 297146,
                       'possible_cl_asbs_rs': 349733,
                       'possible_all_asbs_rs': 393246},
              '0.5': {'possible_tf_asbs': 468555,  # FIXME
                      'possible_cl_asbs': 585750,
                      'possible_all_asbs': 1054305,
                      'possible_tf_asbs_rs': 297146,
                      'possible_cl_asbs_rs': 349733,
                      'possible_all_asbs_rs': 253630}
              }

# possible_tf_asbs = 390916
# possible_cl_asbs = 569814
# possible_all_asbs = 960730
total_tf_candidates = 14630545  # FIXME
total_cl_candidates = 9791542
total_all_candidates = 24422087

# possible_tf_asbs_rs = 249216
# possible_cl_asbs_rs = 349850
# possible_all_asbs_rs = 377221
total_tf_candidates_rs = 3690724
total_cl_candidates_rs = 3690724
total_all_candidates_rs = 3690724
