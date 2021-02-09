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

default_fdr_tr = 0.05

stats_dict = {0.01: {'possible_tf_asbs': 192020,
  'possible_cl_asbs': 328992,
  'possible_all_asbs': 240533,
  'possible_tf_asbs_rs': 139573,
  'possible_cl_asbs_rs': 224246,
  'possible_all_asbs_rs': 240503},
 0.05: {'possible_tf_asbs': 390916,
  'possible_cl_asbs': 569814,
  'possible_all_asbs': 377281,
  'possible_tf_asbs_rs': 249216,
  'possible_cl_asbs_rs': 349850,
  'possible_all_asbs_rs': 377221},
 0.1: {'possible_tf_asbs': 612656,
  'possible_cl_asbs': 796319,
  'possible_all_asbs': 501195,
  'possible_tf_asbs_rs': 354261,
  'possible_cl_asbs_rs': 460328,
  'possible_all_asbs_rs': 501114},
 0.25: {'possible_tf_asbs': 1693972,
  'possible_cl_asbs': 1671933,
  'possible_all_asbs': 990754,
  'possible_tf_asbs_rs': 803457,
  'possible_cl_asbs_rs': 889199,
  'possible_all_asbs_rs': 990574},
 0.5: {'possible_tf_asbs': 1693972,
  'possible_cl_asbs': 1671933,
  'possible_all_asbs': 990754,
  'possible_tf_asbs_rs': 803457,
  'possible_cl_asbs_rs': 889199,
  'possible_all_asbs_rs': 990574}
}



# possible_tf_asbs = 390916
# possible_cl_asbs = 569814
# possible_all_asbs = 960730
total_tf_candidates = 14630545
total_cl_candidates = 9791542
total_all_candidates = 24422087

# possible_tf_asbs_rs = 249216
# possible_cl_asbs_rs = 349850
# possible_all_asbs_rs = 377221
total_tf_candidates_rs = 3690724
total_cl_candidates_rs = 3690724
total_all_candidates_rs = 3690724
