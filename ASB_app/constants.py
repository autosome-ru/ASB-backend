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

possible_tf_asbs = 147909
possible_cl_asbs = 252469
possible_all_asbs = 269934
possible_tf_candidates = 4976303
possible_cl_candidates = 4976303
possible_all_candidates = 4976303
