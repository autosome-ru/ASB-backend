from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from ASB_app import *
from ASB_app.releases import current_release
from ASB_app.models import CandidateSNP

import pandas as pd
import os

from ASB_app.constants import chromosomes, possible_all_asbs_rs, possible_all_candidates_rs, possible_tf_candidates_rs, possible_cl_candidates_rs, possible_cl_asbs_rs, possible_tf_asbs_rs

SNP = current_release.SNP
Gene = current_release.Gene
TFSNP = current_release.TranscriptionFactorSNP
CLSNP = current_release.CellLineSNP

# possible_all_asbs_rs_chr = dict(zip(chromosomes, [SNP.query.filter(SNP.chromosome == chr).group_by(SNP.rs_id).count() for chr in chromosomes]))
# possible_tf_asbs_rs_chr = dict(zip(chromosomes, [SNP.query.filter(SNP.tf_aggregated_snps.any(), SNP.chromosome == chr).group_by(SNP.rs_id).count() for chr in chromosomes]))
# possible_cl_asbs_rs_chr = dict(zip(chromosomes, [SNP.query.filter(SNP.cl_aggregated_snps.any(), SNP.chromosome == chr).group_by(SNP.rs_id).count() for chr in chromosomes]))
#
# possible_all_candidates_rs_chr = dict(zip(chromosomes, [CandidateSNP.query.filter(CandidateSNP.chromosome == chr).group_by(CandidateSNP.rs_id).count() for chr in chromosomes]))
# possible_tf_candidates_rs_chr = dict(zip(chromosomes, [CandidateSNP.query.filter(CandidateSNP.ag_level == 'TF', CandidateSNP.chromosome == chr).group_by(CandidateSNP.rs_id).count() for chr in chromosomes]))
# possible_cl_candidates_rs_chr = dict(zip(chromosomes, [CandidateSNP.query.filter(CandidateSNP.ag_level == 'CL', CandidateSNP.chromosome == chr).group_by(CandidateSNP.rs_id).count() for chr in chromosomes]))

possible_all_asbs_rs_chr = {'chr1': 37101,
 'chr2': 27983,
 'chr3': 20593,
 'chr4': 16419,
 'chr5': 19947,
 'chr6': 25568,
 'chr7': 24044,
 'chr8': 19555,
 'chr9': 13062,
 'chr10': 16948,
 'chr11': 18070,
 'chr12': 19327,
 'chr13': 7955,
 'chr14': 12425,
 'chr15': 12118,
 'chr16': 14419,
 'chr17': 18927,
 'chr18': 7480,
 'chr19': 16558,
 'chr20': 12783,
 'chr21': 5529,
 'chr22': 7719,
 'chrX': 2691,
 'chrY': 0}
possible_tf_asbs_rs_chr = {'chr1': 23431,
 'chr2': 19272,
 'chr3': 14977,
 'chr4': 11487,
 'chr5': 13970,
 'chr6': 16777,
 'chr7': 13971,
 'chr8': 12844,
 'chr9': 9444,
 'chr10': 11860,
 'chr11': 12354,
 'chr12': 12566,
 'chr13': 5744,
 'chr14': 8051,
 'chr15': 7896,
 'chr16': 9214,
 'chr17': 11427,
 'chr18': 5043,
 'chr19': 10347,
 'chr20': 8411,
 'chr21': 3679,
 'chr22': 4642,
 'chrX': 1809,
 'chrY': 0}
possible_cl_asbs_rs_chr = {'chr1': 34493,
 'chr2': 25721,
 'chr3': 18818,
 'chr4': 15235,
 'chr5': 18462,
 'chr6': 23868,
 'chr7': 22526,
 'chr8': 18075,
 'chr9': 11897,
 'chr10': 15634,
 'chr11': 16706,
 'chr12': 17999,
 'chr13': 7329,
 'chr14': 11429,
 'chr15': 11231,
 'chr16': 13296,
 'chr17': 17744,
 'chr18': 6943,
 'chr19': 15714,
 'chr20': 11835,
 'chr21': 5153,
 'chr22': 7252,
 'chrX': 2490,
 'chrY': 0}

possible_all_candidates_rs_chr = {'chr1': 318694,
 'chr2': 308332,
 'chr3': 251887,
 'chr4': 203621,
 'chr5': 231700,
 'chr6': 240129,
 'chr7': 224487,
 'chr8': 195748,
 'chr9': 147078,
 'chr10': 178428,
 'chr11': 175332,
 'chr12': 182963,
 'chr13': 101788,
 'chr14': 125922,
 'chr15': 106522,
 'chr16': 120852,
 'chr17': 131674,
 'chr18': 84116,
 'chr19': 104707,
 'chr20': 113575,
 'chr21': 47801,
 'chr22': 56046,
 'chrX': 39312,
 'chrY': 10}
possible_tf_candidates_rs_chr = {'chr1': 318694,
 'chr2': 308332,
 'chr3': 251887,
 'chr4': 203621,
 'chr5': 231700,
 'chr6': 240129,
 'chr7': 224487,
 'chr8': 195748,
 'chr9': 147078,
 'chr10': 178428,
 'chr11': 175332,
 'chr12': 182963,
 'chr13': 101788,
 'chr14': 125922,
 'chr15': 106522,
 'chr16': 120852,
 'chr17': 131674,
 'chr18': 84116,
 'chr19': 104707,
 'chr20': 113575,
 'chr21': 47801,
 'chr22': 56046,
 'chrX': 39312,
 'chrY': 10}
possible_cl_candidates_rs_chr = {'chr1': 318694,
 'chr2': 308332,
 'chr3': 251887,
 'chr4': 203621,
 'chr5': 231700,
 'chr6': 240129,
 'chr7': 224487,
 'chr8': 195748,
 'chr9': 147078,
 'chr10': 178428,
 'chr11': 175332,
 'chr12': 182963,
 'chr13': 101788,
 'chr14': 125922,
 'chr15': 106522,
 'chr16': 120852,
 'chr17': 131674,
 'chr18': 84116,
 'chr19': 104707,
 'chr20': 113575,
 'chr21': 47801,
 'chr22': 56046,
 'chrX': 39312,
 'chrY': 10}

gene_data = []
for i, gene in enumerate(Gene.query, 1):
    if i % 1000 == 0:
        print(i)
    snps = SNP.query.filter(SNP.chromosome == gene.chromosome,
                            SNP.position.between(gene.start_pos, gene.end_pos)).group_by(SNP.rs_id).count()
    tf_snps = SNP.query.filter(SNP.tf_aggregated_snps.any(), SNP.chromosome == gene.chromosome,
                            SNP.position.between(gene.start_pos, gene.end_pos)).group_by(SNP.rs_id).count()
    cl_snps = SNP.query.filter(SNP.cl_aggregated_snps.any(), SNP.chromosome == gene.chromosome,
                            SNP.position.between(gene.start_pos, gene.end_pos)).group_by(SNP.rs_id).count()
    cands = CandidateSNP.query.filter(CandidateSNP.chromosome == gene.chromosome,
                                      CandidateSNP.position.between(gene.start_pos, gene.end_pos)).group_by(CandidateSNP.rs_id).count()
    tf_cands = CandidateSNP.query.filter(CandidateSNP.ag_level == 'TF', CandidateSNP.chromosome == gene.chromosome,
                                      CandidateSNP.position.between(gene.start_pos, gene.end_pos)).group_by(CandidateSNP.rs_id).count()
    cl_cands = CandidateSNP.query.filter(CandidateSNP.ag_level == 'CL', CandidateSNP.chromosome == gene.chromosome,
                                      CandidateSNP.position.between(gene.start_pos, gene.end_pos)).group_by(CandidateSNP.rs_id).count()
    all_odds, all_p = fisher_exact(((snps, cands), (possible_all_asbs_rs, possible_all_candidates_rs)), alternative='greater')
    tf_odds, tf_p = fisher_exact(((tf_snps, tf_cands), (possible_tf_asbs_rs, possible_tf_candidates_rs)), alternative='greater')
    cl_odds, cl_p = fisher_exact(((cl_snps, cl_cands), (possible_cl_asbs_rs, possible_cl_candidates_rs)), alternative='greater')
    all_odds_chr, all_p_chr = fisher_exact(((snps, cands), (possible_all_asbs_rs_chr[gene.chromosome], possible_all_candidates_rs_chr[gene.chromosome])), alternative='greater')
    tf_odds_chr, tf_p_chr = fisher_exact(((tf_snps, tf_cands), (possible_tf_asbs_rs_chr[gene.chromosome], possible_tf_candidates_rs_chr[gene.chromosome])), alternative='greater')
    cl_odds_chr, cl_p_chr = fisher_exact(((cl_snps, cl_cands), (possible_cl_asbs_rs_chr[gene.chromosome], possible_cl_candidates_rs_chr[gene.chromosome])), alternative='greater')
    gene_data.append((gene.gene_name, gene.chromosome, gene.start_pos, gene.end_pos, '+' if gene.orientation else '-', int(snps), int(cands), all_p, all_odds, all_p_chr, all_odds_chr, int(cl_snps), int(cl_cands), cl_p, cl_odds, cl_p_chr, cl_odds_chr, int(tf_snps), int(tf_cands), tf_p, tf_odds, tf_p_chr, tf_odds_chr))
t = pd.DataFrame(gene_data)
t.to_csv(os.path.expanduser('~/Desktop/all_gene.csv'), header=False, index=False)
