from ASB_app import *

from ASB_app.releases import ReleaseFord
import numpy as np

release = ReleaseFord
Gene = release.Gene
SNP = release.SNP
TFSNP = release.TranscriptionFactorSNP

session = release.session


def get_filters_by_gene(self, gene):
    if gene.orientation:
        return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos - 5000, 1),
                                                                                 gene.end_pos)
    else:
        return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos, 1),
                                                                                 gene.end_pos + 5000)

q = session.query(Gene, SNP, TFSNP)\
    .join(SNP, Gene.snps_by_target)\
    .join(TFSNP, SNP.tf_aggregated_snps)\
    .filter(TFSNP.best_p_value > np.log10(20))\
    .limit(10)

for (gene, snp, tfsnp) in q:
    print(gene.gene_name, snp.rs_id, tfsnp.p_value_ref)
