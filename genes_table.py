from ASB_app import *

from ASB_app.releases import ReleaseFord
import numpy as np

release = ReleaseFord
Gene = release.Gene
SNP = release.SNP
TFSNP = release.TranscriptionFactorSNP
TF = release.TranscriptionFactor
CL = release.CellLine
ExpSNP = release.ExpSNP
Experiment = release.Experiment

session = release.session
db = release.db


def get_filters_by_gene(self, gene):
    if gene.orientation:
        return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos - 5000, 1),
                                                                                 gene.end_pos)
    else:
        return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos, 1),
                                                                                 gene.end_pos + 5000)


if __name__ == '__main__':
    q = session.query(
            Gene,
            SNP,
            TFSNP,
            TF,
            db.func.group_concat(db.func.distinct(CL.name), separator=', ')
        ).join(
            SNP,
            Gene.snps_by_target
        ).join(
            TFSNP,
            SNP.tf_aggregated_snps
        ).filter(
            TFSNP.best_p_value > np.log10(20)
        ).join(
            TF,
            TFSNP.transcription_factor
        ).join(
            ExpSNP,
            TFSNP.exp_snps
        ).filter(
            (ExpSNP.p_value_ref - ExpSNP.p_value_alt) * (
                        TFSNP.log_p_value_alt - TFSNP.log_p_value_ref) > 0
        ).join(
            Experiment,
            ExpSNP.experiment,
        ).join(
            CL,
            Experiment.cell_line,
        ).limit(
            10
        ).from_self(
        ).group_by(
            Gene.gene_id,
            SNP.chromosome,
            SNP.position,
            SNP.alt,
            TFSNP.tf_snp_id,
        )


for (gene, snp, tfsnp, tf, cl_names) in q:
    print(gene.gene_name, 'rs' + str(snp.rs_id), tfsnp.log_p_value_ref, tfsnp.log_p_value_ref, tf.name, cl_names)
