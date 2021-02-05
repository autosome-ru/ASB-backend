import os
import sys

from ASB_app import *

from ASB_app.releases import ReleaseFord
import numpy as np

release = ReleaseFord
Gene = release.Gene
SNP = release.SNP
TFSNP = release.TranscriptionFactorSNP
CLSNP = release.CellLineSNP
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
    if sys.argv[1] == 'TF':
        AGSNP = TFSNP
        AG = TF
        OTHER = CL
    elif sys.argv[1] == 'CL':
        AGSNP = CLSNP
        AG = CL
        OTHER = TF
    else:
        raise ValueError

    q_target = session.query(
            Gene,
            SNP,
            AGSNP,
            AG,
            OTHER.name
        ).join(
            SNP,
            Gene.snps_by_target
        ).join(
            AGSNP,
            getattr(SNP, 'tf_aggregated_snps' if AG == TF else 'cl_aggregated_snps')
        ).filter(
            AGSNP.best_p_value > np.log10(20)
        ).join(
            AG,
            getattr(AGSNP, 'transcription_factor' if AG == TF else 'cell_line')
        ).join(
            ExpSNP,
            AGSNP.exp_snps
        # ).filter(
        #     (ExpSNP.p_value_ref - ExpSNP.p_value_alt) * (
        #                 AGSNP.log_p_value_alt - AGSNP.log_p_value_ref) > 0
        ).join(
            Experiment,
            ExpSNP.experiment,
        ).join(
            OTHER,
            getattr(Experiment, 'transcription_factor' if AG == CL else 'cell_line')
        )

    q_promoter = session.query(
            Gene,
            SNP,
            AGSNP,
            AG,
            OTHER.name
        ).join(
            SNP,
            (SNP.chromosome == Gene.chromosome) &
            (
                (Gene.orientation & SNP.position.between(Gene.start_pos - 5000,Gene.end_pos)) |
                (~Gene.orientation & SNP.position.between(Gene.start_pos, Gene.end_pos + 5000))
            )
        ).join(
            AGSNP,
            getattr(SNP, 'tf_aggregated_snps' if AG == TF else 'cl_aggregated_snps')
        ).filter(
            AGSNP.best_p_value > np.log10(20)
        ).join(
            AG,
            getattr(AGSNP, 'transcription_factor' if AG == TF else 'cell_line')
        ).join(
            ExpSNP,
            AGSNP.exp_snps
        # ).filter(
        #     (ExpSNP.p_value_ref - ExpSNP.p_value_alt) * (
        #                 AGSNP.log_p_value_alt - AGSNP.log_p_value_ref) > 0
        ).join(
            Experiment,
            ExpSNP.experiment,
        ).join(
            OTHER,
            getattr(Experiment, 'transcription_factor' if AG == CL else 'cell_line')
        )

    promoter_dict = {}
    target_dict = {}
    for q, q_dict in (q_promoter, promoter_dict), (q_target, target_dict):
        for (gene, snp, agsnp, ag, other_name) in q:
            if gene.gene_id in q_dict:
                q_dict[gene.gene_id][-1].append(other_name)
            else:
                q_dict[gene.gene_id] = [gene.gene_name, snp.chromosome, snp.position,
                                        'rs' + str(snp.rs_id), snp.ref, snp.alt, ag.name,
                                        '{} ({})'.format(*(('ref', snp.ref) if agsnp.log_p_value_ref > agsnp.log_p_value_alt else ('alt', snp.alt))),
                                        max(agsnp.log_p_value_ref, agsnp.log_p_value_alt),
                                        agsnp.es_ref if agsnp.log_p_value_ref > agsnp.log_p_value_alt else agsnp.es_alt,
                                        [other_name]]

    all_genes = list(set(promoter_dict.keys()) | set(target_dict.keys()))

    with open(os.path.expanduser('~/{}_genes_all.tsv'.format('tf' if AG == TF else 'cl')), 'w') as out:
        out.write(
            '\t'.join(
                map(str, [
                    'Gene_name',
                    'Chromosome',
                    'Position'
                    'rs_ID',
                    'Ref',
                    'Alt',
                    'TF' if AG == TF else 'Cells',
                    'Preferred_allele',
                    'Log10_p_value',
                    'Effect_size(log2)',
                    '{}'.format('TFs' if AG == CL else 'Cell_types'), ###
                    'eQTL',
                    'promoter_SNP',
                ])
            )
        )
        for gene_id in all_genes:
            if gene_id in target_dict:
                data = target_dict[gene_id]
            else:
                data = promoter_dict[gene_id]
            out.write(
                '\t'.join(
                    map(str, data[:-1] + ['|'.join(data[-1])] + [gene_id in target_dict, gene_id in promoter_dict])
                )
            )
