import os
import sys

from ASB_app import *

from ASB_app.releases import ReleaseBillCipher
import numpy as np

release = ReleaseBillCipher
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
        ).filter(
            (ExpSNP.p_value_ref - ExpSNP.p_value_alt) * (
                        AGSNP.log_p_value_alt - AGSNP.log_p_value_ref) > 0
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
            Gene.proximal_promoter_snps
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
        ).filter(
            (ExpSNP.p_value_ref - ExpSNP.p_value_alt) * (
                        AGSNP.log_p_value_alt - AGSNP.log_p_value_ref) > 0
        ).join(
            Experiment,
            ExpSNP.experiment,
        ).join(
            OTHER,
            getattr(Experiment, 'transcription_factor' if AG == CL else 'cell_line')
        )

    promoter_dict = {}
    target_dict = {}
    for i, (q, q_dict) in enumerate(((q_promoter, promoter_dict), (q_target, target_dict))):
        print(i)
        for (gene, snp, agsnp, ag, other_name) in q:
            if (gene.gene_id, snp.rs_id, snp.alt) in q_dict:
                q_dict[(gene.gene_id, snp.rs_id, snp.alt)][-1].add(other_name)
            else:
                q_dict[(gene.gene_id, snp.rs_id, snp.alt)] = [gene.gene_name, gene.gene_id, snp.chromosome, snp.position,
                                        'rs' + str(snp.rs_id), snp.ref, snp.alt, snp.position - gene.start_pos, ag.name,
                                        '{} ({})'.format(*(('ref', snp.ref) if agsnp.log_p_value_ref > agsnp.log_p_value_alt else ('alt', snp.alt))),
                                        max(agsnp.log_p_value_ref, agsnp.log_p_value_alt),
                                        agsnp.es_ref if agsnp.log_p_value_ref > agsnp.log_p_value_alt else agsnp.es_alt,
                                        {other_name}]

    all_keys = list(set(promoter_dict.keys()) | set(target_dict.keys()))

    with open(os.path.expanduser('~/{}_genes_all_{}.tsv'.format('tf' if AG == TF else 'cl', release.name)), 'w') as out:
        out.write(
            '\t'.join(
                map(str, [
                    'Gene_name',
                    'Gene_id'
                    'Chromosome',
                    'Position',
                    'rs_ID',
                    'Ref',
                    'Alt',
                    'Distance_to_TSS',
                    'TF' if AG == TF else 'Cell_type',
                    'Preferred_allele',
                    'Log10_FDR',
                    'Effect_size(log2)',
                    'Supporting_{}'.format('TFs' if AG == CL else 'Cell_types'),
                    'eQTL',
                    'promoter_SNP',
                ])
            ) + '\n'
        )
        for key in all_keys:
            if key in target_dict:
                data = target_dict[key]
            else:
                data = promoter_dict[key]
            out.write(
                '\t'.join(
                    map(str, data[:-1] + ['|'.join(data[-1])] + [key in target_dict, key in promoter_dict])
                ) + '\n'
            )
