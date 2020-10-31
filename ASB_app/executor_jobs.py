import os
import re
from datetime import datetime
from ASB_app import logger, executor
from ASB_app.constants import possible_tf_asbs, possible_cl_asbs, possible_cl_candidates, possible_all_asbs, \
    possible_all_candidates, possible_tf_candidates
from ASB_app.service import ananastra_service
from ASB_app.utils import pack, process_row
from sqlalchemy.orm import aliased
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from ASB_app.models import CandidateSNP

from ASB_app.releases import current_release

session = current_release.session
db = current_release.db

TranscriptionFactor, TranscriptionFactorSNP, CellLine, CellLineSNP, \
SNP, ExpSNP, Phenotype, PhenotypeSNPCorrespondence = \
    current_release.TranscriptionFactor, current_release.TranscriptionFactorSNP, current_release.CellLine, current_release.CellLineSNP, \
    current_release.SNP, current_release.ExpSNP, current_release.Phenotype, current_release.PhenotypeSNPCorrespondence


def convert_rs_to_int(rs_str):
    if not re.match(r'^rs\d+$', rs_str):
        raise ValueError(rs_str)
    return int(rs_str[2:])


def get_tf_query(rs_ids):
    grasp = aliased(Phenotype, name='grasp')
    ebi = aliased(Phenotype, name='ebi')
    clinvar = aliased(Phenotype, name='clinvar')
    finemapping = aliased(Phenotype, name='finemapping')
    qtl = aliased(Phenotype, name='qtl')
    phewas = aliased(Phenotype, name='phewas')

    return session.query(
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.chromosome)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.position)),
        db.func.group_concat(db.func.distinct(SNP.rs_id)),
        db.func.group_concat(db.func.distinct(SNP.ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactor.name)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.peak_calls)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.mean_bad)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.log_p_value_ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.log_p_value_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.es_ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.es_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_log_p_ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_log_p_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_log_2_fc)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_position)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_orientation)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_concordance)),
        db.func.group_concat(db.func.distinct(CellLine.name)),
        db.func.group_concat(db.func.distinct(qtl.phenotype_name)),
        db.func.group_concat(db.func.distinct(ebi.phenotype_name)),
        db.func.group_concat(db.func.distinct(phewas.phenotype_name)),
        db.func.group_concat(db.func.distinct(finemapping.phenotype_name)),
        db.func.group_concat(db.func.distinct(grasp.phenotype_name)),
        db.func.group_concat(db.func.distinct(clinvar.phenotype_name)),
    ).join(
        SNP,
        TranscriptionFactorSNP.snp
    ).filter(
        SNP.rs_id.in_(rs_ids)
    ).join(
        TranscriptionFactor,
        TranscriptionFactorSNP.transcription_factor
    ).join(
        PhenotypeSNPCorrespondence,
        (SNP.chromosome == PhenotypeSNPCorrespondence.chromosome) &
        (SNP.position == PhenotypeSNPCorrespondence.position) &
        (SNP.alt == PhenotypeSNPCorrespondence.alt)
    ).join(
        ExpSNP,
        TranscriptionFactorSNP.exp_snps
    ).join(
        CellLine,
        ExpSNP.cell_line
    ).join(
        qtl,
        (PhenotypeSNPCorrespondence.phenotype_id == qtl.phenotype_id) &
        (qtl.db_name == 'QTL'),
        isouter=True
    ).join(
        ebi,
        (PhenotypeSNPCorrespondence.phenotype_id == ebi.phenotype_id) &
        (ebi.db_name == 'ebi'),
        isouter=True
    ).join(
        phewas,
        (PhenotypeSNPCorrespondence.phenotype_id == phewas.phenotype_id) &
        (phewas.db_name == 'phewas'),
        isouter=True
    ).join(
        finemapping,
        (PhenotypeSNPCorrespondence.phenotype_id == finemapping.phenotype_id) &
        (finemapping.db_name == 'finemapping'),
        isouter=True
    ).join(
        grasp,
        (PhenotypeSNPCorrespondence.phenotype_id == grasp.phenotype_id) &
        (grasp.db_name == 'grasp'),
        isouter=True
    ).join(
        clinvar,
        (PhenotypeSNPCorrespondence.phenotype_id == clinvar.phenotype_id) &
        (clinvar.db_name == 'clinvar'),
        isouter=True
    ).group_by(TranscriptionFactorSNP.tf_snp_id)


def get_cl_query(rs_ids):
    grasp = aliased(Phenotype, name='grasp')
    ebi = aliased(Phenotype, name='ebi')
    clinvar = aliased(Phenotype, name='clinvar')
    finemapping = aliased(Phenotype, name='finemapping')
    qtl = aliased(Phenotype, name='qtl')
    phewas = aliased(Phenotype, name='phewas')

    return session.query(
        db.func.group_concat(db.func.distinct(CellLineSNP.chromosome)),
        db.func.group_concat(db.func.distinct(CellLineSNP.position)),
        db.func.group_concat(db.func.distinct(SNP.rs_id)),
        db.func.group_concat(db.func.distinct(SNP.ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.alt)),
        db.func.group_concat(db.func.distinct(CellLine.name)),
        db.func.group_concat(db.func.distinct(CellLineSNP.peak_calls)),
        db.func.group_concat(db.func.distinct(CellLineSNP.mean_bad)),
        db.func.group_concat(db.func.distinct(CellLineSNP.log_p_value_ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.log_p_value_alt)),
        db.func.group_concat(db.func.distinct(CellLineSNP.es_ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.es_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactor.name)),
        db.func.group_concat(db.func.distinct(qtl.phenotype_name)),
        db.func.group_concat(db.func.distinct(ebi.phenotype_name)),
        db.func.group_concat(db.func.distinct(phewas.phenotype_name)),
        db.func.group_concat(db.func.distinct(finemapping.phenotype_name)),
        db.func.group_concat(db.func.distinct(grasp.phenotype_name)),
        db.func.group_concat(db.func.distinct(clinvar.phenotype_name)),
    ).join(
        SNP,
        CellLineSNP.snp
    ).filter(
        SNP.rs_id.in_(rs_ids)
    ).join(
        CellLine,
        CellLineSNP.cell_line
    ).join(
        PhenotypeSNPCorrespondence,
        (SNP.chromosome == PhenotypeSNPCorrespondence.chromosome) &
        (SNP.position == PhenotypeSNPCorrespondence.position) &
        (SNP.alt == PhenotypeSNPCorrespondence.alt)
    ).join(
        ExpSNP,
        CellLineSNP.exp_snps
    ).join(
        TranscriptionFactor,
        ExpSNP.transcription_factor
    ).join(
        qtl,
        (PhenotypeSNPCorrespondence.phenotype_id == qtl.phenotype_id) &
        (qtl.db_name == 'QTL'),
        isouter=True
    ).join(
        ebi,
        (PhenotypeSNPCorrespondence.phenotype_id == ebi.phenotype_id) &
        (ebi.db_name == 'ebi'),
        isouter=True
    ).join(
        phewas,
        (PhenotypeSNPCorrespondence.phenotype_id == phewas.phenotype_id) &
        (phewas.db_name == 'phewas'),
        isouter=True
    ).join(
        finemapping,
        (PhenotypeSNPCorrespondence.phenotype_id == finemapping.phenotype_id) &
        (finemapping.db_name == 'finemapping'),
        isouter=True
    ).join(
        grasp,
        (PhenotypeSNPCorrespondence.phenotype_id == grasp.phenotype_id) &
        (grasp.db_name == 'grasp'),
        isouter=True
    ).join(
        clinvar,
        (PhenotypeSNPCorrespondence.phenotype_id == clinvar.phenotype_id) &
        (clinvar.db_name == 'clinvar'),
        isouter=True
    ).group_by(CellLineSNP.cl_snp_id)


def get_tf_asbs(rs_ids):
    return SNP.query.filter(
        SNP.rs_id.in_(rs_ids),
        SNP.tf_aggregated_snps.any()
    ).group_by(SNP.rs_id).count()


def get_cl_asbs(rs_ids):
    return SNP.query.filter(
        SNP.rs_id.in_(rs_ids),
        SNP.cl_aggregated_snps.any()
    ).group_by(SNP.rs_id).count()


def get_all_asbs(rs_ids):
    return SNP.query.filter(
        SNP.rs_id.in_(rs_ids)
    ).group_by(SNP.rs_id).count()


def get_tf_candidates(rs_ids):
    return CandidateSNP.query.filter(
        CandidateSNP.ag_level == 'TF',
        CandidateSNP.rs_id.in_(rs_ids)
    ).group_by(CandidateSNP.rs_id).count()


def get_cl_candidates(rs_ids):
    return CandidateSNP.query.filter(
        CandidateSNP.ag_level == 'CL',
        CandidateSNP.rs_id.in_(rs_ids)
    ).group_by(CandidateSNP.rs_id).count()


def get_all_candidates(rs_ids):
    return CandidateSNP.query.filter(
        CandidateSNP.rs_id.in_(rs_ids)
    ).group_by(CandidateSNP.rs_id).count()


def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def divide_query(get_query, rs_ids):
    for chunk in divide_chunks(rs_ids, 990):
        yield get_query(chunk)


@executor.job
def process_snp_file(ticket_id, annotate_tf=True, annotate_cl=True):
    input_file_name = ananastra_service.get_path_by_ticket_id(ticket_id)
    ticket = ananastra_service.get_ticket(ticket_id)
    ticket.status = 'Processing'
    try:
        data = pd.read_table(input_file_name, sep='\t', header=None)
        assert len(data.columns) == 1
        rs_ids = data[0].apply(convert_rs_to_int).unique().tolist()

        common_header_1 = ['CHROMOSOME', 'POSITION', 'RS_ID', 'REF', 'ALT']
        common_header_2 = ['PEAK_CALLS', 'MEAN_BAD', 'LOG10_FDR_REF', 'LOG10_FDR_ALT',
                           'EFFECT_SIZE_REF', 'EFFECT_SIZE_ALT']
        common_header_3 = ['GTEX EQTL', 'EBI', 'PHEWAS', 'FINEMAPPING', 'GRASP', 'CLINVAR']
        cl_header = common_header_1 + ['CELL_TYPE'] + common_header_2 + ['AGGREGATED_TFS'] + common_header_3
        tf_header = common_header_1 + ['TRANSCRIPTION_FACTOR'] + common_header_2 + \
                    ['MOTIF_LOG_P_REF', 'MOTIF_LOG_P_ALT', 'MOTIF_LOG2_FC', 'MOTIF_POSITION',
                     'MOTIF_ORIENTATION', 'MOTIF_CONCORDANCE', 'AGGREGATED_CELL_TYPES'] + common_header_3

        if annotate_tf:
            ananastra_service.create_processed_path(ticket_id, 'tf')
            tf_path = ananastra_service.get_path_by_ticket_id(ticket_id, path_type='tf', ext='.tsv')

            with open(tf_path, 'w') as out:
                out.write(pack(tf_header))

            for q_tf in divide_query(get_tf_query, rs_ids):
                with open(tf_path, 'a') as out:
                    for tup in q_tf:
                        out.write(pack(process_row(tup, 'TF', tf_header)))

            ananastra_service.create_processed_path(ticket_id, 'tf_sum')

            tf_table = pd.read_table(tf_path)
            tf_table['BEST_FDR'] = tf_table[['LOG10_FDR_REF', 'LOG10_FDR_ALT']].max(axis=1)
            idx = tf_table.groupby(['RS_ID', 'ALT'])['BEST_FDR'].transform(max) == tf_table['BEST_FDR']
            tf_table.drop(columns=['BEST_FDR'], inplace=True)
            tf_table[idx].to_csv(ananastra_service.get_path_by_ticket_id(ticket_id, 'tf_sum'), sep='\t', index=False)

        if annotate_cl:
            ananastra_service.create_processed_path(ticket_id, 'cl')
            cl_path = ananastra_service.get_path_by_ticket_id(ticket_id, path_type='cl', ext='.tsv')

            with open(cl_path, 'w') as out:
                out.write(pack(cl_header))

            for q_cl in divide_query(get_cl_query, rs_ids):
                with open(cl_path, 'a') as out:
                    for tup in q_cl:
                        out.write(pack(process_row(tup, 'CL', cl_header)))

            ananastra_service.create_processed_path(ticket_id, 'cl_sum')

            cl_table = pd.read_table(cl_path)
            cl_table['BEST_FDR'] = cl_table[['LOG10_FDR_REF', 'LOG10_FDR_ALT']].max(axis=1)
            idx = cl_table.groupby(['RS_ID', 'ALT'])['BEST_FDR'].transform(max) == cl_table['BEST_FDR']
            cl_table.drop(columns=['BEST_FDR'], inplace=True)
            cl_table[idx].to_csv(ananastra_service.get_path_by_ticket_id(ticket_id, 'cl_sum'), sep='\t', index=False)

        all_rs = len(rs_ids)
        tf_asbs = sum(divide_query(get_tf_asbs, rs_ids))
        cl_asbs = sum(divide_query(get_cl_asbs, rs_ids))
        all_asbs = sum(divide_query(get_all_asbs, rs_ids))
        tf_candidates = sum(divide_query(get_tf_candidates, rs_ids))
        cl_candidates = sum(divide_query(get_cl_candidates, rs_ids))
        all_candidates = sum(divide_query(get_all_candidates, rs_ids))

        if tf_candidates:
            tf_odds, tf_p = fisher_exact(((tf_asbs, tf_candidates), (possible_tf_asbs, possible_tf_candidates)))
        else:
            tf_odds, tf_p = 0, 1

        if cl_candidates:
            cl_odds, cl_p = fisher_exact(((cl_asbs, cl_candidates), (possible_cl_asbs, possible_cl_candidates)))
        else:
            cl_odds, cl_p = 0, 1

        if all_candidates:
            all_odds, all_p = fisher_exact(((all_asbs, all_candidates), (possible_all_asbs, possible_all_candidates)))
        else:
            all_odds, all_p = 0, 1

    except Exception as e:
        logger.error(e, exc_info=True)
        ticket.status = 'Failed'
        session.commit()
        return

    ticket.status = 'Processed'
    ticket.meta_info = {
        'processing_time': str(datetime.now() - ticket.date_created),
        'all_rs': all_rs,
        'tf_asbs': tf_asbs,
        'cl_asbs': cl_asbs,
        'all_asbs': all_asbs,
        'tf_candidates': tf_candidates,
        'cl_canidates': cl_candidates,
        'all_candidates': all_candidates,
        'tf_odds': tf_odds,
        'tf_log10_p_value': -np.log10(tf_p),
        'cl_odds': cl_odds,
        'cl_log10_p_value': -np.log10(cl_p),
        'all_odds': all_odds,
        'all_log10_p_value': -np.log10(all_p),
    }
    session.commit()
