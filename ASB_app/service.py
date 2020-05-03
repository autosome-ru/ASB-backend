from ASB_app import db
from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, SNP, TranscriptionFactor, CellLine
from ASB_app.exceptions import ParsingError
from ASB_app.utils.aggregates import db_name_property_dict, TsvDialect
from sqlalchemy import not_
import csv
import tempfile
from flask import send_file


def get_filters_by_rs_id(rs_id):
    return (SNP.rs_id == rs_id, )


def get_full_snp(rs_id, alt):
    return SNP.query.filter(
        (SNP.rs_id == rs_id) &
        (SNP.alt == alt)
    ).one()


def get_full_snp_tsv(what_for, rs_id, alt, headers):
    agg_snps = getattr(get_full_snp(rs_id, alt), what_for + '_aggregated_snps')

    file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
    csv_writer = csv.writer(file, dialect=TsvDialect)

    csv_writer.writerow(headers)

    getter = lambda snp, header: getter(getattr(snp, header[:header.find('.')]),
                                        header[header.find('.') + 1:]) if '.' in header else getattr(snp, header)

    for snp in agg_snps:
        csv_writer.writerow([getter(snp, header) for header in headers])

    file.flush()

    return send_file(
        file.name,
        cache_timeout=0,
        mimetype="text/tsv",
        as_attachment=True
    )


def get_filters_by_genome_position(chr, pos1, pos2):
    return SNP.chromosome == chr, SNP.position.between(pos1, pos2)


def construct_advanced_filters(filters_object):
    filters = []
    if filters_object['transcription_factors']:
        filters += [SNP.tf_aggregated_snps.any(
            TranscriptionFactorSNP.tf_id == getattr(TranscriptionFactor.query.filter(
                TranscriptionFactor.name == tf_name
            ).one_or_none(), 'tf_id', None))
            for tf_name in filters_object['transcription_factors']]

    if filters_object['cell_types']:
        filters += [SNP.cl_aggregated_snps.any(
            CellLineSNP.cl_id == getattr(CellLine.query.filter(
                CellLine.name == cl_name
            ).one_or_none(), 'cl_id', None))
            for cl_name in filters_object['cell_types']]

    if filters_object['chromosome']:
        if not filters_object['start'] or not filters_object['end']:
            filters += [SNP.chromosome == filters_object['chromosome']]
        else:
            filters += [SNP.chromosome == filters_object['chromosome'],
                        SNP.position.between(filters_object['start'], filters_object['end'])]

    if filters_object['phenotype_databases']:
        filters += [getattr(SNP, db_name_property_dict[phenotype_db])
                    for phenotype_db in filters_object['phenotype_databases']]

    return filters


def get_snps_by_advanced_filters_tsv(filters_object):
    found_snps = SNP.query.filter(*construct_advanced_filters(filters_object)).all()

    file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
    csv_writer = csv.writer(file, dialect=TsvDialect)

    headers = ['chromosome', 'position', 'ref', 'alt']

    csv_writer.writerow(headers)

    for snp in found_snps:
        csv_writer.writerow([getattr(snp, header) for header in headers])

    file.flush()
    return send_file(
        file.name,
        cache_timeout=0,
        mimetype="text/tsv",
        as_attachment=True
    )

# import numpy as np
# def save_for_sarus():
#     found_snps = TranscriptionFactorSNP.query.filter(
#         TranscriptionFactorSNP.tf_id != 75,
#         db.func.max(TranscriptionFactorSNP.motif_log_p_alt, TranscriptionFactorSNP.motif_log_p_ref) >= -np.log10(0.0001),
#         db.func.max(TranscriptionFactorSNP.motif_log_2_fc, -1 * TranscriptionFactorSNP.motif_log_2_fc) >= 2,
#     )
#     import os
#     file = open(os.path.expanduser('~/Top10TFs/All_but_CTCF.tsv'), 'w')
#     csv_writer = csv.writer(file, dialect=TsvDialect)
#
#     headers = ['chromosome', 'position', 'snp.ref', 'alt', 'log_p_value_ref', 'log_p_value_alt', 'motif_log_2_fc']
#
#     csv_writer.writerow(['chr', 'pos', 'ref', 'alt', 'log_p_value_ref', 'log_p_value_alt', 'motif_log_2_fc'])
#
#     getter = lambda snp, header: getter(getattr(snp, header[:header.find('.')]),
#                                         header[header.find('.') + 1:]) if '.' in header else getattr(snp, header)
#
#     for snp in found_snps:
#         csv_writer.writerow([getter(snp, header) for header in headers])
#
#     file.flush()
#     file.close()


def get_hints(what_for, in_str, used_options):
    cls = {'TF': TranscriptionFactor, 'CL': CellLine}[what_for]
    filters = (((cls.name.like(in_str),) if in_str else ()) +
               ((not_(cls.name.in_(used_options)),) if used_options else ()) +
               (cls.aggregated_snps_count,))
    return cls.query.filter(*filters).order_by(cls.aggregated_snps_count.desc()).limit(3).all()


def get_overall_statistics():
    return {
        'transcription_factors_count': TranscriptionFactor.query.filter(TranscriptionFactor.aggregated_snps_count > 0).count(),
        'cell_types_count': CellLine.query.filter(CellLine.aggregated_snps_count > 0).count(),
        'snps_count': SNP.query.count(),  # FIXME: consider .group_by(SNP.rs_id)
    }
