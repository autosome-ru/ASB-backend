from sqlalchemy.orm import aliased

from ASB_app import db, session
from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, SNP, TranscriptionFactor, CellLine, Gene
from ASB_app.exceptions import ParsingError
from ASB_app.utils.aggregates import db_name_property_dict, TsvDialect
from sqlalchemy import not_
import csv
import tempfile
from flask import send_file


def get_filters_by_rs_id(rs_id):
    return (SNP.rs_id == rs_id,)


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


def get_gene_by_name(gene_name):
    return Gene.query.filter_by(gene_name=gene_name).one_or_none()


def get_gene_by_id(gene_id):
    return Gene.query.get(gene_id)


def get_filters_by_gene(gene):
    return SNP.chromosome == gene.chromosome, SNP.position.between(max(gene.start_pos - 1000, 1), gene.end_pos)


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

    if filters_object['motif_concordance']:
        if filters_object['transcription_factors']:
            filters += [SNP.tf_aggregated_snps.any(TranscriptionFactorSNP.motif_concordance.in_(filters_object['motif_concordance']) &
                                                  TranscriptionFactorSNP.transcription_factor.has(TranscriptionFactor.name.in_(filters_object['transcription_factors'])))]
        else:
            filters += [SNP.tf_aggregated_snps.any(TranscriptionFactorSNP.motif_concordance.in_(filters_object['motif_concordance']))]

    return filters


def get_snps_by_advanced_filters_tsv(filters_object):
    file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
    csv_writer = csv.writer(file, dialect=TsvDialect)

    headers = ['Chromosome', 'Position', 'Ref', 'Alt']
    names = dict(zip(headers, ['chromosome', 'position', 'ref', 'alt']))

    join_tuples = []
    additional_columns = []
    query_args = [SNP]
    aliases = dict()

    for what_for in ('TF', 'CL'):
        aggregation_class = {'TF': TranscriptionFactor, 'CL': CellLine}[what_for]
        aggregated_snp_class = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[what_for]
        id_field = {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
        query_args.append(db.func.group_concat(aggregation_class.name.distinct()))
        headers.append('ASB in {}'.format({'TF': 'transcription factors', 'CL': 'cell types'}[what_for]))
        join_tuples += [
            (
                aggregated_snp_class,
                (aggregated_snp_class.chromosome == SNP.chromosome) &
                (aggregated_snp_class.position == SNP.position) &
                (aggregated_snp_class.alt == SNP.alt),
            ),
            (
                aggregation_class,
                getattr(aggregation_class, id_field) == getattr(aggregated_snp_class, id_field),
            )
        ]

    for what_for in ('TF', 'CL'):
        filter_object_key = {'TF': 'transcription_factors', 'CL': 'cell_types'}[what_for]
        aggregation_class = {'TF': TranscriptionFactor, 'CL': CellLine}[what_for]
        aggregated_snp_class = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[what_for]
        id_field = {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
        if filters_object[filter_object_key]:
            for name in filters_object[filter_object_key]:
                aliases[name] = aliased(aggregated_snp_class)
                aggregation_entity = aggregation_class.query.filter_by(name=name).one_or_none()
                if not aggregation_entity:
                    continue
                aggregation_id = getattr(aggregation_entity, id_field)
                join_tuples.append(
                    (aliases[name],
                     (getattr(
                         aliases[name],
                         {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
                     ) == aggregation_id) &
                     (aliases[name].chromosome == SNP.chromosome) &
                     (aliases[name].position == SNP.position) &
                     (aliases[name].alt == SNP.alt))
                )
                for field, label in [
                    ('log_p_value_ref', '{}_FDR_Ref'.format(name)),
                    ('log_p_value_alt', '{}_FDR_Alt'.format(name)),
                    ('es_ref', '{}_Effect_Size_Ref'.format(name)),
                    ('es_alt', '{}_Effect_Size_Alt'.format(name)),
                ]:
                    headers.append(label)
                    additional_columns.append(getattr(aggregated_snp_class, field).label(label))

    found_snps = session.query(*query_args)
    found_snps = found_snps.filter(*construct_advanced_filters(filters_object))
    for cls, condition in join_tuples:
        found_snps = found_snps.join(cls, condition)
    for column in additional_columns:
        found_snps = found_snps.add_column(column)
    found_snps = found_snps.group_by(SNP)

    csv_writer.writerow(headers)
    for tup in found_snps:
        snp = tup[0]
        columns = tup[1:]
        csv_writer.writerow([getattr(snp, names[header]) for header in headers[:4]] + list(columns))

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


def get_hints_for_gene_name(in_str):
    filters = (Gene.gene_name.like(in_str),) if in_str else ()
    return Gene.query.filter(*filters).order_by(Gene.gene_name).limit(3).all()


def get_overall_statistics():
    return {
        'transcription_factors_count': TranscriptionFactor.query.filter(
            TranscriptionFactor.aggregated_snps_count > 0).count(),
        'cell_types_count': CellLine.query.filter(CellLine.aggregated_snps_count > 0).count(),
        'snps_count': db.session.query(SNP.rs_id).distinct().count(),
    }
