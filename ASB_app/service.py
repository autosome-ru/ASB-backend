from ASB_app import session, logger, db
from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, SNP, TranscriptionFactor, CellLine, Phenotype
from ASB_app.exceptions import ParsingError
from sqlalchemy import not_
import csv
import tempfile
from flask import send_file
from ASB_app.utils import TsvDialect


def get_snps_by_rs_id(rs_id):
    return SNP.query.filter(SNP.rs_id == rs_id).all()


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


def get_snps_by_genome_position(chr, pos1, pos2):
    return SNP.query.filter(SNP.chromosome == chr, SNP.position.between(pos1, pos2)).all()


def get_snps_by_advanced_filters_or(filters_object):
    if filters_object['transcription_factors']:
        tf_filters = (SNP.tf_aggregated_snps.any(TranscriptionFactorSNP.tf_id.in_(
            [getattr(TranscriptionFactor.query.filter(TranscriptionFactor.name == tf_name).one_or_none(),
                     'tf_id', None)
             for tf_name in filters_object['transcription_factors']])),)
    else:
        tf_filters = ()

    if filters_object['cell_types']:
        cl_filters = (SNP.cl_aggregated_snps.any(CellLineSNP.cl_id.in_(
            [getattr(CellLine.query.filter(CellLine.name == cl_name).one_or_none(),
                     'cl_id', None)
             for cl_name in filters_object['cell_types']])),)
    else:
        cl_filters = ()

    if filters_object['chromosome']:
        if not filters_object['start'] or not filters_object['end']:
            raise ParsingError
        chrpos_filters = (SNP.chromosome == filters_object['chromosome'],
                          SNP.position.between(filters_object['start'], filters_object['end']))
    else:
        chrpos_filters = ()

    return SNP.query.filter(*(tf_filters + cl_filters + chrpos_filters)).all()


def get_snps_by_advanced_filters(filters_object):
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

    for phenotype_db in filters_object['phenotype_databases']:
        filters += [SNP.phenotypes.any(Phenotype.db_name == phenotype_db)]

    return SNP.query.filter(*filters).all()


def get_snps_by_advanced_filters_tsv(filters_object):
    found_snps = get_snps_by_advanced_filters(filters_object)

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


def get_hints(what_for, in_str, used_options):
    cls = {'TF': TranscriptionFactor, 'CL': CellLine}[what_for]
    filters = (((cls.name.like(in_str),) if in_str else ()) +
               ((not_(cls.name.in_(used_options)),) if used_options else ()) +
               (cls.aggregated_snps_count,))
    return cls.query.filter(*filters).order_by(cls.aggregated_snps_count.desc()).limit(3).all()
