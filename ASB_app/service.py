from ASB_app import session, logger, db
from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, SNP, TranscriptionFactor, CellLine
from ASB_app.exceptions import ParsingError


def get_snps_by_rs_id(rs_id):
    return SNP.query.filter(SNP.rs_id == rs_id).all()


def get_full_snp(rs_id, alt):
    return SNP.query.filter(
        (SNP.rs_id == rs_id) &
        (SNP.alt == alt)
    ).one()


def get_snps_by_genome_position(chr, pos1, pos2):
    return SNP.query.filter(SNP.chromosome == chr, SNP.position.between(pos1, pos2)).all()


def get_snps_by_advanced_filters(filters_object):
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
