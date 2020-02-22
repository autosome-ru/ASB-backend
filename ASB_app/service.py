from ASB_app import session, logger, db
from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, SNP


def get_snps_by_rs_id(rs_id):
    return SNP.query.filter(SNP.rs_id == rs_id).all()


def get_full_snp(rs_id, alt):
    return SNP.query.filter(
        (SNP.rs_id == rs_id) &
        (SNP.alt == alt)
    ).one_or_none()
