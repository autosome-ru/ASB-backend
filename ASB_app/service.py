from ASB_app import session, logger
from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, SNP


def get_snps_by_rs(rs_id, param):
    Class = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[param]

    snps = (session.query(SNP.rs_id, Class)
            .filter(SNP.rs_id == rs_id)
            .join(Class)).all()

    if not snps:
        return []

    print(snps)

    snp_ids, ret_snps = zip(*snps)
    return list(ret_snps)


def get_full_snps_by_rs(rs_id):
    return SNP.query.filter(SNP.rs_id == rs_id).all()
