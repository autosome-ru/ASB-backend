from ASB_app import session, logger
from sqlalchemy_utils.aggregates import manager

import numpy as np

from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, TranscriptionFactor, CellLine


class TsvDialect:
    delimiter = '\t'
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r\n'
    quoting = 0


def update_all_aggregated_fields():
    """
    Updates all columns with @aggregated decorator
    :return: None
    """
    for cls in [TranscriptionFactorSNP, CellLineSNP]:
        logger.info('Updating aggregates, depending on {}'.format(cls.__name__))
        count = cls.query.count()
        offset = 0
        max_count = 999
        while count > max_count:
            print(count)
            objects = cls.query.order_by({TranscriptionFactorSNP: TranscriptionFactorSNP.tf_id,
                                          CellLineSNP: CellLineSNP.cl_id}
                                         [cls]).offset(offset).limit(max_count).all()
            manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
            session.commit()
            session.close()
            offset += max_count
            count -= max_count
        objects = cls.query.offset(offset).all()
        manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
        session.commit()
        session.close()


def update_aggregated_snp_count():
    objects = []
    for cls in [TranscriptionFactor, CellLine]:
        for item in cls.query:
            print(item.name)
            if cls == TranscriptionFactor:
                objects.append(TranscriptionFactorSNP.query.filter(
                    TranscriptionFactorSNP.tf_id == item.tf_id
                ).first())
            elif cls == CellLine:
                objects.append(CellLineSNP.query.filter(
                    CellLineSNP.cl_id == item.cl_id
                ).first())
        manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
        session.commit()
        session.close()


def update_motif_concordance():
    query = TranscriptionFactorSNP.query
    count = query.count()
    offset = 0
    max_count = 999
    while count > 0:
        print(count)
        for snp in query.offset(offset).limit(max_count):
            if snp.motif_log_p_ref:
                snp.motif_log_2_fc = (snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)
            passes_filters = (snp.best_p_value >= 1 + np.log10(2)  # 0.05
                              and (abs((snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)) >= 4
                                   if snp.motif_log_p_ref
                                   else False))
            if not passes_filters:
                snp.motif_concordance = None
                continue

            snp.motif_log_2_fc = (snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)
            snp.motif_concordance = (snp.motif_log_p_alt - snp.motif_log_p_ref) * \
                                    (snp.log_p_value_alt - snp.log_p_value_ref) > 0
        session.commit()
        session.close()
        offset += max_count
        count -= max_count
