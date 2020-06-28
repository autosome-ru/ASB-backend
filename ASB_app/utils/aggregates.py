from sqlalchemy import distinct

from ASB_app import session, logger
from sqlalchemy_utils.aggregates import manager

import numpy as np

from ASB_app.constants import db_name_property_dict
from ASB_app.models import TranscriptionFactorSNP, CellLineSNP, TranscriptionFactor, CellLine, Phenotype, SNP, \
    PhenotypeSNPCorrespondence, Experiment


class TsvDialect:
    delimiter = '\t'
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r\n'
    quoting = 0


def update_aggregated_fields(mappers=(TranscriptionFactorSNP, CellLineSNP)):
    """
    Updates all columns with @aggregated decorator, depending on mappers
    :param mappers: list of db.Model subclasses
    :return: None
    """
    for cls in mappers:
        logger.info('Updating aggregates, depending on {}'.format(cls.__name__))
        count = cls.query.count()
        offset = 0
        max_count = 999
        while count > 0:
            print(count)
            objects = cls.query.order_by({TranscriptionFactorSNP: TranscriptionFactorSNP.tf_id,
                                          CellLineSNP: CellLineSNP.cl_id,
                                          Phenotype: Phenotype.phenotype_id}
                                         [cls]).offset(offset).limit(max_count).all()
            manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
            session.commit()
            session.close()
            offset += max_count
            count -= max_count


def update_aggregated_snp_count():
    for cls in [TranscriptionFactor, CellLine]:
        objects = []
        query = cls.query
        count = cls.query.count()
        offset = 0
        max_count = 999
        while count > 0:
            for item in query.offset(offset).limit(max_count):
                if cls == TranscriptionFactor:
                    objects.append(TranscriptionFactorSNP.query.filter(
                        TranscriptionFactorSNP.tf_id == item.tf_id
                    ).first())
                    print(objects[-1])
                elif cls == CellLine:
                    objects.append(CellLineSNP.query.filter(
                        CellLineSNP.cl_id == item.cl_id
                    ).first())
            manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
            session.commit()
            session.close()
            offset += max_count
            count -= max_count


def update_experiments_count():
    for cls in [TranscriptionFactor, CellLine]:
        objects = []
        query = cls.query
        count = cls.query.count()
        offset = 0
        max_count = 999
        while count > 0:
            for item in query.offset(offset).limit(max_count):
                if cls == TranscriptionFactor:
                    objects.append(Experiment.query.filter(
                        Experiment.tf_id == item.tf_id
                    ).first())
                    print(objects[-1])
                elif cls == CellLine:
                    objects.append(Experiment.query.filter(
                        Experiment.cl_id == item.cl_id
                    ).first())
            manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
            session.commit()
            session.close()
            offset += max_count
            count -= max_count


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
                passes_fdr_filters = snp.best_p_value >= 1 + np.log10(2)  # 0.05
                passes_motif_filters = snp.motif_log_p_ref >= 3 + np.log10(2)  # 0.0005
            else:
                assert not snp.motif_log_p_alt
                assert not snp.motif_log_2_fc
                continue
            if not passes_motif_filters:
                snp.motif_concordance = 'No Hit'
                snp.motif_log_2_fc = None
                continue

            snp.motif_log_2_fc = (snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)

            if passes_fdr_filters:
                if abs((snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)) >= 2:
                    snp.motif_concordance = 'Concordant' if (snp.motif_log_p_alt - snp.motif_log_p_ref) * \
                                                            (snp.log_p_value_alt - snp.log_p_value_ref) > 0 \
                        else 'Discordant'
                else:
                    snp.motif_concordance = 'Weak Concordant' if (snp.motif_log_p_alt - snp.motif_log_p_ref) * \
                                                            (snp.log_p_value_alt - snp.log_p_value_ref) > 0 \
                        else 'Weak Discordant'
            else:
                snp.motif_concordance = None
        session.commit()
        session.close()
        offset += max_count
        count -= max_count


def update_phenotype_associations():
    q = session.query(SNP, Phenotype.db_name).join(
        PhenotypeSNPCorrespondence,
        (SNP.chromosome == PhenotypeSNPCorrespondence.chromosome) &
        (SNP.position == PhenotypeSNPCorrespondence.position) &
        (SNP.alt == PhenotypeSNPCorrespondence.alt)
    ).join(
        Phenotype,
        PhenotypeSNPCorrespondence.phenotype_id == Phenotype.phenotype_id
    ).group_by(SNP.rs_id, SNP.alt, Phenotype.db_name)
    count = q.count()
    offset = 0
    max_count = 999
    while count > 0:
        print(count)
        for snp, db_name in q.order_by(SNP.rs_id).limit(max_count).offset(offset):
            setattr(snp, db_name_property_dict[db_name], True)
        session.commit()
        session.close()
        offset += max_count
        count -= max_count
