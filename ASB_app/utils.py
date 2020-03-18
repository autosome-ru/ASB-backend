from ASB_app import session, logger
from sqlalchemy_utils.aggregates import manager

from ASB_app.models import TranscriptionFactorSNP, CellLineSNP


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
