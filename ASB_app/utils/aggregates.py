from ASB_app import logger
from sqlalchemy_utils.aggregates import manager
from sqlalchemy import update

import numpy as np

from ASB_app.constants import db_name_property_dict, fdr_classes, fdr_choices

from ASB_app.releases import current_release
from ASB_app.models import CandidateSNP, CandidateRS
from ASB_app.utils.statistics import get_fdr_class

session = current_release.session
db = current_release.db

chunk_size = 10000


if current_release.name != 'dnase':
    TranscriptionFactorSNP, CellLineSNP, TranscriptionFactor, CellLine, Phenotype, SNP, \
    PhenotypeSNPCorrespondence, Experiment, Gene = \
        current_release.TranscriptionFactorSNP, current_release.CellLineSNP, current_release.TranscriptionFactor, \
        current_release.CellLine, current_release.Phenotype, current_release.SNP, current_release.PhenotypeSNPCorrespondence, \
        current_release.Experiment, current_release.Gene
else:
    CellLineSNP, CellLine, Phenotype, SNP, \
    PhenotypeSNPCorrespondence, Experiment, Gene = \
        current_release.CellLineSNP, \
        current_release.CellLine, current_release.Phenotype, current_release.SNP, current_release.PhenotypeSNPCorrespondence, \
        current_release.Experiment, current_release.Gene


class TsvDialect:
    delimiter = '\t'
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r\n'
    quoting = 0


def update_aggregated_fields(mappers=(CellLineSNP) if current_release.name == 'dnase' else (TranscriptionFactorSNP, CellLineSNP)):
    """
    Updates all columns with @aggregated decorator, depending on mappers
    :param mappers: list of db.Model subclasses
    :return: None
    """
    for cls in mappers:
        logger.info('Updating aggregates, depending on {}'.format(cls.__name__))
        count = cls.query.count()
        offset = 0
        max_count = chunk_size
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
        max_count = chunk_size
        while count > 0:
            print(count)
            for item in query.offset(offset).limit(max_count):
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
            offset += max_count
            count -= max_count


def update_experiments_count():
    for cls in [TranscriptionFactor, CellLine]:
        objects = []
        query = cls.query
        count = cls.query.count()
        offset = 0
        max_count = chunk_size
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
    max_count = chunk_size * 100
    while count > 0:
        print(count)
        for snp in query.offset(offset).limit(max_count):
            if snp.motif_log_p_ref:
                snp.motif_log_2_fc = (snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)
                # passes_fdr_filters = snp.best_p_value >= 1 + np.log10(2)  # 0.05
                passes_motif_filters = ((snp.motif_log_p_ref >= 3 + np.log10(2)) or
                                        (snp.motif_log_p_alt >= 3 + np.log10(2)))  # 0.0005
            else:
                assert not snp.motif_log_p_alt
                assert not snp.motif_log_2_fc
                continue
            if not passes_motif_filters:
                snp.motif_concordance = 'No Hit'
                snp.motif_log_2_fc = None
                continue

            snp.motif_log_2_fc = (snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)

            if abs((snp.motif_log_p_alt - snp.motif_log_p_ref) / np.log10(2)) >= 2:
                snp.motif_concordance = 'Concordant' if (snp.motif_log_p_alt - snp.motif_log_p_ref) * \
                                                        (snp.log_p_value_alt - snp.log_p_value_ref) > 0 \
                    else 'Discordant'
            else:
                snp.motif_concordance = 'Weak Concordant' if (snp.motif_log_p_alt - snp.motif_log_p_ref) * \
                                                             (snp.log_p_value_alt - snp.log_p_value_ref) > 0 \
                    else 'Weak Discordant'
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
    ).group_by(SNP.chromosome, SNP.position, SNP.rs_id, SNP.ref, SNP.alt, Phenotype.db_name)
    count = q.count()
    offset = 0
    max_count = chunk_size
    while count > 0:
        print(count)
        for snp, db_name in q.order_by(SNP.rs_id).limit(max_count).offset(offset):
            setattr(snp, db_name_property_dict[db_name], True)
        session.commit()
        session.close()
        offset += max_count
        count -= max_count


def update_has_concordance():
    q = session.query(SNP, db.func.coalesce(TranscriptionFactorSNP.tf_snp_id)).join(
        TranscriptionFactorSNP,
        (SNP.chromosome == TranscriptionFactorSNP.chromosome) &
        (SNP.position == TranscriptionFactorSNP.position) &
        (SNP.alt == TranscriptionFactorSNP.alt)
    ).filter(
        TranscriptionFactorSNP.motif_concordance.in_({'Concordant', 'Weak Concordant'})
    ).group_by(SNP, db.func.coalesce(TranscriptionFactorSNP.tf_snp_id))
    count = q.count()
    print(count)
    offset = 0
    max_count = chunk_size
    while count > 0:
        print(count)
        for snp, tf_snp_id in q.order_by(SNP.rs_id).limit(max_count).offset(offset):
            setattr(snp, 'has_concordance', True)
        session.commit()
        session.close()
        offset += max_count
        count -= max_count


def update_snp_best_p_value():
    q = session.query(SNP, db.func.max(TranscriptionFactorSNP.best_p_value), db.func.max(CellLineSNP.best_p_value)).join(
        TranscriptionFactorSNP,
        SNP.tf_aggregated_snps,
        isouter=True,
    ).join(
        CellLineSNP,
        SNP.cl_aggregated_snps,
        isouter=True,
    ).group_by(SNP)

    for snp, tfsnp_bp, clsnp_bp in q:
        if tfsnp_bp is None:
            snp.best_p_value = clsnp_bp
        elif clsnp_bp is None:
            snp.best_p_value = tfsnp_bp
        else:
            snp.best_p_value = max(tfsnp_bp, clsnp_bp)
    session.commit()


def update_fdr_class(model):
    table = model.__table__
    for fdr in fdr_choices[::-1]:
        print(fdr)
        session.execute(update(table).where(table.c.best_p_value >= -np.log10(float(fdr))).values(fdr_class=fdr))
        session.commit()
        session.close()


def update_all_fdr_class():
    for model in SNP, TranscriptionFactorSNP, CellLineSNP, CandidateSNP:
        print(model)
        update_fdr_class(model)


def migrate_genes():
    print(Gene.query.count())
    for i, gene in enumerate(Gene.query.all(), 1):
        if i % 10000 == 0:
            print(i)
        if gene.start_pos == 1 and gene.end_pos == 1:
            assert gene.gene_id == gene.gene_name
            gene.snps_count = None
            continue

        if gene.orientation:
            gene.start_pos += 5000
        else:
            gene.end_pos -= 5000

        if gene.orientation:
            filters = SNP.chromosome == gene.chromosome, SNP.position.between(gene.start_pos - 5000, gene.end_pos)
        else:
            filters = SNP.chromosome == gene.chromosome, SNP.position.between(gene.start_pos, gene.end_pos + 5000)

        gene.snps_count = SNP.query.filter(*filters).count()
    session.commit()


def update_gene_snps_count():
    for i, gene in enumerate(Gene.query.all(), 1):
        if i % 100 == 0:
            print(i)
        if gene.orientation:
            filters = SNP.chromosome == gene.chromosome, SNP.position.between(gene.start_pos - 5000, gene.end_pos)
        else:
            filters = SNP.chromosome == gene.chromosome, SNP.position.between(gene.start_pos, gene.end_pos + 5000)

        snps = gene.proximal_promoter_snps
        gene.snps_count = len(snps)
        gene.snps_count005 = len([x for x in gene.proximal_promoter_snps if x.fdr_class in ('0.01', '0.05')])
        # gene.snps_count = SNP.query.filter(*filters).count()
        # gene.snps_count005 = SNP.query.filter(*filters, SNP.fdr_class.in_(('0.01', '0.05'))).count()
        # gene.eqtl_snps_count = SNP.query.join(Gene, SNP.target_genes).filter(Gene.gene_id == gene.gene_id).count()
        # gene.eqtl_snps_count005 = SNP.query.filter(SNP.fdr_class.in_(('0.01', '0.05'))).join(Gene, SNP.target_genes).filter(Gene.gene_id == gene.gene_id).count()
    session.commit()


def update_best_p_value():
    pass


def update_gene_promoter_snp_correspondence():
    with open('D:\\Sashok\\Desktop\\genes_promoter_snps.sql', 'w') as f:
        f.write('INSERT INTO adastra_dan.genes_promoter_snps (chromosome, position, alt, pair_id, gene_id) VALUES\n')
        ai = 1
        for i, gene in enumerate(Gene.query.all(), 1):
            if i == 1:
                session.close()
            if i % 100 == 0:
                print(i)
            if gene.orientation:
                filters = SNP.chromosome == gene.chromosome, SNP.position.between(gene.start_pos - 5000, gene.end_pos)
            else:
                filters = SNP.chromosome == gene.chromosome, SNP.position.between(gene.start_pos, gene.end_pos + 5000)
            snps = session.query(SNP.chromosome, SNP.position, SNP.alt).filter(*filters).all()
            if snps:
                gene.snps_count = len(snps)
                for k, (c, p, a) in enumerate(snps, 1):
                    if i != 1 or k != 1:
                        f.write(', ')
                    # gene.proximal_promoter_snps = snps
                    f.write("('{}', {}, '{}', {}, '{}')".format(c, p, a, ai, gene.gene_id))
                    ai += 1
            if i % 10000 == 0:
                session.close()
        f.write(';')


def update_all():
    update_phenotype_associations()
    update_has_concordance()
