from ASB_app import logger
from sqlalchemy_utils.aggregates import manager
from sqlalchemy import update, case

import numpy as np

from ASB_app.constants import db_name_property_dict, fdr_classes, fdr_choices, es_choices, es_classes

from ASB_app.releases import current_release
from ASB_app.models import CandidateSNP, CandidateRS
from ASB_app.utils.statistics import get_fdr_class, get_es_class

session = current_release.session
db = current_release.db

chunk_size = 10000

FaireSNP, Faire, DnaseSNP, Dnase, AtacSNP, Atac, Phenotype, SNP, \
PhenotypeSNPCorrespondence, Experiment, Gene = current_release.FaireSNP, current_release.Faire, \
    current_release.DnaseSNP, current_release.Dnase, \
    current_release.AtacSNP, current_release.Atac, \
    current_release.Phenotype, current_release.SNP, current_release.PhenotypeSNPCorrespondence, \
    current_release.Experiment, current_release.Gene


class TsvDialect:
    delimiter = '\t'
    quotechar = '"'
    escapechar = None
    doublequote = True
    skipinitialspace = False
    lineterminator = '\r\n'
    quoting = 0


def update_aggregated_fields(
        mappers=(AtacSNP, FaireSNP, DnaseSNP)):
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
            objects = cls.query.order_by({AtacSNP: AtacSNP.cl_id,
                                          DnaseSNP: DnaseSNP.cl_id,
                                          FaireSNP: FaireSNP.cl_id,
                                          Phenotype: Phenotype.phenotype_id}
                                         [cls]).offset(offset).limit(max_count).all()
            manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
            session.commit()
            session.close()
            offset += max_count
            count -= max_count


def update_aggregated_snp_count():
    for snp_ne_obyect, cls in [(AtacSNP, Atac), (FaireSNP, Faire), (DnaseSNP, Dnase)]:
        objects = []
        query = cls.query
        count = cls.query.count()
        offset = 0
        max_count = chunk_size
        while count > 0:
            print(count)
            for item in query.offset(offset).limit(max_count):
                objects.append(snp_ne_obyect.query.filter(
                    snp_ne_obyect.cl_id == item.cl_id
                ).first())
            manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
            session.commit()
            session.close()
            offset += max_count
            count -= max_count


def update_experiments_count():
    for cls in [AtacSNP, FaireSNP, DnaseSNP]:
        objects = []
        query = cls.query
        count = cls.query.count()
        offset = 0
        max_count = chunk_size
        while count > 0:
            for item in query.offset(offset).limit(max_count):
                objects.append(Experiment.query.filter(
                    Experiment.cl_id == item.cl_id
                ).first())
            manager.construct_aggregate_queries(session, ...)  # Второй параметр не используется
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
            db = db_name if db_name != 'qtltissues' else 'qtl'
            setattr(snp, db_name_property_dict[db], True)
        session.commit()
        session.close()
        offset += max_count
        count -= max_count


def update_snp_best_p_value():
    q = session.query(SNP, db.func.max(DnaseSNP.best_p_value),
                      db.func.max(AtacSNP.best_p_value), db.func.max(FaireSNP.best_p_value)).join(
        DnaseSNP,
        SNP.tf_aggregated_snps,
        isouter=True,
    ).join(
        FaireSNP,
        SNP.cl_aggregated_snps,
        isouter=True,
    ).join(
        AtacSNP,
        SNP.cl_aggregated_snps,
        isouter=True,
    ).group_by(SNP)

    for snp, dnasebp, atacbp, fairebp in q:
        for c in (dnasebp, atacbp, fairebp):
            if c is not None:
                snp.best_p_value = c
                break
    session.commit()


def update_fdr_class(model):
    table = model.__table__
    session.execute(
        update(table).values(
            fdr_class=case(
                [(table.c.best_p_value >= -np.log10(float(fdr)), fdr)
                 for fdr in fdr_choices],
                else_=fdr_classes[-1]
            )
        )
    )
    session.commit()
    session.close()


def update_es_class(model):
    table = model.__table__
    session.execute(
        update(table).values(
            es_class=case(
                [(table.c.best_es >= float(es), es)
                 for es in es_choices],
                else_=es_classes[-1]
            )
        )
    )
    session.commit()
    session.close()


def update_all_fdr_class():
    for model in SNP, TranscriptionFactorSNP, CellLineSNP, CandidateSNP:
        print(model)
        update_fdr_class(model)


def update_all_es_class():
    for model in SNP, TranscriptionFactorSNP, CellLineSNP, CandidateSNP:
        print(model)
        update_es_class(model)


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
        # gene.snps_count005 = len([x for x in gene.proximal_promoter_snps if x.fdr_class in ('0.01', '0.05')])
        # gene.snps_count = SNP.query.filter(*filters).count()
        # gene.snps_count005 = SNP.query.filter(*filters, SNP.fdr_class.in_(('0.01', '0.05'))).count()
        # gene.eqtl_snps_count = SNP.query.join(Gene, SNP.target_genes).filter(Gene.gene_id == gene.gene_id).count()
        # gene.eqtl_snps_count005 = SNP.query.filter(SNP.fdr_class.in_(('0.01', '0.05'))).join(Gene, SNP.target_genes).filter(Gene.gene_id == gene.gene_id).count()
    session.commit()


def update_best_p_value():
    q = session.query(SNP, db.func.greatest(
            db.func.coalesce(db.func.max(AtacSNP.best_p_value), 0),
            db.func.coalesce(db.func.max(DnaseSNP.best_p_value), 0),
            db.func.coalesce(db.func.max(FaireSNP.best_p_value), 0)
    )) \
        .join(FaireSNP, SNP.faire_aggregated_snps, isouter=True) \
        .join(AtacSNP, SNP.atac_aggregated_snps, isouter=True) \
        .join(DnaseSNP, SNP.dnase_aggregated_snps, isouter=True) \
        .group_by(SNP)
    for i, (snp, best_p) in enumerate(q, 1):
        snp.best_p_value = best_p
        snp.fdr_class = get_fdr_class(best_p)
    session.commit()


def update_best_es():
    q = session.query(SNP, db.func.greatest(
        db.func.coalesce(db.func.max(AtacSNP.best_es), -1000),
        db.func.coalesce(db.func.max(DnaseSNP.best_es), -1000),
        db.func.coalesce(db.func.max(FaireSNP.best_es), -1000),
    )) \
        .join(FaireSNP, SNP.faire_aggregated_snps, isouter=True) \
        .join(AtacSNP, SNP.atac_aggregated_snps, isouter=True) \
        .join(DnaseSNP, SNP.dnase_aggregated_snps, isouter=True) \
        .group_by(SNP)
    for i, (snp, best_es) in enumerate(q, 1):
        snp.best_es = None if best_es == -1000 else best_es
        snp.es_class = get_es_class(best_es)
    session.commit()


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
