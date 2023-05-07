from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy_utils.aggregates import aggregated

from ASB_app.constants import chromosomes, nucleotides, bads, fdr_classes, es_classes

from ASB_app.releases import Release
from .placeholders import abstract_models

import numpy as np

__all__ = []

for release in Release.__subclasses__():
    db = release.db

    class Faire(db.Model):
        __tablename__ = 'faire'
        __bind_key__ = release.name

        cl_id = db.Column(db.Integer, primary_key=True)
        name = db.Column(db.String(200), nullable=False, index=True)

        non_input_experiments = db.relationship(
            'Experiment',
            primaryjoin='Experiment.faire_id == Faire.cl_id'
        )

        @aggregated('cl_aggregated_snps', db.Column(db.Integer))
        def aggregated_snps_count(self):
            return db.func.count(FaireSNP.cl_snp_id)

        @aggregated('non_input_experiments', db.Column(db.Integer))
        def experiments_count(self):
            return db.func.count(Experiment.exp_id)

        cl_aggregated_snps005 = db.relationship(
            'FaireSNP',
            primaryjoin='(Faire.cl_id == FaireSNP.cl_id) & (FaireSNP.best_p_value >= {})'.format(
                -np.log10(0.05))
        )

        @aggregated('cl_aggregated_snps005', db.Column(db.Integer))
        def aggregated_snps_count005(self):
            return db.func.count(FaireSNP.cl_snp_id)

        cl_aggregated_snps010 = db.relationship(
            'FaireSNP',
            primaryjoin='(Faire.cl_id == FaireSNP.cl_id) & (FaireSNP.best_p_value >= {})'.format(
                -np.log10(0.1))
        )

        @aggregated('cl_aggregated_snps010', db.Column(db.Integer))
        def aggregated_snps_count010(self):
            return db.func.count(FaireSNP.cl_snp_id)

        def __repr__(self):
            return '<CellLine #{0.cl_id}, {0.name}>'.format(self)

    class Atac(db.Model):
        __tablename__ = 'atac'
        __bind_key__ = release.name

        cl_id = db.Column(db.Integer, primary_key=True)
        name = db.Column(db.String(200), nullable=False, index=True)

        non_input_experiments = db.relationship(
            'Experiment',
            primaryjoin='Experiment.atac_id == Atac.cl_id'
        )

        @aggregated('cl_aggregated_snps', db.Column(db.Integer))
        def aggregated_snps_count(self):
            return db.func.count(AtacSNP.cl_snp_id)

        @aggregated('non_input_experiments', db.Column(db.Integer))
        def experiments_count(self):
            return db.func.count(Experiment.exp_id)

        cl_aggregated_snps005 = db.relationship(
            'AtacSNP',
            primaryjoin='(Atac.cl_id == AtacSNP.cl_id) & (AtacSNP.best_p_value >= {})'.format(
                -np.log10(0.05))
        )

        @aggregated('cl_aggregated_snps005', db.Column(db.Integer))
        def aggregated_snps_count005(self):
            return db.func.count(AtacSNP.cl_snp_id)

        cl_aggregated_snps010 = db.relationship(
            'AtacSNP',
            primaryjoin='(Atac.cl_id == AtacSNP.cl_id) & (AtacSNP.best_p_value >= {})'.format(
                -np.log10(0.1))
        )

        @aggregated('cl_aggregated_snps010', db.Column(db.Integer))
        def aggregated_snps_count010(self):
            return db.func.count(AtacSNP.cl_snp_id)

        def __repr__(self):
            return '<CellLine #{0.cl_id}, {0.name}>'.format(self)

    class Dnase(db.Model):
        __tablename__ = 'dnase'
        __bind_key__ = release.name

        cl_id = db.Column(db.Integer, primary_key=True)
        name = db.Column(db.String(200), nullable=False, index=True)

        non_input_experiments = db.relationship(
            'Experiment',
            primaryjoin='Experiment.dnase_id == Dnase.cl_id'
        )

        @aggregated('cl_aggregated_snps', db.Column(db.Integer))
        def aggregated_snps_count(self):
            return db.func.count(DnaseSNP.cl_snp_id)

        @aggregated('non_input_experiments', db.Column(db.Integer))
        def experiments_count(self):
            return db.func.count(Experiment.exp_id)

        cl_aggregated_snps005 = db.relationship(
            'DnaseSNP',
            primaryjoin='(Dnase.cl_id == DnaseSNP.cl_id) & (DnaseSNP.best_p_value >= {})'.format(
                -np.log10(0.05))
        )

        @aggregated('cl_aggregated_snps005', db.Column(db.Integer))
        def aggregated_snps_count005(self):
            return db.func.count(DnaseSNP.cl_snp_id)

        cl_aggregated_snps010 = db.relationship(
            'DnaseSNP',
            primaryjoin='(Dnase.cl_id == DnaseSNP.cl_id) & (DnaseSNP.best_p_value >= {})'.format(
                -np.log10(0.1))
        )

        @aggregated('cl_aggregated_snps010', db.Column(db.Integer))
        def aggregated_snps_count010(self):
            return db.func.count(DnaseSNP.cl_snp_id)

        def __repr__(self):
            return '<CellLine #{0.cl_id}, {0.name}>'.format(self)

    class Experiment(db.Model):
        __tablename__ = 'experiments'
        __bind_key__ = release.name
        __table_args__ = (
            db.Index('align_index', 'align'),
        )
        exp_id = db.Column(db.String(10), primary_key=True)
        align = db.Column(db.String(13), nullable=False)

        dnase_id = db.Column(db.Integer, db.ForeignKey('Dnase.cl_id'))
        atac_id = db.Column(db.Integer, db.ForeignKey('Atac.cl_id'))
        faire_id = db.Column(db.Integer, db.ForeignKey('Faire.cl_id'))
        geo_gse = db.Column(db.String(26))
        encode = db.Column(db.String(30))

        bad_group_id = db.Column(db.Integer, db.ForeignKey('bad_groups.bad_group_id'))

        bad_group = db.relationship('BADGroup', backref='experiments')
        cell_line = db.relationship('CellLine', backref='experiments')

        def __repr__(self):
            return '<Experiment #{0.exp_id}>'.format(self)


    class BADGroup(db.Model):
        __tablename__ = 'bad_groups'
        __bind_key__ = release.name

        bad_group_id = db.Column(db.Integer, primary_key=True)
        bad_group_name = db.Column(db.String(200), nullable=False, unique=True)

        def __repr__(self):
            return '<BAD Group #{0.bad_group_id}: {0.bad_group_name}>'.format(self)


    class ExpSNP(db.Model):
        __tablename__ = 'exp_snps'
        __bind_key__ = release.name

        exp_snp_id = db.Column(db.Integer, primary_key=True)
        ref_readcount = db.Column(db.Integer, nullable=False)
        alt_readcount = db.Column(db.Integer, nullable=False)
        p_value_ref = db.Column(db.Float)
        p_value_alt = db.Column(db.Float)
        bad = db.Column(db.Enum(*bads))
        atac_snp_id = db.Column(db.Integer, db.ForeignKey('atac_snps.cl_snp_id'))
        exp_id = db.Column(db.String(10), db.ForeignKey('experiments.exp_id'), nullable=False, index=True)

        atac_aggregated_snp = db.relationship('AtacSNP', backref='exp_snps_atac')
        dnase_aggregated_snp = db.relationship('DnaseSNP', backref='exp_snps_dnase')
        faire_aggregated_snp = db.relationship('FaireSNP', backref='exp_snps_faire')
        experiment = db.relationship('Experiment', backref='exp_snps')

        atac = db.relationship(
            'Atac',
            secondary='atac_snps',
            backref='exp_snps_atac'
        )
        dnase = db.relationship(
            'Dnase',
            secondary='dnase_snps',
            backref='exp_snps_dnase'
        )
        faire = db.relationship(
            'Faire',
            secondary='faire_snps',
            backref='exp_snps_faire'
        )

        def __repr__(self):
            return '<ExpSNP #{0.exp_snp_id}, {0.exp_id}, {0.tf_snp_id}, {0.cl_snp_id}>'.format(self)


    class GenomePolymorphismLocation(db.Model):
        __abstract__ = True

        chromosome = db.Column(db.Enum(*chromosomes), nullable=False)
        position = db.Column(db.Integer, nullable=False)
        alt = db.Column(db.Enum(*nucleotides), nullable=False)


    class SNP(GenomePolymorphismLocation):
        __tablename__ = 'snps'
        __bind_key__ = release.name
        __table_args__ = (
            db.PrimaryKeyConstraint('chromosome', 'position', 'alt'),
            db.Index('rs_index', 'rs_id'),
        )

        ref = db.Column(db.Enum(*nucleotides), nullable=False)
        rs_id = db.Column(db.Integer, nullable=False)

        cl_aggregated_snps = db.relationship('CellLineSNP',
                                             order_by='CellLineSNP.best_p_value.desc()',
                                             back_populates='snp')

        has_clinvar_associations = db.Column(db.Boolean)
        has_phewas_associations = db.Column(db.Boolean)
        has_ebi_associations = db.Column(db.Boolean)
        has_qtl_associations = db.Column(db.Boolean)
        has_grasp_associations = db.Column(db.Boolean)
        has_finemapping_associations = db.Column(db.Boolean)
        best_p_value = db.Column(db.Float, index=True)
        fdr_class = db.Column(db.Enum(*fdr_classes), index=True)
        best_es = db.Column(db.Float, index=True)
        es_class = db.Column(db.Enum(*es_classes), index=True)

        context = db.Column(db.String(51))

        has_concordance = db.Column(db.Boolean, index=True, server_default='0', nullable=False)

        def __repr__(self):
            return '<SNP rs{0.rs_id}, {0.chromosome}, {0.position}, {0.ref}, {0.alt}>'.format(self)


    class AggregatedSNP(GenomePolymorphismLocation):
        __abstract__ = True

        log_p_value_ref = db.Column(db.Float)
        log_p_value_alt = db.Column(db.Float)
        is_asb = db.Column(db.Boolean, nullable=False)
        es_ref = db.Column(db.Float)
        es_alt = db.Column(db.Float)
        mean_bad = db.Column(db.Float)
        peak_calls = db.Column(db.Integer)
        peak_callers = db.Column(db.Integer)
        best_p_value = db.Column(db.Float)
        fdr_class = db.Column(db.Enum(*fdr_classes))
        best_es = db.Column(db.Float)
        es_class = db.Column(db.Enum(*es_classes))

    class FaireSNP(AggregatedSNP):
        __tablename__ = 'faire_snps'
        __bind_key__ = release.name
        __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                                  ['snps.chromosome', 'snps.position', 'snps.alt']),
                          db.UniqueConstraint('chromosome', 'position', 'alt', 'cl_id',
                                              name='faire_unique_mutation'),
                          db.Index('faire_id_index', 'cl_id'),
                          db.Index('ix_faire_snps_best_p_value', 'best_p_value'),
                          db.Index('ix_faire_snps_fdr_class', 'fdr_class'),
                          db.Index('ix_faire_snps_best_es', 'best_es'),
                          db.Index('ix_faire_snps_es_class', 'es_class'),
                          )

        cl_snp_id = db.Column(db.Integer, primary_key=True)
        cl_id = db.Column(db.Integer, db.ForeignKey('faire.cl_id'), nullable=False)

        snp = db.relationship('SNP', back_populates='faire_aggregated_snps')
        cell_line = db.relationship('Faire', backref='faire_aggregated_snps')

        def __repr__(self):
            return '<CellLineSNP #{0.cl_snp_id} at {0.chromosome} {0.position} {0.alt}>'.format(self)

    class DnaseSNP(AggregatedSNP):
        __tablename__ = 'dnase_snps'
        __bind_key__ = release.name
        __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                                  ['snps.chromosome', 'snps.position', 'snps.alt']),
                          db.UniqueConstraint('chromosome', 'position', 'alt', 'cl_id',
                                              name='dnase_unique_mutation'),
                          db.Index('dnase_id_index', 'cl_id'),
                          db.Index('ix_dnase_snps_best_p_value', 'best_p_value'),
                          db.Index('ix_dnase_snps_fdr_class', 'fdr_class'),
                          db.Index('ix_dnase_snps_best_es', 'best_es'),
                          db.Index('ix_dnase_snps_es_class', 'es_class'),
                          )

        cl_snp_id = db.Column(db.Integer, primary_key=True)
        cl_id = db.Column(db.Integer, db.ForeignKey('dnase.cl_id'), nullable=False)

        snp = db.relationship('SNP', back_populates='dnase_aggregated_snps')
        cell_line = db.relationship('Dnase', backref='dnase_aggregated_snps')

        def __repr__(self):
            return '<CellLineSNP #{0.cl_snp_id} at {0.chromosome} {0.position} {0.alt}>'.format(self)

    class AtacSNP(AggregatedSNP):
        __tablename__ = 'atac_snps'
        __bind_key__ = release.name
        __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                                  ['snps.chromosome', 'snps.position', 'snps.alt']),
                          db.UniqueConstraint('chromosome', 'position', 'alt', 'cl_id',
                                              name='atac_unique_mutation'),
                          db.Index('atac_id_index', 'cl_id'),
                          db.Index('ix_atac_snps_best_p_value', 'best_p_value'),
                          db.Index('ix_atac_snps_fdr_class', 'fdr_class'),
                          db.Index('ix_atac_snps_best_es', 'best_es'),
                          db.Index('ix_atac_snps_es_class', 'es_class'),
                          )

        cl_snp_id = db.Column(db.Integer, primary_key=True)
        cl_id = db.Column(db.Integer, db.ForeignKey('atac.cl_id'), nullable=False)

        snp = db.relationship('SNP', back_populates='atac_aggregated_snps')
        cell_line = db.relationship('Atac', backref='atac_aggregated_snps')

        def __repr__(self):
            return '<CellLineSNP #{0.cl_snp_id} at {0.chromosome} {0.position} {0.alt}>'.format(self)


    class Phenotype(db.Model):
        __tablename__ = 'phenotypes'
        __bind_key__ = release.name
        __tableargs__ = (db.UniqueConstraint('phenotype_name', 'db_name',
                                             name='unique phenotype in db'),
                         )

        phenotype_id = db.Column(db.Integer, primary_key=True)
        db_name = db.Column(db.String(100), index=True)
        phenotype_name = db.Column(db.String(200))

        snps = db.relationship(
            'SNP',
            secondary='phenotypes_SNPs',
            backref='phenotypes'
        )


    class PhenotypeSNPCorrespondence(GenomePolymorphismLocation):
        __tablename__ = 'phenotypes_SNPs'
        __bind_key__ = release.name
        __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                                  ['snps.chromosome', 'snps.position', 'snps.alt']),
                          db.UniqueConstraint('chromosome', 'position', 'alt', 'phenotype_id',
                                              name='unique_phenotype_snp_pair'),
                          db.Index('snp_index', 'chromosome', 'position', 'alt'),
                          db.Index('phenotype_index', 'phenotype_id')
                          )

        pair_id = db.Column(db.Integer, primary_key=True)
        phenotype_id = db.Column(db.Integer, db.ForeignKey('phenotypes.phenotype_id'),
                                 nullable=False)


    class Gene(db.Model):
        __tablename__ = 'genes'
        __bind_key__ = release.name

        gene_id = db.Column(db.String(30), primary_key=True)
        gene_name = db.Column(db.String(30), index=True, nullable=False)
        chromosome = db.Column(db.Enum(*chromosomes), nullable=False)
        start_pos = db.Column(db.Integer, nullable=False)
        end_pos = db.Column(db.Integer, nullable=False)
        orientation = db.Column(db.Boolean, nullable=False)
        if int(release.version) >= 3:
            snps_count = db.Column(db.Integer)
            snps_count010 = db.Column(db.Integer)
            eqtl_snps_count = db.Column(db.Integer)
            eqtl_snps_count010 = db.Column(db.Integer)
        else:
            snps_count = db.Column(db.Integer)

        snps_by_target = db.relationship(
            'SNP',
            secondary='genes_SNPs',
            backref='target_genes'
        )

        proximal_promoter_snps = db.relationship(
            'SNP',
            secondary='genes_promoter_SNPs',
            backref='genes_by_proximal_promoter'
        )


    class GeneSNPCorrespondence(GenomePolymorphismLocation):
        __tablename__ = 'genes_SNPs'
        __bind_key__ = release.name
        __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                                  ['snps.chromosome', 'snps.position', 'snps.alt']),
                          db.UniqueConstraint('chromosome', 'position', 'alt', 'gene_id',
                                              name='unique_gene_snp_pair'),
                          db.Index('snp_genes_SNPs_index', 'chromosome', 'position', 'alt'),
                          db.Index('gene_index', 'gene_id')
                          )

        pair_id = db.Column(db.Integer, primary_key=True)
        gene_id = db.Column(db.String(30), db.ForeignKey('genes.gene_id'),
                            nullable=False)


    class ProximalPromoterSNPCorrespondence(GenomePolymorphismLocation):
        __tablename__ = 'genes_promoter_SNPs'
        __bind_key__ = release.name
        __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                                  ['snps.chromosome', 'snps.position', 'snps.alt']),
                          db.UniqueConstraint('chromosome', 'position', 'alt', 'gene_id',
                                              name='unique_gene_promoter_snp_pair'),
                          db.Index('snp_genes_promoter_SNPs_index', 'chromosome', 'position', 'alt'),
                          db.Index('gene_promoter_index', 'gene_id')
                          )

        pair_id = db.Column(db.Integer, primary_key=True)
        gene_id = db.Column(db.String(30), db.ForeignKey('genes.gene_id'),
                            nullable=False)


    models = [
        Faire,
        FaireSNP,
        Dnase,
        DnaseSNP,
        Atac,
        AtacSNP,
        Experiment,
        ExpSNP,
        SNP,
        Phenotype,
        PhenotypeSNPCorrespondence,
        BADGroup,
        Gene,
        GeneSNPCorrespondence,
    ]
    for abstract_model, model in zip(abstract_models, models):
        setattr(release, abstract_model.__name__, model)
