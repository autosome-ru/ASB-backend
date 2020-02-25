from sqlalchemy.ext.hybrid import hybrid_property

from ASB_app import db
from ASB_app.constants import chromosomes, nucleotides, bads


class TranscriptionFactor(db.Model):
    __tablename__ = 'transcription_factors'
    __table_args__ = (
        db.Index('tf_uniprot_ac_index', 'name'),
    )

    tf_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), nullable=False)
    motif_legnth = db.Column(db.Integer)

    def __repr__(self):
        return '<TranscriptionFactor #{0.tf_id}, {0.name}>'.format(self)


class CellLine(db.Model):
    __tablename__ = 'cell_lines'
    __table_args__ = (
        db.Index('cell_title_index', 'name'),
    )

    cl_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), nullable=False)

    def __repr__(self):
        return '<CellLine #{0.cl_id}, {0.name}>'.format(self)


class Experiment(db.Model):
    __tablename__ = 'experiments'
    __table_args__ = (
        db.Index('align_index', 'align'),
    )

    exp_id = db.Column(db.Integer, primary_key=True)
    align = db.Column(db.Integer, nullable=False)
    tf_id = db.Column(db.Integer, db.ForeignKey('transcription_factors.tf_id'), nullable=False)
    cl_id = db.Column(db.Integer, db.ForeignKey('cell_lines.cl_id'), nullable=False)
    geo_gse = db.Column(db.String(10))
    encode = db.Column(db.String(11))

    transcription_factor = db.relationship('TranscriptionFactor', backref='experiments')
    cell_line = db.relationship('CellLine', backref='experiments')

    def __repr__(self):
        return '<Experiment #{0.exp_id}>'.format(self)


class ExpSNP(db.Model):
    __tablename__ = 'exp_snps'
    __table_args__ = (
        db.UniqueConstraint('exp_id', 'tf_snp_id',
                            name='unique_tf_aggregated_snp'),
        db.UniqueConstraint('exp_id', 'cl_snp_id',
                            name='unique_cl_aggregated_snp'),
        db.Index('ag_snp_index', 'exp_id', 'tf_snp_id', 'cl_snp_id'),
    )

    exp_snp_id = db.Column(db.Integer, primary_key=True)
    ref_readcount = db.Column(db.Integer, nullable=False)
    alt_readcount = db.Column(db.Integer, nullable=False)
    p_value_ref = db.Column(db.Float)
    p_value_alt = db.Column(db.Float)
    bad = db.Column(db.Enum(*bads))
    tf_snp_id = db.Column(db.Integer, db.ForeignKey('tf_snps.tf_snp_id'))
    cl_snp_id = db.Column(db.Integer, db.ForeignKey('cl_snps.cl_snp_id'))
    exp_id = db.Column(db.Integer, db.ForeignKey('experiments.exp_id'), nullable=False)

    tf_aggregated_snp = db.relationship('TranscriptionFactorSNP', backref='exp_snps')
    cl_aggregated_snp = db.relationship('CellLineSNP', backref='exp_snps')
    experiment = db.relationship('Experiment', backref='exp_snps')

    def __repr__(self):
        return '<ExpSNP #{0.exp_snp_id}, {0.exp_id}, {0.tf_snp_id}, {0.cl_snp_id}>'.format(self)


class GenomePolymorphismLocation(db.Model):
    __abstract__ = True

    chromosome = db.Column(db.Enum(*chromosomes), nullable=False)
    position = db.Column(db.Integer, nullable=False)
    alt = db.Column(db.Enum(*nucleotides), nullable=False)


class SNP(GenomePolymorphismLocation):
    __tablename__ = 'snps'
    __table_args__ = (
        db.PrimaryKeyConstraint('chromosome', 'position', 'alt'),
        db.Index('rs_index', 'rs_id'),
    )

    ref = db.Column(db.Enum(*nucleotides), nullable=False)
    rs_id = db.Column(db.Integer, nullable=False)

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
    motif_log_p_ref = db.Column(db.Float)
    motif_log_p_alt = db.Column(db.Float)
    motif_log_2_fc = db.Column(db.Float)
    motif_orientation = db.Column(db.Boolean)
    motif_position = db.Column(db.Integer)
    motif_concordance = db.Column(db.Boolean)

    @hybrid_property
    def best_p_value(self):
        return min(self.log_p_value_alt, self.log_p_value_ref)

    @best_p_value.expression
    def best_p_value(cls):
        return db.func.max(cls.log_p_value_alt, cls.log_p_value_ref)


class TranscriptionFactorSNP(AggregatedSNP):
    __tablename__ = 'tf_snps'
    __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                              ['snps.chromosome', 'snps.position', 'snps.alt']),
                      db.Index('unique_tf_mutation_index', 'chromosome', 'position', 'alt', 'tf_id'),
                      )

    tf_snp_id = db.Column(db.Integer, primary_key=True)
    tf_id = db.Column(db.Integer, db.ForeignKey('transcription_factors.tf_id'), nullable=False)

    snp = db.relationship('SNP',
                          backref=db.backref('tf_aggregated_snps',
                                             order_by='TranscriptionFactorSNP.best_p_value.desc()'))
    transcription_factor = db.relationship('TranscriptionFactor', backref='tf_aggregated_snps')

    def __repr__(self):
        return '<TranscriptionFactorSNP #{0.tf_snp_id} at {0.chromosome} {0.position} {0.alt}>'.format(self)


class CellLineSNP(AggregatedSNP):
    __tablename__ = 'cl_snps'
    __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                              ['snps.chromosome', 'snps.position', 'snps.alt']),
                      db.UniqueConstraint('chromosome', 'position', 'alt', 'cl_id',
                                          name='cell_line_unique_mutation'),
                      db.Index('unique_cl_mutation_index', 'chromosome', 'position', 'alt', 'cl_id'),
                      )

    cl_snp_id = db.Column(db.Integer, primary_key=True)
    cl_id = db.Column(db.Integer, db.ForeignKey('cell_lines.cl_id'), nullable=False)

    snp = db.relationship('SNP', backref=db.backref('cl_aggregated_snps',
                                                    order_by='CellLineSNP.best_p_value.desc()'))
    cell_line = db.relationship('CellLine', backref='cl_agregated_snps')

    def __repr__(self):
        return '<CellLineSNP #{0.cl_snp_id} at {0.chromosome} {0.position} {0.alt}>'.format(self)


class Phenotype(db.Model):
    __tablename__ = 'phenotypes'
    __tableargs__ = (db.UniqueConstraint('phenotype_name', 'db_name',
                                         name='unique phenotype in db'),
                     )

    phenotype_id = db.Column(db.Integer, primary_key=True)
    db_name = db.Column(db.String(100))
    phenotype_name = db.Column(db.String(100))

    snps = db.relationship(
        'SNP',
        secondary='phenotypes_SNPs',
        backref='phenotypes'
    )


class PhenotypeSNPCorrespondence(GenomePolymorphismLocation):
    __tablename__ = 'phenotypes_SNPs'
    __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                              ['snps.chromosome', 'snps.position', 'snps.alt']),
                      db.UniqueConstraint('chromosome', 'position', 'alt', 'phenotype_id',
                                          name='unique_phenotype_snp_pair'),
                      )

    pair_id = db.Column(db.Integer, primary_key=True)
    phenotype_id = db.Column(db.Integer, db.ForeignKey('phenotypes.phenotype_id'),
                             nullable=False)


db.create_all()
db.session.commit()
