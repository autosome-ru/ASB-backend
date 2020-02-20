from ASB_app import db
from ASB_app.constants import chromosomes, nucleotides, bads


class TranscriptionFactor(db.Model):
    __tablename__ = 'transcription_factors'

    tf_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), nullable=False)

    def __repr__(self):
        return '<TranscriptionFactor #{0.tf_id}, {0.name}>'.format(self)


class CellLine(db.Model):
    __tablename__ = 'cell_lines'

    cl_id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(50), nullable=False)

    def __repr__(self):
        return '<CellLine #{0.cl_id}, {0.name}>'.format(self)


class Experiment(db.Model):
    __tablename__ = 'experiments'

    exp_id = db.Column(db.Integer, primary_key=True)
    align = db.Column(db.Integer, nullable=False)
    tf_id = db.Column(db.Integer, db.ForeignKey('transcription_factors.tf_id'), nullable=False)
    cl_id = db.Column(db.Integer, db.ForeignKey('cell_lines.cl_id'), nullable=False)

    transcription_factor = db.relationship('TranscriptionFactor', backref='experiments')
    cell_line = db.relationship('CellLine', backref='experiments')

    def __repr__(self):
        return '<Experiment #{0.exp_id}>'.format(self)


class ExpSNP(db.Model):
    __tablename__ = 'exp_snps'

    exp_snp_id = db.Column(db.Integer, primary_key=True)
    ref_readcount = db.Column(db.Integer, nullable=False)
    alt_readcount = db.Column(db.Integer, nullable=False)
    p_value_ref = db.Column(db.Float)
    p_value_alt = db.Column(db.Float)
    bad = db.Column(db.Enum(*bads))
    tf_snp_id = db.Column(db.Integer, db.ForeignKey('tf_snps.tf_snp_id'), nullable=False)
    cl_snp_id = db.Column(db.Integer, db.ForeignKey('cl_snps.cl_snp_id'), nullable=False)
    exp_id = db.Column(db.Integer, db.ForeignKey('experiments.exp_id'), nullable=False)

    tf_aggregated_snp = db.relationship('TranscriptionFactorSNP', backref='exp_snps')
    cl_aggregated_snp = db.relationship('CellLineSNP', backref='exp_snps')
    experiment = db.relationship('Experiment', backref='exp_snps')

    def __repr__(self):
        return '<ExpSNP #{0.exp_snp_id}, {0.chromosome}, {0.position}, {0.ref}, {0.alt}>'.format(self)


class GenomePolymorphismLocation(db.Model):
    __abstract__ = True

    chromosome = db.Column(db.Enum(*chromosomes))
    position = db.Column(db.Integer)
    ref = db.Column(db.Enum(*nucleotides))
    alt = db.Column(db.Enum(*nucleotides))


class SNP(GenomePolymorphismLocation):
    __tablename__ = 'snps'
    __table_args__ = (
        db.PrimaryKeyConstraint('chromosome', 'position', 'alt'),
    )

    rs_id = db.Column(db.Integer, nullable=False)

    def __repr__(self):
        return '<SNP {0.rs_id}, {0.chromosome}, {0.position}, {0.ref}, {0.alt}>'.format(self)


class AggregatedSNP(GenomePolymorphismLocation):
    __abstract__ = True

    p_value_ref = db.Column(db.Float)
    p_value_alt = db.Column(db.Float)
    is_asb = db.Column(db.Boolean, nullable=False)
    es_ref = db.Column(db.Float)
    es_alt = db.Column(db.Float)


class TranscriptionFactorSNP(AggregatedSNP):
    __tablename__ = 'tf_snps'
    __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'alt'],
                                              ['snps.chromosome', 'snps.position', 'snps.alt']),
                      db.UniqueConstraint('chromosome', 'position', 'alt', 'tf_id',
                                          name='transcription_factor_unique_mutation'),
                      )

    tf_snp_id = db.Column(db.Integer, primary_key=True)
    tf_id = db.Column(db.Integer, db.ForeignKey('transcription_factors.tf_id'), nullable=False)

    snp = db.relationship('SNP', backref='tf_aggregated_snps')
    transcription_factor = db.relationship('TranscriptionFactor', backref='tf_aggregated_snps')

    def __repr__(self):
        return '<TranscriptionFactorSNP #{0.tf_snp_id}>'.format(self)


class CellLineSNP(AggregatedSNP):
    __tablename__ = 'cl_snps'
    __table_args__ = (db.ForeignKeyConstraint(['chromosome', 'position', 'ref', 'alt'],
                                              ['snps.chromosome', 'snps.position', 'snps.ref', 'snps.alt']),
                      db.UniqueConstraint('chromosome', 'position', 'alt', 'cl_id',
                                          name='cell_line_unique_mutation'),
                      )

    cl_snp_id = db.Column(db.Integer, primary_key=True)
    cl_id = db.Column(db.Integer, db.ForeignKey('cell_lines.cl_id'), nullable=False)

    snp = db.relationship('SNP', backref='cl_aggregated_snps')
    cell_line = db.relationship('CellLine', backref='cl_agregated_snps')

    def __repr__(self):
        return '<CellLineSNP #{0.cl_snp_id}>'.format(self)


db.create_all()
db.session.commit()
