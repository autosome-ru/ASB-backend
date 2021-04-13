from datetime import datetime

from ASB_app.constants import nucleotides, chromosomes, fdr_classes, es_classes
from ASB_app.releases import current_release

db = current_release.db


class Ticket(db.Model):
    __tablename__ = 'tickets'
    __bind_key__ = 'tickets'

    ticket_id = db.Column(db.String(50), primary_key=True)
    status = db.Column(db.Enum('Created', 'Processing', 'Processed', 'Failed'), server_default='Created', nullable=False)
    date_created = db.Column(db.DATETIME, nullable=False)
    expiration_date = db.Column(db.DATETIME)
    meta_info = db.Column(db.JSON, default={})
    user_id = db.Column(db.String(36))
    fdr = db.Column(db.Enum(*fdr_classes))
    es = db.Column(db.Enum(*es_classes))

    def __repr__(self):
        return '<AnanastraTicket {0.ticket_id}, created {0.date_created}, {0.status}>'.format(self)


class GenomePolymorphismLocation(db.Model):
    __abstract__ = True

    chromosome = db.Column(db.Enum(*chromosomes), nullable=False)
    position = db.Column(db.Integer, nullable=False)
    alt = db.Column(db.Enum(*nucleotides), nullable=False)


class CandidateSNP(GenomePolymorphismLocation):
    __tablename__ = 'candidate_snps'
    __bind_key__ = 'candidates'
    __table_args__ = (
        db.PrimaryKeyConstraint('chromosome', 'position', 'alt', 'ag_level', 'ag_id'),
        db.Index('rs_index', 'rs_id'),
        db.Index('ag_level_index', 'ag_level'),
    )

    ref = db.Column(db.Enum(*nucleotides), nullable=False)
    rs_id = db.Column(db.Integer, nullable=False)
    ag_level = db.Column(db.Enum('TF', 'CL'), nullable=False)
    ag_id = db.Column(db.Integer)
    best_p_value = db.Column(db.Float, index=True)
    fdr_class = db.Column(db.Enum(*fdr_classes), index=True)
    best_es = db.Column(db.Float, index=True)
    es_class = db.Column(db.Enum(*es_classes), index=True)

    def __repr__(self):
        return '<CandidateSNP rs{0.rs_id}, {0.alt}, {0.ag_level}, {0.ag_id}, {0.fdr_class}, {0.es_class}>'.format(self)


class CandidateRS(db.Model):
    __tablename__ = 'candidate_rs_snps'
    __bind_key__ = 'candidates'
    __table_args__ = (
        db.PrimaryKeyConstraint('rs_id'),
    )

    rs_id = db.Column(db.Integer, nullable=False)
    best_p_value = db.Column(db.Float, index=True)
    fdr_class = db.Column(db.Enum(*fdr_classes), index=True)
    best_es = db.Column(db.Float, index=True)
    es_class = db.Column(db.Enum(*es_classes), index=True)

    def __repr__(self):
        return '<CandidateRS rs{0.rs_id}, {0.fdr_class}, {0.es_class}>'.format(self)


class CandidateTFRS(db.Model):
    __tablename__ = 'candidate_tf_rs_snps'
    __bind_key__ = 'candidates'
    __table_args__ = (
        db.PrimaryKeyConstraint('rs_id'),
    )

    rs_id = db.Column(db.Integer, nullable=False)
    best_p_value = db.Column(db.Float, index=True)
    fdr_class = db.Column(db.Enum(*fdr_classes), index=True)
    best_es = db.Column(db.Float, index=True)
    es_class = db.Column(db.Enum(*es_classes), index=True)

    def __repr__(self):
        return '<CandidateTFRS rs{0.rs_id}, {0.fdr_class}, {0.es_class}>'.format(self)


class CandidateCLRS(db.Model):
    __tablename__ = 'candidate_cl_rs_snps'
    __bind_key__ = 'candidates'
    __table_args__ = (
        db.PrimaryKeyConstraint('rs_id'),
    )

    rs_id = db.Column(db.Integer, nullable=False)
    best_p_value = db.Column(db.Float, index=True)
    fdr_class = db.Column(db.Enum(*fdr_classes), index=True)
    best_es = db.Column(db.Float, index=True)
    es_class = db.Column(db.Enum(*es_classes), index=True)

    def __repr__(self):
        return '<CandidateCLRS rs{0.rs_id}, {0.fdr_class}, {0.es_class}>'.format(self)
