from flask_restplus import fields
from ASB_app import api
from ASB_app.constants import chromosomes, bads, nucleotides

genome_polymorphism_location_model = api.model('Genome Mutation Position and id', {
    'chromosome': fields.String(enumerate=chromosomes),
    'position': fields.Integer,
    'ref': fields.String(enumerate=nucleotides),
    'alt': fields.String(enumerate=nucleotides),
    'rs_id': fields.Integer,
})

aggregated_snp_model = api.model('Agregated SNP (no genome info) ', {
    'p_value_ref': fields.Float,
    'p_value_alt': fields.Float,
    'is_asb': fields.Boolean,
})

transcription_factor_model = api.model('Transcription factor', {
    'tf_id': fields.Integer(readonly=True),
    'name': fields.String,
})

cell_line_model = api.model('Cell line', {
    'cl_id': fields.Integer(readonly=True),
    'name': fields.String,
})

tf_snp_model = api.inherit('Transcription Factor SNP (no genome info)', aggregated_snp_model, {
    'tf_snp_id': fields.Integer(readonly=True),
    'transcription_factor': fields.Nested(transcription_factor_model),
})

cl_snp_model = api.inherit('Cell line SNP (no genome info)', aggregated_snp_model, {
    'cl_snp_id': fields.Integer(readonly=True),
    'cell_line': fields.Nested(cell_line_model),
})

rs_snp_model = api.inherit('rs-SNP info for search', genome_polymorphism_location_model, {
    'tf_aggregated_snps': fields.List(fields.Nested(tf_snp_model)),
    'cl_aggregated_snps': fields.List(fields.Nested(cl_snp_model)),
})

exp_snp_model = api.model('Experiment SNP', {
    'exp_snp_id': fields.Integer(readonly=True),
    'ref_readcount': fields.Integer,
    'alt_readcount': fields.Integer,
    'p_value_ref': fields.Float,
    'p_value_alt': fields.Float,
    'bad': fields.String(enumerate=bads)
})

aggregated_snp_model_full = api.inherit('Agregated SNP (with exp snps)', aggregated_snp_model, {
    'p_value_ref': fields.Float,
    'p_value_alt': fields.Float,
    'is_asb': fields.Boolean,
    'es_ref': fields.Float,
    'es_alt': fields.Float,
    'exp_snps': fields.List(fields.Nested(exp_snp_model))
})

tf_snp_model_full = api.inherit('Transcription Factor SNP (with exp snps)', aggregated_snp_model_full, {
    'tf_snp_id': fields.Integer(readonly=True),
    'transcription_factor': fields.Nested(transcription_factor_model),
})

cl_snp_model_full = api.inherit('Cell line SNP (with exp snps)', aggregated_snp_model_full, {
    'cl_snp_id': fields.Integer(readonly=True),
    'cell_line': fields.Nested(cell_line_model),
})

rs_snp_model_full = api.inherit('Complete rs-SNP info (with exp snps)', genome_polymorphism_location_model, {
    'tf_aggregated_snps': fields.List(fields.Nested(tf_snp_model_full)),
    'cl_aggregated_snps': fields.List(fields.Nested(cl_snp_model_full)),
})
