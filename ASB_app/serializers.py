from flask_restplus import fields
from ASB_app import api
from ASB_app.constants import chromosomes, bads, nucleotides

genome_polymorphism_location_model = api.model('Genome Mutation Position', {
    'chromosome': fields.String(enumerate=chromosomes),
    'position': fields.Integer,
    'ref': fields.String(enumerate=nucleotides),
    'alt': fields.String(enumerate=nucleotides),
})

rs_snp_model_short = api.clone('Short rs-SNP info', genome_polymorphism_location_model, {
    'rs_id': fields.Integer,
})

aggregated_snp_model_short = api.model('Agregated SNP (no genome info) ', {
    'p_value_ref': fields.Float,
    'p_value_alt': fields.Float,
    'is_asb': fields.Boolean,
    'es_ref': fields.Float,
    'es_alt': fields.Float,
})

transcription_factor_model = api.model('Transcription factor', {
    'tf_id': fields.Integer(readonly=True),
    'name': fields.String,
})

cell_line_model = api.model('Cell line', {
    'cl_id': fields.Integer(readonly=True),
    'name': fields.String,
})

tf_snp_model_short = api.inherit('Transcription Factor SNP (no genome info)', aggregated_snp_model_short, {
    'tf_snp_id': fields.Integer(readonly=True),
    'transcription_factor': fields.Nested(transcription_factor_model),
})

cl_snp_model_short = api.inherit('Cell line SNP (no genome info)', aggregated_snp_model_short, {
    'cl_snp_id': fields.Integer(readonly=True),
    'cell_line': fields.Nested(cell_line_model),
})

rs_snp_model = api.clone('Complete rs-SNP info', genome_polymorphism_location_model, {
    'rs_id': fields.Integer,
    'tf_aggregated_snps': fields.List(fields.Nested(tf_snp_model_short)),
    'cl_aggregated_snps': fields.List(fields.Nested(cl_snp_model_short)),
})
