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
    'log_p_value_ref': fields.Float,
    'log_p_value_alt': fields.Float,
    'is_asb': fields.Boolean,
})

transcription_factor_model = api.model('Transcription factor', {
    'tf_id': fields.Integer(readonly=True),
    'name': fields.String,
    'aggregated_snps_count': fields.Integer(readonly=True),
})

cell_line_model = api.model('Cell line', {
    'cl_id': fields.Integer(readonly=True),
    'name': fields.String,
    'aggregated_snps_count': fields.Integer(readonly=True),
})

tf_snp_model = api.inherit('Transcription Factor SNP (no genome info)', aggregated_snp_model, {
    'tf_snp_id': fields.Integer(readonly=True),
    'motif_log_p_ref': fields.Float,
    'motif_log_p_alt': fields.Float,
    'motif_log_2_fc': fields.Float,
    'motif_orientation': fields.Boolean,
    'motif_position': fields.Integer,
    'motif_concordance': fields.Boolean,
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

search_results_model = api.model('SNP search results model', {
    'total': fields.Integer,
    'results': fields.List(fields.Nested(rs_snp_model))
})

exp_model_short = api.model('Experiment (short info)', {
    'exp_id': fields.Integer,
    'align': fields.Integer,
    'tf_name': fields.String(attribute='transcription_factor.name'),
    'cl_name': fields.String(attribute='cell_line.name'),
})

exp_snp_model = api.model('Experiment SNP', {
    'exp_snp_id': fields.Integer(readonly=True),
    'ref_readcount': fields.Integer,
    'alt_readcount': fields.Integer,
    'p_value_ref': fields.Float,
    'p_value_alt': fields.Float,
    'bad': fields.String(enumerate=bads),
    'experiment': fields.Nested(exp_model_short)
})

aggregated_snp_model_full = api.inherit('Agregated SNP (with exp snps)', aggregated_snp_model, {
    'es_ref': fields.Float,
    'es_alt': fields.Float,
    'mean_bad': fields.Float,
    'exp_snps': fields.List(fields.Nested(exp_snp_model))
})

tf_snp_model_full = api.inherit('Transcription Factor SNP (with exp snps)', aggregated_snp_model_full, {
    'tf_snp_id': fields.Integer(readonly=True),
    'motif_log_p_ref': fields.Float,
    'motif_log_p_alt': fields.Float,
    'motif_log_2_fc': fields.Float,
    'motif_orientation': fields.Boolean,
    'motif_position': fields.Integer,
    'motif_concordance': fields.Boolean,
    'transcription_factor': fields.Nested(transcription_factor_model),
})

cl_snp_model_full = api.inherit('Cell line SNP (with exp snps)', aggregated_snp_model_full, {
    'cl_snp_id': fields.Integer(readonly=True),
    'cell_line': fields.Nested(cell_line_model),
})

phenotype_model = api.model('Phenotype', {
    'db_name': fields.String,
    'phenotype_name': fields.String,
})

rs_snp_model_full = api.inherit('Complete rs-SNP info (with exp snps)', genome_polymorphism_location_model, {
    'tf_aggregated_snps': fields.List(fields.Nested(tf_snp_model_full)),
    'cl_aggregated_snps': fields.List(fields.Nested(cl_snp_model_full)),
    'phenotypes': fields.List(fields.Nested(phenotype_model))
})
