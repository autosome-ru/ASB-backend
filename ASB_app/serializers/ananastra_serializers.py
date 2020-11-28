from flask_restplus import fields
from ASB_app.releases import current_release

api = current_release.api

asb_data_model = api.model('ASB data entity', {
    'name': fields.String,
    'asbs': fields.Integer,
    'candidates': fields.Integer,
    'odds': fields.Float,
    'log10_p_value': fields.Float,
    'log10_fdr': fields.Float,
})

asb_count_model = api.model('ASB count entity', {
    'name': fields.String,
    'count': fields.Integer,
    'odds': fields.Float,
})

concordant_asb_model = api.model('Motif-concordant ASB', {
    'tf_name': fields.String,
    'rs_id': fields.String,
    'alt': fields.String(enum=('A', 'C', 'G', 'T')),
    'concordance': fields.String(enum=('Concordant', 'Discordant', 'Weak Concordant', 'Weak Discordant'))
})

meta_info_model = api.model('Ticket meta info', {
    'processing_time': fields.String,
    'all_rs': fields.Integer,
    'tf_asbs': fields.Integer,
    'cl_asbs': fields.Integer,
    'all_asbs': fields.Integer,
    'tf_candidates': fields.Integer,
    'cl_candidates': fields.Integer,
    'all_candidates': fields.Integer,
    'tf_odds': fields.Float,
    'tf_log10_p_value': fields.Float,
    'cl_odds': fields.Float,
    'cl_log10_p_value': fields.Float,
    'all_odds': fields.Float,
    'all_log10_p_value': fields.Float,
    'tf_asb_counts': fields.List(fields.Nested(asb_count_model)),
    'cl_asb_counts': fields.List(fields.Nested(asb_count_model)),
    'tf_asb_data': fields.List(fields.Nested(asb_data_model)),
    'cl_asb_data': fields.List(fields.Nested(asb_data_model)),
    'concordant_asbs': fields.List(fields.Nested(concordant_asb_model)),
})

ticket_model = api.model('ANANASTRA ticket', {
    'ticket_id': fields.String,
    'date_created': fields.DateTime,
    'expiration_date': fields.DateTime,
    'status': fields.String(enum=('Created', 'Processing', 'Processed', 'Failed')),
    'meta_info': fields.Nested(meta_info_model),
})
