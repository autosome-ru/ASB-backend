from flask_restplus import fields
from ASB_app.releases import current_release

api = current_release.api

meta_info_model = api.model('Ticket meta info', {
    'processing_time': fields.String,
    'tf_asbs': fields.Integer,
    'cl_asbs': fields.Integer,
    'all_asbs': fields.Integer,
    'tf_candidates': fields.Integer,
    'cl_canidates': fields.Integer,
    'all_candidates': fields.Integer,
    'tf_odds': fields.Float,
    'tf_log10_p_value': fields.Float,
    'cl_odds': fields.Float,
    'cl_log10_p_value': fields.Float,
    'all_odds': fields.Float,
    'all_log10_p_value': fields.Float,
})

ticket_model = api.model('ANANASTRA ticket', {
    'ticket_id': fields.String,
    'date_created': fields.DateTime,
    'status': fields.String(enum=('Created', 'Processing', 'Processed', 'Failed')),
    'meta_info': fields.Nested(meta_info_model),
})
