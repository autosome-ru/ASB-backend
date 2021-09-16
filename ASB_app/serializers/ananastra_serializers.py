from flask_restplus import fields
from ASB_app.releases import current_release

api = current_release.api

asb_data_model = api.model('ASB data entity', {
    'name': fields.String,
    'asbs': fields.Integer,
    'candidates': fields.Integer,
    'asbs_rs': fields.Integer,
    'candidates_rs': fields.Integer,
    'odds': fields.String,
    'log10_p_value': fields.String,
    'log10_fdr': fields.String,
})

asb_count_model = api.model('ASB count entity', {
    'name': fields.String,
    'count': fields.Integer,
})

concordant_asb_model = api.model('Motif-concordant ASB', {
    'tf_name': fields.String,
    'rs_id': fields.String,
    'alt': fields.String(enum=('A', 'C', 'G', 'T')),
    'concordance': fields.String(enum=('Concordant', 'Discordant', 'Weak Concordant', 'Weak Discordant'))
})

meta_info_model = api.model('Ticket meta info', {
    'processing_time': fields.String,
    'processing_started_at': fields.String,
    'last_status_update_at': fields.String,
    'status_details': fields.String,
    'all_rs': fields.Integer,
    'tf_asbs': fields.Integer,
    'cl_asbs': fields.Integer,
    'all_asbs': fields.Integer,
    'tf_candidates': fields.Integer,
    'cl_candidates': fields.Integer,
    'all_candidates': fields.Integer,
    'tf_asbs_rs': fields.Integer,
    'cl_asbs_rs': fields.Integer,
    'all_asbs_rs': fields.Integer,
    'undefined_rs': fields.Integer,
    'tf_candidates_rs': fields.Integer,
    'cl_candidates_rs': fields.Integer,
    'all_candidates_rs': fields.Integer,
    'tf_odds': fields.String,
    'tf_log10_p_value': fields.String,
    'cl_odds': fields.String,
    'cl_log10_p_value': fields.String,
    'all_odds': fields.String,
    'all_log10_p_value': fields.String,
    'tf_odds_rs': fields.String,
    'tf_log10_p_value_rs': fields.String,
    'cl_odds_rs': fields.String,
    'cl_log10_p_value_rs': fields.String,
    'all_odds_rs': fields.String,
    'all_log10_p_value_rs': fields.String,
    'expected_fraction_all': fields.Float,
    'expected_fraction_tf': fields.Float,
    'expected_fraction_cl': fields.Float,
    'tf_asb_counts': fields.List(fields.Nested(asb_count_model)),
    'cl_asb_counts': fields.List(fields.Nested(asb_count_model)),
    'tf_asb_counts_top': fields.List(fields.Nested(asb_count_model)),
    'cl_asb_counts_top': fields.List(fields.Nested(asb_count_model)),
    'tf_asb_data': fields.List(fields.Nested(asb_data_model)),
    'cl_asb_data': fields.List(fields.Nested(asb_data_model)),
    'chr_asb_data': fields.List(fields.Nested(asb_data_model)),
    'chr_log10_p_value_rs': fields.String,
    'concordant_asbs': fields.List(fields.Nested(concordant_asb_model)),
})

ticket_model = api.model('ANANASTRA ticket', {
    'ticket_id': fields.String,
    'date_created': fields.DateTime,
    'expiration_date': fields.DateTime,
    'status': fields.String(enum=('Created', 'Processing', 'Processed', 'Failed')),
    'fdr': fields.String,
    'es': fields.String,
    'meta_info': fields.Nested(meta_info_model),
})

ticket_model_short = api.model('ANANASTRA ticket (short)', {
    'ticket_id': fields.String,
    'date_created': fields.DateTime,
    'expiration_date': fields.DateTime,
    'status': fields.String(enum=('Created', 'Processing', 'Processed', 'Failed')),
    'elapsed_time': fields.Integer,
    'status_details': fields.String,
    'processing_started_at': fields.String,
    'position_in_queue': fields.Integer,
    'fdr': fields.String,
    'es': fields.String
})
