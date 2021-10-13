from flask_restplus import fields

from ASB_app.constants import background_choices
from ASB_app.releases import current_release

api = current_release.api

asb_data_model = api.model('ASB data entity', {
    'name': fields.String,
    'asbs': fields.Integer,
    'negatives': fields.Integer,
    'asbs_rs': fields.Integer,
    'negatives_rs': fields.Integer,
    'expected_asbs_rs': fields.Integer,
    'expected_negatives_rs': fields.Integer,
    'odds': fields.String,
    'log10_p_value': fields.String,
    'log10_fdr': fields.String,
})


chr_asb_data_model = api.clone('ASB data by chromosome entity', asb_data_model, {
    'tf_asbs': fields.Integer,
    'tf_negatives': fields.Integer,
    'tf_asbs_rs': fields.Integer,
    'tf_negatives_rs': fields.Integer,
    'expected_tf_asbs_rs': fields.Integer,
    'expected_tf_negatives_rs': fields.Integer,
    'tf_odds': fields.String,
    'tf_log10_p_value': fields.String,

    'cl_asbs': fields.Integer,
    'cl_negatives': fields.Integer,
    'cl_asbs_rs': fields.Integer,
    'cl_negatives_rs': fields.Integer,
    'expected_cl_asbs_rs': fields.Integer,
    'expected_cl_negatives_rs': fields.Integer,
    'cl_odds': fields.String,
    'cl_log10_p_value': fields.String,
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


meta_info_for_snps = api.model('Ticket meta info stats', {
    'asbs': fields.Integer,
    'asbs_rs': fields.Integer,
    'negatives': fields.Integer,
    'negatives_rs': fields.Integer,
    'expected_asbs_rs': fields.Integer,
    'expected_negatives_rs': fields.Integer,
    'odds_rs': fields.String,
    'log10_p_value_rs': fields.String,
})


meta_info_for_aggregation_level = api.clone('Ticket meta info for aggregation level', meta_info_for_snps, {
    'asb_counts': fields.List(fields.Nested(asb_count_model)),
    'asb_counts_top': fields.List(fields.Nested(asb_count_model)),
    'asb_data': fields.List(fields.Nested(asb_data_model)),
})


meta_info_for_chromosomes = api.model('Ticket meta info by chromosome', {
    'log10_p_value_rs': fields.String,
    'tf_log10_p_value_rs': fields.String,
    'cl_log10_p_value_rs': fields.String,
    'asb_data': fields.List(fields.Nested(chr_asb_data_model)),
})


meta_info_model = api.model('Ticket meta info', {
    'processing_time': fields.String,
    'all_rs': fields.Integer,
    'submitted_snps_count': fields.Integer,
    'unique_submitted_snps_count': fields.Integer,
    'undefined_rs': fields.Integer,
    'tf': fields.Nested(meta_info_for_aggregation_level),
    'cl': fields.Nested(meta_info_for_aggregation_level),
    'all': fields.Nested(meta_info_for_snps),
    'chr': fields.Nested(meta_info_for_chromosomes),
})

ticket_model = api.model('ANANASTRA ticket', {
    'ticket_id': fields.String,
    'date_created': fields.DateTime,
    'expiration_date': fields.DateTime,
    'status': fields.String(enum=('Created', 'Processing', 'Processed', 'Failed')),
    'fdr': fields.String,
    'es': fields.String,
    'background': fields.String(enum=background_choices),
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
    'es': fields.String,
    'background': fields.String(enum=background_choices),
})
