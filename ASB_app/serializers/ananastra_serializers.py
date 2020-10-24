from flask_restplus import fields
from ASB_app.releases import current_release

api = current_release.api

ticket_model = api.model('ANANASTRA ticket', {
    'ticket_id': fields.String,
    'date_created': fields.DateTime,
    'status': fields.String(enum=('Created', 'Processing', 'Processed', 'Failed')),
})
