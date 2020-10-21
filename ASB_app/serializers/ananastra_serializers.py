from flask_restplus import fields
from ASB_app import api

ticket_model = api.model('ANANASTRA ticket', {
    'ticket_id': fields.String,
    'status': fields.String(enum=('Created', 'Processing', 'Processed', 'Failed')),
})
