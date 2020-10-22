import os
import tempfile

from flask import request

from ASB_app import api
from flask_restplus import Resource

from ASB_app.executor_jobs import process_snp_file
from ASB_app.serializers import ticket_model
from ASB_app.service import ananastra_service
from ASB_app.service import get_ticket_id_from_path, get_tickets_dir
from ASB_app.utils import PaginationMixin
from ASB_app.models import Ticket

from ASB_app.routes import file_parser, pagination_parser

ananastra_nsp = api.namespace('ANANASTRA web-service', path='/ananastra', description='SNP annotation by ananastra')


@ananastra_nsp.route('/commit')
class CommitFile(Resource):
    @api.expect(file_parser)
    @api.marshal_with(ticket_model)
    def post(self):
        """
        Accepts a file with rs-ids of SNPs to annotate
        """
        if 'file' not in request.files:
            api.abort(400, 'No files')
        else:
            fd, filename = tempfile.mkstemp(suffix='.tsv', dir=get_tickets_dir('accepted'))
            request.files['file'].save(filename)
            os.close(fd)
            ticket_id = get_ticket_id_from_path(filename)
            return ananastra_service.create_ticket(ticket_id)


@ananastra_nsp.route('/process/<string:ticket_id>')
class ProcessTicket(Resource):
    @api.response(202, 'Ticket accepted for processing')
    @api.response(404, 'Ticket not found')
    def post(self, ticket_id):
        """
        Submits a ticket for processing
        """
        ananastra_service.update_ticket_status(ticket_id, 'Processing')
        process_snp_file.submit(ticket_id)

        return {'message': 'success'}, 202


@ananastra_nsp.route('/ticket/<string:ticket_id>')
class TicketItem(Resource):
    @api.marshal_with(ticket_model)
    def get(self, ticket_id):
        """
        Get ticket info
        """
        return ananastra_service.get_ticket(ticket_id)

    @api.response(403, 'File is processing')
    def delete(self, ticket_id):
        """
        Delete ticket and corresponding files
        """
        if not ananastra_service.delete_ticket(ticket_id):
            return {'message': 'file is processing'}, 403
        return {'message': 'success'}, 200


@ananastra_nsp.route('/result/<string:ticket_id>/<string:result_param>')
class ProcessingResult(Resource):
    def get(self, ticket_id, result_param):
        """
        Get first 1000 rows of result
        """
        if result_param not in ('tf', 'cl'):
            api.abort(400, 'Wrong result param')
        ok, result = ananastra_service.get_result(ticket_id, result_param)
        if not ok:
            return {'message': 'file is not processed'}, 403
        return result, 200


@ananastra_nsp.route('/ticket')
class TicketCollection(Resource, PaginationMixin):
    BaseEntity = Ticket

    @api.marshal_list_with(ticket_model)
    @api.expect(pagination_parser)
    def get(self):
        """
        Get all tickets
        """
        return self.paginate(pagination_parser.parse_args())
