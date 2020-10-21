import os
import tempfile

from flask import request
from werkzeug.datastructures import FileStorage

from ASB_app import api
from flask_restplus import Resource

from ASB_app.executor_jobs import process_snp_file
from ASB_app.serializers import ticket_model
from ASB_app.service import ananastra_service
from ASB_app.service import get_ticket_id_from_path

ananastra_nsp = api.namespace('ANANASTRA web-service', path='/ananastra', description='SNP annotation by ananastra')

file_parser = api.parser()
file_parser.add_argument('file', type=FileStorage, location='files')


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
            fd, filename = tempfile.mkstemp(suffix='.tsv', dir=os.path.expanduser('~/adastra/tickets/accepted'))
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


@ananastra_nsp.route('/info/<string:ticket_id>')
class TicketInfo(Resource):
    @api.marshal_with(ticket_model)
    def get(self, ticket_id):
        """
        Get ticket info
        """
        return ananastra_service.get_ticket(ticket_id)


# @ananastra_nsp.route('/result/tf/<string:ticket_id>')
# class ProcessingResult(Resource):
#     @api.marshal_list_with()
