import os
import tempfile
from datetime import datetime

from flask import request

from flask_restplus import Resource
from ASB_app import executor
from ASB_app.constants import fdr_choices, background_choices
from ASB_app.executor_jobs import process_snp_file
from ASB_app.serializers import ticket_model, ticket_model_short
from ASB_app.service import ananastra_service, FileNotProcessed, get_path_by_ticket_id
from ASB_app.service import get_ticket_id_from_path, get_tickets_dir
from ASB_app.utils import PaginationMixin
from ASB_app.models import Ticket

from ASB_app.releases import current_release

from ASB_app.routes import file_parser, pagination_parser, result_param_parser, thresholds_parser

api = current_release.api

ananastra_nsp = api.namespace('ANANASTRA web-service', path='/ananastra', description='SNP annotation by ananastra')


commit_parser = file_parser.copy()
commit_parser.add_argument('user_id')


@ananastra_nsp.route('/commit')
class CommitFile(Resource):
    @api.expect(commit_parser)
    @api.marshal_with(ticket_model)
    def post(self):
        """
        Accepts a file with rs-ids of SNPs to annotate
        """
        if 'file' not in request.files:
            api.abort(400, 'No files')
        else:
            fd, filename = tempfile.mkstemp(prefix=current_release.name + '_', suffix='.tsv',
                                            dir=get_tickets_dir('accepted'))
            request.files['file'].save(filename)
            os.close(fd)
            ticket_id = get_ticket_id_from_path(filename)
            return ananastra_service.create_ticket(ticket_id, commit_parser.parse_args()['user_id'])


@ananastra_nsp.route('/process/<string:ticket_id>')
class ProcessTicket(Resource):
    @api.response(202, 'Ticket accepted for processing')
    @api.response(404, 'Ticket not found')
    @api.expect(thresholds_parser)
    def post(self, ticket_id):
        """
        Submits a ticket for processing
        """
        args = thresholds_parser.parse_args()
        ananastra_service.start_processing_ticket(ticket_id)
        process_snp_file.submit_stored(ticket_id, ticket_id, fdr_class=args['fdr'], background=args['background'])
        return {'message': 'success'}, 202


@ananastra_nsp.route('/ticket/ping/<string:ticket_id>')
class TicketPingItem(Resource):
    @api.marshal_with(ticket_model_short)
    def get(self, ticket_id):
        """
        Get ticket info
        """
        ticket = ananastra_service.get_ticket(ticket_id)
        processing_start = ticket.meta_info.get('processing_started_at')
        if processing_start:
            ticket.elapsed_time = round((datetime.now() - datetime.strptime(processing_start, '%Y-%m-%d %H:%M:%S.%f')).total_seconds())
        else:
            ticket.elapsed_time = None
        ticket.status_details = ticket.meta_info.get('status_details')
        ticket.processing_started_at = ticket.meta_info.get('processing_started_at')
        if not executor.futures.done(ticket_id):
            ticket.position_in_queue = ananastra_service.get_ticket_position_in_queue(ticket_id)
            if ticket.position_in_queue == 0:
                ticket.position_in_queue = None
        else:
            ticket.position_in_queue = None
        return ticket


user_id_parser = api.parser()
user_id_parser.add_argument('user_id')


@ananastra_nsp.route('/ticket/<string:ticket_id>')
class TicketItem(Resource):
    @api.marshal_with(ticket_model)
    def get(self, ticket_id):
        """
        Get ticket info
        """
        return ananastra_service.get_ticket(ticket_id)

    @api.expect(user_id_parser)
    @api.response(403, 'File is processing')
    def delete(self, ticket_id):
        """
        Delete ticket and corresponding files
        """
        if not ananastra_service.get_ticket(ticket_id).user_id == user_id_parser.parse_args()['user_id']:
            return {'message': 'user ids do not match'}, 403
        if not ananastra_service.delete_ticket(ticket_id):
            return {'message': 'file is processing'}, 403
        return {'message': 'success'}, 200


@ananastra_nsp.route('/result/<string:ticket_id>')
class ProcessingResult(Resource):
    @api.expect(result_param_parser)
    def get(self, ticket_id):
        """
        Get first <limit> rows of result
        """
        args = result_param_parser.parse_args()
        result_param = args['result_param']
        format = args['format']
        limit = args['limit']
        try:
            if format == 'json':
                if result_param in ('all', 'not_found'):
                    return {'message': 'JSON format is not available for specified options'}, 403
                return ananastra_service.get_result(ticket_id, result_param, limit, format), 200
            else:
                return ananastra_service.get_result(ticket_id, result_param, limit, format)
        except FileNotProcessed:
            return {'message': 'file is not processed'}, 403


