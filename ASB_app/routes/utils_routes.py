import os
import tempfile

from flask import request

from flask_restplus import Resource
from ASB_app.service import utils_service


from ASB_app.releases import ReleaseTest

release = ReleaseTest

from ASB_app.routes import file_parser

api = release.api

utils_nsp = api.namespace('ADASTRA utils', path='/utils', description='SNP annotation utils')

context_parser = file_parser.copy()
context_parser.add_argument('field', default='RS-ID', help='Name of column corresponding to SNP rs id in tsv file')


@utils_nsp.route('/context')
class ContextFile(Resource):
    @api.expect(context_parser)
    def post(self):
        """
        Accepts a tsv file with rs-ids of SNPs to annotate with context
        """
        if 'file' not in request.files:
            api.abort(400, 'No files')
        else:
            fd, filename = tempfile.mkstemp(suffix='.tsv')
            request.files['file'].save(filename)
            os.close(fd)
            return utils_service.annotate_with_context(release, filename, context_parser.parse_args()['field'])
