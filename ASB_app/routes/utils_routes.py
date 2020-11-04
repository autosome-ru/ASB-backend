import os
import tempfile

from flask import request

from flask_restplus import Resource
from ASB_app.service.utils_service import UtilsService
from ASB_app.routes import file_parser
from ASB_app.releases import Release

context_parser = file_parser.copy()
context_parser.add_argument('field', default='RS-ID', help='Name of column corresponding to SNP rs id in tsv file')


def set_utils_service(utils_service):
    def wrapper(cls):
        cls.utils_service = utils_service
        return cls
    return wrapper


for release in Release.__subclasses__():
    api = release.api
    utils_nsp = api.namespace('ADASTRA utils', path='/utils', description='SNP annotation utils')

    utils_service = UtilsService(release)

    @utils_nsp.route('/context')
    @set_utils_service(utils_service)
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
                ok, result = self.utils_service.annotate_with_context(filename, context_parser.parse_args()['field'])
                if not ok:
                    return 400, 'Bad file'
                return result
