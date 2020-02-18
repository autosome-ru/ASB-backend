from ASB_app import api, logger, service
from ASB_app.serializers import tf_snp_model, cl_snp_model, rs_snp_model
from flask import request, jsonify, g
from flask_restplus import Resource, inputs

snp_nsp = api.namespace('SNPs', path='/snps', description='Access to Single Nucleotide Polymorphisms')


@snp_nsp.route('/<int:rs_id>')
class SNPItem(Resource):
    @api.marshal_list_with(rs_snp_model)
    def get(self, rs_id):
        """
        Get complete imformation about an SNP by rs-ID
        """
        return service.get_full_snp(rs_id)
