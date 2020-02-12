from ASB_app import api, logger, service
from ASB_app.serializers import tf_snp_model, cl_snp_model, rs_snp_model
from flask import request, jsonify, g
from flask_restplus import Resource, inputs

snp_nsp = api.namespace('SNPs', path='/snps', description='Access to Single Nucleotide Polymorphisms')


@snp_nsp.route('/TF/<int:rs_id>')
class TranscriptionFactorSNPCollection(Resource):
    @api.marshal_list_with(tf_snp_model)
    def get(self, rs_id):
        """
        Get list of all TF-aggregated SNPs by rs-ID
        """
        return service.get_snps_by_rs(rs_id, 'TF')


@snp_nsp.route('/CL/<int:rs_id>')
class CellLineSNPCollection(Resource):
    @api.marshal_list_with(cl_snp_model)
    def get(self, rs_id):
        """
        Get list of all CL-aggregated SNPs by rs-ID
        """
        return service.get_snps_by_rs(rs_id, 'CL')


@snp_nsp.route('/<int:rs_id>')
class RsSNPItem(Resource):
    @api.marshal_list_with(rs_snp_model)
    def get(self, rs_id):
        """
        Get complete imformation about an SNP by rs-ID
        """
        return service.get_full_snps_by_rs(rs_id)
