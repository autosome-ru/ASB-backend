from ASB_app import api, logger, service
from ASB_app.serializers import rs_snp_model, rs_snp_model_full
from flask import request, jsonify, g
from flask_restplus import Resource, inputs
from sqlalchemy.orm.exc import NoResultFound

snp_nsp = api.namespace('SNPs', path='/snps', description='Access to Single Nucleotide Polymorphisms')
search_nsp = api.namespace('Search', path='/search', description='Search SNPs')


@snp_nsp.route('/<int:rs_id>/<string:alt>')
class SNPItem(Resource):
    @api.marshal_with(rs_snp_model_full)
    def get(self, rs_id, alt):
        """
        Get complete imformation about an SNP by rs-ID and alt allele
        """
        try:
            return service.get_full_snp(rs_id, alt)
        except NoResultFound:
            api.abort(404)


@search_nsp.route('/snps/rs/<int:rs_id>')
class SNPSearchSNPByIdCollection(Resource):
    @api.marshal_list_with(rs_snp_model)
    def get(self, rs_id):
        """
        Get all SNPs by rs-ID short info
        """
        return service.get_snps_by_rs_id(rs_id)
