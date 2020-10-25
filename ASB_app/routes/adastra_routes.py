from ASB_app.service import ReleaseService
from ASB_app.serializers import ReleaseSerializers
from ASB_app.exceptions import ParsingError
from flask_restplus import Resource
from sqlalchemy.orm.exc import NoResultFound

from ASB_app.utils import PaginationMixin
from ASB_app.routes import search_parser, csv_columns_parser, used_hints_parser, pagination_parser

from ASB_app.releases import Release

for release in Release.__subclasses__():
    api = release.api

    snp_nsp = api.namespace('SNPs', path='/snps', description='Access to Single Nucleotide Polymorphisms')
    search_nsp = api.namespace('Search', path='/search', description='Search SNPs')
    browse_nsp = api.namespace('Browse', path='/browse', description='Catalog of TFs and Cell types')

    release_service = ReleaseService(release)
    release_serializers = ReleaseSerializers(release)

    @snp_nsp.route('/<int:rs_id>/<string:alt>')
    class SNPItem(Resource):
        @api.marshal_with(release_serializers.rs_snp_model_full)
        def get(self, rs_id, alt):
            """
            Get complete imformation about an SNP by rs-ID and alt allele
            """
            try:
                return release_service.get_full_snp(rs_id, alt)
            except NoResultFound:
                api.abort(404)


    @search_nsp.route('/snps/rs/<int:rs_id>')
    class SearchSNPByIdCollection(Resource, PaginationMixin):
        BaseEntity = release_service.SNP

        @api.marshal_with(release_serializers.search_results_model, PaginationMixin)
        @api.expect(pagination_parser)
        def get(self, rs_id):
            """
            Get all SNPs by rs-ID short info
            """
            all_args = pagination_parser.parse_args()
            filters = release_service.get_filters_by_rs_id(rs_id)
            result = self.paginate(all_args, extra_filters=filters)
            return {'results': result, 'total': self.items_count(extra_filters=filters)}


    @search_nsp.route('/snps/gene_id/<string:gene_id>')
    class SearchSNPByGeneIdCollection(Resource, PaginationMixin):
        BaseEntity = release_service.SNP

        @api.marshal_with(release_serializers.search_results_model)
        @api.expect(pagination_parser)
        def get(self, gene_id):
            """
            Get all SNPs by genome position short info
            """
            all_args = pagination_parser.parse_args()
            gene = release_service.get_gene_by_id(gene_id)
            if gene is None:
                return {'results': [], 'gene': None, 'total': 0}

            filters = release_service.get_filters_by_gene(release_service.get_gene_by_id(gene_id))
            result = self.paginate(all_args, extra_filters=filters)

            return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


    @search_nsp.route('/snps/gene_name/<string:gene_name>')
    class SearchSNPByGeneNameCollection(Resource, PaginationMixin):
        BaseEntity = release_service.SNP

        @api.marshal_with(release_serializers.search_results_model)
        @api.expect(pagination_parser)
        def get(self, gene_name):
            """
            Get all SNPs by genome position short info
            """
            all_args = pagination_parser.parse_args()
            gene = release_service.get_gene_by_name(gene_name)
            if gene is None:
                return {'results': [], 'gene': None, 'total': 0}
            filters = release_service.get_filters_by_gene(gene)
            result = self.paginate(all_args, extra_filters=filters)

            return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


    @search_nsp.route('/snps/advanced')
    class AdvancedSearchSNP(Resource, PaginationMixin):
        BaseEntity = release_service.SNP

        @api.marshal_with(release_serializers.search_results_model)
        @api.expect(search_parser)
        def get(self):
            """
            Get all SNPs with advanced filters short info
            """
            all_args = search_parser.parse_args()
            try:
                filters = release_service.construct_advanced_filters(all_args)
                result = self.paginate(all_args, extra_filters=filters)
                return {'results': result, 'total': self.items_count(extra_filters=filters)}
            except ParsingError:
                api.abort(400)


    @search_nsp.route('/snps/advanced/tsv')
    class AdvancedSearchSNPCSV(Resource):
        @api.expect(search_parser)
        def get(self):
            """
            Get all SNPs with advanced filters short info in tsv file
            """
            try:
                return release_service.get_snps_by_advanced_filters_tsv(search_parser.parse_args())
            except ParsingError:
                api.abort(400)


    @browse_nsp.route('/tf')
    class TransctiptionFactorBrowse(Resource, PaginationMixin):
        BaseEntity = release_service.TranscriptionFactor

        @api.marshal_list_with(release_serializers.transcription_factor_model)
        @api.expect(pagination_parser)
        def get(self):
            """
            Get the list of transcription factors available in the database
            """
            return self.paginate(pagination_parser.parse_args(),
                                 extra_filters=(self.BaseEntity.aggregated_snps_count > 0, ))


    @browse_nsp.route('/cl')
    class CellLineBrowse(Resource, PaginationMixin):
        BaseEntity = release_service.CellLine

        @api.marshal_list_with(release_serializers.cell_line_model)
        @api.expect(pagination_parser)
        def get(self):
            """
            Get the list of cell types available in the database
            """
            return self.paginate(pagination_parser.parse_args(),
                                 extra_filters=(self.BaseEntity.aggregated_snps_count > 0, ))


    @browse_nsp.route('/total')
    @api.hide
    class FrontPageStatistics(Resource):
        @api.marshal_with(release_serializers.frontpage_statistics_model)
        def get(self):
            return release_service.get_overall_statistics()


    @search_nsp.route('/tf/hint')
    @api.hide
    class TranscriptionFactorHint(Resource):
        @api.expect(used_hints_parser)
        @api.marshal_list_with(release_serializers.transcription_factor_model)
        def get(self):
            args = used_hints_parser.parse_args()
            return release_service.get_hints('TF', args.get('search', ''), args.get('options', []))


    @search_nsp.route('/cl/hint')
    @api.hide
    class CellLineHint(Resource):
        @api.expect(used_hints_parser)
        @api.marshal_list_with(release_serializers.cell_line_model)
        def get(self):
            args = used_hints_parser.parse_args()
            return release_service.get_hints('CL', args.get('search', ''), args.get('options', []))


    @search_nsp.route('/gene_name/hint')
    @api.hide
    class GeneNameHint(Resource):
        @api.expect(used_hints_parser)
        @api.marshal_list_with(release_serializers.gene_model)
        def get(self):
            args = used_hints_parser.parse_args()
            return release_service.get_hints_for_gene_name(args.get('search', ''))


    @snp_nsp.route('/<int:rs_id>/<string:alt>/<string:what_for>/tsv')
    class SNPItemCSV(Resource):
        @api.expect(csv_columns_parser)
        def get(self, rs_id, alt, what_for):
            """
            Get complete information about a SNP by rs-ID and Alternative allele in TSV file
            """
            args = csv_columns_parser.parse_args()
            try:
                return release_service.get_full_snp_tsv(what_for, rs_id, alt, args['columns'])
            except NoResultFound:
                api.abort(404)