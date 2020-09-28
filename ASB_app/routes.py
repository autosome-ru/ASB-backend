from ASB_app import api, logger, service
from ASB_app.serializers import rs_snp_model_full, transcription_factor_model, cell_line_model, \
    search_results_model, frontpage_statistics_model, gene_model
from ASB_app.constants import chromosomes
from ASB_app.exceptions import ParsingError
from flask_restplus import Resource, inputs
from sqlalchemy.orm.exc import NoResultFound

from ASB_app.service import get_gene_by_id, get_gene_by_name
from ASB_app.utils import PaginationMixin

snp_nsp = api.namespace('SNPs', path='/snps', description='Access to Single Nucleotide Polymorphisms')
search_nsp = api.namespace('Search', path='/search', description='Search SNPs')
browse_nsp = api.namespace('Browse', path='/browse', description='Catalog of TFs and Cell types')

pagination_parser = api.parser()
pagination_parser.add_argument('page', type=inputs.positive, help='Page number', default=1)
pagination_parser.add_argument('size', type=inputs.natural, help='Items per page or 0 for all items', default=0)
pagination_parser.add_argument('offset', type=inputs.natural, help='Skip first N items', default=0)
pagination_parser.add_argument('filter', help='Comma-separated filters: field1 EQ "value1", field2 GE "value2"')
pagination_parser.add_argument('order_by', help='ORDER BY criterion: "field1", "-field2"')


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
class SearchSNPByIdCollection(Resource, PaginationMixin):
    BaseEntity = service.SNP

    @api.marshal_with(search_results_model, PaginationMixin)
    @api.expect(pagination_parser)
    def get(self, rs_id):
        """
        Get all SNPs by rs-ID short info
        """
        all_args = pagination_parser.parse_args()
        filters = service.get_filters_by_rs_id(rs_id)
        result = self.paginate(all_args, extra_filters=filters)
        return {'results': result, 'total': self.items_count(extra_filters=filters)}


@search_nsp.route('/snps/gene_id/<string:gene_id>')
class SearchSNPByGeneIdCollection(Resource, PaginationMixin):
    BaseEntity = service.SNP

    @api.marshal_with(search_results_model)
    @api.expect(pagination_parser)
    def get(self, gene_id):
        """
        Get all SNPs by genome position short info
        """
        all_args = pagination_parser.parse_args()
        gene = get_gene_by_id(gene_id)
        if gene is None:
            return {'results': [], 'gene': None, 'total': 0}

        filters = service.get_filters_by_gene(get_gene_by_id(gene_id))
        result = self.paginate(all_args, extra_filters=filters)

        return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


@search_nsp.route('/snps/gene_name/<string:gene_name>')
class SearchSNPByGeneNameCollection(Resource, PaginationMixin):
    BaseEntity = service.SNP

    @api.marshal_with(search_results_model)
    @api.expect(pagination_parser)
    def get(self, gene_name):
        """
        Get all SNPs by genome position short info
        """
        all_args = pagination_parser.parse_args()
        gene = get_gene_by_name(gene_name)
        if gene is None:
            return {'results': [], 'gene': None, 'total': 0}
        filters = service.get_filters_by_gene(gene)
        result = self.paginate(all_args, extra_filters=filters)

        return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


search_parser = pagination_parser.copy()
search_parser.add_argument('cell_types', action='split', help='Comma-separated list of cell types, search SNPs ASB for every cell type scpecified')
search_parser.add_argument('transcription_factors', action='split', help='Comma-separated list of cell types, search SNPs ASB for every cell type scpecified')
search_parser.add_argument('chromosome', choices=chromosomes, help='Search only SNPs on the specified chromosome')
search_parser.add_argument('start', type=inputs.positive, help='Search SNPs in interval from specified position, Requires "chromosome" and "end"')
search_parser.add_argument('end', type=inputs.positive, help='Search SNPs in interval to specified position, Requiers "chromosome" and "start"')
search_parser.add_argument('phenotype_databases', action='split', help='Comma-separated list of databases, possible choices {grasp, ebi, clinvar, phewas, finemapping, QTL}, earch SNPs that have phenotype associations in all specified databases')
search_parser.add_argument('motif_concordance', action='split', help='Comma-separated list of motif concordance values, possible choices {Concordant. Discordant, Weak Concordant, Weak Discordant}, if no TF specified will search SNPs with any TF having any of the specified concordance valuese, else only SNPs ASB for the specified TFs and any of the specified concordance values for these TFs')


@search_nsp.route('/snps/advanced')
class AdvancedSearchSNP(Resource, PaginationMixin):
    BaseEntity = service.SNP

    @api.marshal_with(search_results_model)
    @api.expect(search_parser)
    def get(self):
        """
        Get all SNPs with advanced filters short info
        """
        all_args = search_parser.parse_args()
        try:
            filters = service.construct_advanced_filters(all_args)
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
            return service.get_snps_by_advanced_filters_tsv(search_parser.parse_args())
        except ParsingError:
            api.abort(400)


@browse_nsp.route('/tf')
class TransctiptionFactorBrowse(Resource, PaginationMixin):
    BaseEntity = service.TranscriptionFactor

    @api.marshal_list_with(transcription_factor_model)
    @api.expect(pagination_parser)
    def get(self):
        """
        Get the list of transcription factors available in the database
        """
        return self.paginate(pagination_parser.parse_args(),
                             extra_filters=(self.BaseEntity.aggregated_snps_count > 0, ))


@browse_nsp.route('/cl')
class CellLineBrowse(Resource, PaginationMixin):
    BaseEntity = service.CellLine

    @api.marshal_list_with(cell_line_model)
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
    @api.marshal_with(frontpage_statistics_model)
    def get(self):
        return service.get_overall_statistics()


used_hints_parser = api.parser()
used_hints_parser.add_argument('options', action='split')
used_hints_parser.add_argument('search')


@search_nsp.route('/tf/hint')
@api.hide
class TransctiptionFactorHint(Resource):
    @api.expect(used_hints_parser)
    @api.marshal_list_with(transcription_factor_model)
    def get(self):
        args = used_hints_parser.parse_args()
        print(args)
        return service.get_hints('TF', args.get('search', ''), args.get('options', []))


@search_nsp.route('/cl/hint')
@api.hide
class CellLineHint(Resource):
    @api.expect(used_hints_parser)
    @api.marshal_list_with(cell_line_model)
    def get(self):
        args = used_hints_parser.parse_args()
        return service.get_hints('CL', args.get('search', ''), args.get('options', []))


csv_columns_parser = api.parser()
csv_columns_parser.add_argument('columns', action='split', required=True)
csv_columns_parser.add_argument('filter')


@search_nsp.route('/gene_name/hint')
@api.hide
class GeneNameHint(Resource):
    @api.expect(used_hints_parser)
    @api.marshal_list_with(gene_model)
    def get(self):
        args = used_hints_parser.parse_args()
        return service.get_hints_for_gene_name(args.get('search', ''))


@snp_nsp.route('/<int:rs_id>/<string:alt>/<string:what_for>/tsv')
class SNPItemCSV(Resource):
    @api.expect(csv_columns_parser)
    def get(self, rs_id, alt, what_for):
        """
        Get complete information about a SNP by rs-ID and Alternative allele in TSV file
        """
        args = csv_columns_parser.parse_args()
        try:
            return service.get_full_snp_tsv(what_for, rs_id, alt, args['columns'])
        except NoResultFound:
            api.abort(404)
