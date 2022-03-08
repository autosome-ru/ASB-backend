from flask import render_template
from sqlalchemy import desc

from ASB_app.constants import default_fdr_tr
from ASB_app.service import ReleaseService
from ASB_app.serializers import ReleaseSerializers
from ASB_app.exceptions import ParsingError, ReleaseNotFound
from flask_restplus import Resource
from sqlalchemy.orm.exc import NoResultFound

from ASB_app.utils import PaginationMixin
from ASB_app.routes import search_parser, csv_columns_parser, used_hints_parser, pagination_parser, search_parser_tsv,\
    browse_parser

from ASB_app.releases import Release, get_release_by_version

from ASB_app import app


def set_release_service(release_service):
    def wrapper(cls):
        cls.release_service = release_service
        return cls
    return wrapper


@app.route('/api', strict_slashes=False)
def get_api_page():
    return render_template('api_page.html', releases=[release for release in Release.__subclasses__()])


@app.route('/sitemap/v<version>/tfs', strict_slashes=False)
def get_tf_page(version):
    try:
        url_release = get_release_by_version(version)
    except ReleaseNotFound as e:
        return {'message': '{}'.format(e)}, 404
    release_service = ReleaseService(url_release)
    return render_template('tf_page.html', tfs=release_service.get_tf_links(), release_name=url_release.name)


@app.route('/sitemap/v<version>/snps', strict_slashes=False)
@app.route('/sitemap/v<version>/snps/<int:page>', strict_slashes=False)
def get_snp_page(version, page=None):
    if page is None:
        page = 0
    if page < 0:
        return 'Page not found', 404
    try:
        url_release = get_release_by_version(version)
    except ReleaseNotFound as e:
        return {'message': '{}'.format(e)}, 404
    release_service = ReleaseService(url_release)
    snp_links, pages = release_service.get_snp_links(page)
    if not snp_links:
        return 'Page not found', 404
    return render_template('snp_page.html', snps=snp_links, release=url_release, pages=pages, current_page=page)


for release in Release.__subclasses__():
    api = release.api

    snp_nsp = api.namespace('SNPs', path='/snps', description='Access to Single Nucleotide Polymorphisms')
    search_nsp = api.namespace('Search', path='/search', description='Search SNPs')
    browse_nsp = api.namespace('Browse', path='/browse', description='Catalog of TFs and Cell types')

    release_service = ReleaseService(release)
    release_serializers = ReleaseSerializers(release)

    @snp_nsp.route('/<int:rs_id>/<string:alt>')
    @set_release_service(release_service)
    class SNPItem(Resource):
        @api.marshal_with(release_serializers.rs_snp_model_full)
        def get(self, rs_id, alt):
            """
            Get complete imformation about an SNP by rs-ID and alt allele
            """
            try:
                return self.release_service.get_full_snp(rs_id, alt)
            except NoResultFound:
                api.abort(404)


    @snp_nsp.route('/<string:chromosome>/<int:position>/<string:alt>/TF/<int:ag_id>')
    @set_release_service(release_service)
    class TFAgrSnpItem(Resource):
        @api.marshal_with(release_serializers.tf_snp_model_full)
        def get(self, chromosome, position, alt, ag_id):
            """
            Get individual SNPs that support a given SNP in aggregation by TF or cell type
            """
            try:
                return self.release_service.get_aggregated_snp(chromosome, position, alt, 'TF', ag_id)
            except NoResultFound:
                api.abort(404)


    @snp_nsp.route('/<string:chromosome>/<int:position>/<string:alt>/CL/<int:ag_id>')
    @set_release_service(release_service)
    class CLAgrSNPItem(Resource):
        @api.marshal_with(release_serializers.cl_snp_model_full)
        def get(self, chromosome, position, alt, ag_id):
            """
            Get individual SNPs that support a given SNP in aggregation by cell type
            """
            try:
                return self.release_service.get_aggregated_snp(chromosome, position, alt, 'CL', ag_id)
            except NoResultFound:
                api.abort(404)


    @search_nsp.route('/snps/rs/<int:rs_id>')
    @set_release_service(release_service)
    class SearchSNPByIdCollection(Resource, PaginationMixin):
        BaseEntity = release_service.SNP
        used_release = release

        @api.marshal_with(release_serializers.search_results_model, PaginationMixin)
        @api.expect(pagination_parser)
        def get(self, rs_id):
            """
            Get all SNPs by rs-ID short info
            """
            all_args = pagination_parser.parse_args()
            filters = self.release_service.get_filters_by_rs_id(rs_id) + \
                self.release_service.get_filters_by_fdr(default_fdr_tr(int(self.used_release.version)))
            result = self.paginate(all_args, extra_filters=filters)
            return {'results': result, 'total': self.items_count(extra_filters=filters)}


    @search_nsp.route('/snps/gene_id/<string:gene_id>')
    @set_release_service(release_service)
    class SearchSNPByGeneIdCollection(Resource, PaginationMixin):
        BaseEntity = release_service.SNP
        used_release = release

        @api.marshal_with(release_serializers.gene_search_results_model)
        @api.expect(pagination_parser)
        def get(self, gene_id):
            """
            Get all SNPs by gene encode id short info
            """
            all_args = pagination_parser.parse_args()
            gene = self.release_service.get_gene_by_id(gene_id)
            if gene is None:
                return {'results': [], 'gene': None, 'total': 0}
            gene.locus_start, gene.locus_end = self.release_service.get_gene_locus(gene)
            filters = self.release_service.get_filters_by_gene(gene) + \
                self.release_service.get_filters_by_fdr(default_fdr_tr(int(self.used_release.version)))
            result = self.paginate(all_args, extra_filters=filters)

            return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


    @search_nsp.route('/snps/gene_name/<string:gene_name>')
    @set_release_service(release_service)
    class SearchSNPByGeneNameCollection(Resource, PaginationMixin):
        BaseEntity = release_service.SNP
        used_release = release

        @api.marshal_with(release_serializers.gene_search_results_model)
        @api.expect(pagination_parser)
        def get(self, gene_name):
            """
            Get all SNPs by gene name short info
            """
            all_args = pagination_parser.parse_args()
            gene = self.release_service.get_gene_by_name(gene_name)
            if gene is None:
                return {'results': [], 'gene': None, 'total': 0}
            gene.locus_start, gene.locus_end = self.release_service.get_gene_locus(gene)
            filters = self.release_service.get_filters_by_gene(gene) + \
                self.release_service.get_filters_by_fdr(default_fdr_tr(int(self.used_release.version)))
            result = self.paginate(all_args, extra_filters=filters)

            return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


    if int(release.version) >= 3:
        @search_nsp.route('/snps/eqtl_gene_id/<string:gene_id>')
        @set_release_service(release_service)
        class SearchEQTLSNPByGeneIdCollection(Resource, PaginationMixin):
            BaseEntity = release_service.SNP
            used_release = release

            @api.marshal_with(release_serializers.gene_search_results_model)
            @api.expect(pagination_parser)
            def get(self, gene_id):
                """
                Get all SNPs by eqtl target gene encode id short info
                """
                all_args = pagination_parser.parse_args()
                gene = self.release_service.get_gene_by_id(gene_id)
                if gene is None:
                    return {'results': [], 'gene': None, 'total': 0}
                gene.locus_start, gene.locus_end = self.release_service.get_gene_locus(gene)
                filters = self.release_service.get_filters_by_eqtl_gene(gene) + \
                    self.release_service.get_filters_by_fdr(default_fdr_tr(int(self.used_release.version)))
                result = self.paginate(all_args, extra_filters=filters)

                return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


        @search_nsp.route('/snps/eqtl_gene_name/<string:gene_name>')
        @set_release_service(release_service)
        class SearchEQTLSNPByGeneNameCollection(Resource, PaginationMixin):
            BaseEntity = release_service.SNP
            used_release = release

            @api.marshal_with(release_serializers.gene_search_results_model)
            @api.expect(pagination_parser)
            def get(self, gene_name):
                """
                Get all SNPs by eqtl target gene name short info
                """
                all_args = pagination_parser.parse_args()
                gene = self.release_service.get_gene_by_name(gene_name)
                if gene is None:
                    print('None')
                    return {'results': [], 'gene': None, 'total': 0}
                gene.locus_start, gene.locus_end = self.release_service.get_gene_locus(gene)
                filters = self.release_service.get_filters_by_eqtl_gene(gene) + \
                    self.release_service.get_filters_by_fdr(default_fdr_tr(int(self.used_release.version)))
                result = self.paginate(all_args, extra_filters=filters)

                return {'results': result, 'gene': gene, 'total': self.items_count(extra_filters=filters)}


    @search_nsp.route('/snps/advanced')
    @set_release_service(release_service)
    class AdvancedSearchSNP(Resource, PaginationMixin):
        BaseEntity = release_service.SNP
        used_release = release

        @api.marshal_with(release_serializers.search_results_model)
        @api.expect(search_parser)
        def get(self):
            """
            Get all SNPs with advanced filters short info
            """
            all_args = search_parser.parse_args()
            try:
                filters = self.release_service.construct_advanced_filters(all_args)
                result = self.paginate(all_args, extra_filters=filters)
                return {'results': result, 'total': self.items_count(extra_filters=filters)}
            except ParsingError:
                api.abort(400)


    @search_nsp.route('/snps/advanced/tsv')
    @set_release_service(release_service)
    class AdvancedSearchSNPCSV(Resource):
        @api.expect(search_parser_tsv)
        def get(self):
            """
            Get all SNPs with advanced filters short info in tsv file
            """
            all_args = search_parser_tsv.parse_args()
            try:
                return self.release_service.get_snps_by_advanced_filters_tsv(all_args)
            except ParsingError:
                api.abort(400)


    if int(release.version) >= 3:
        @search_nsp.route('/snps/advanced/targets')
        @set_release_service(release_service)
        class AdvancedSearchTargets(Resource):
            @api.expect(search_parser)
            def get(self):
                """
                Get all SNPs with advanced filters short info in tsv file annotated with eQTL targets
                """
                all_args = search_parser.parse_args()
                try:
                    return self.release_service.get_snps_by_advanced_filters_tsv_with_targets(all_args)
                except ParsingError:
                    api.abort(400)

    if release.name != 'dnase':
        @browse_nsp.route('/tf')
        @set_release_service(release_service)
        class TransctiptionFactorBrowse(Resource, PaginationMixin):
            BaseEntity = release_service.TranscriptionFactor
            used_release = release

            @api.marshal_with(release_serializers.transcription_factor_browse_model)
            @api.expect(browse_parser)
            def get(self):
                """
                Get the list of transcription factors available in the database
                """
                args = pagination_parser.parse_args()
                filters = (self.BaseEntity.aggregated_snps_count > 0, )
                if 'filter' in args:
                    filter_exp = args.pop('filter')
                    filters += (self.BaseEntity.name.like(filter_exp))
                result = self.paginate(args, extra_filters=filters, default_order_clauses=(desc(self.BaseEntity.aggregated_snps_count), ))
                total = self.items_count(extra_filters=filters)
                return {'results': result, 'total': total}


    @browse_nsp.route('/cl')
    @set_release_service(release_service)
    class CellLineBrowse(Resource, PaginationMixin):
        BaseEntity = release_service.CellLine
        used_release = release

        @api.marshal_with(release_serializers.cell_line_browse_model)
        @api.expect(browse_parser)
        def get(self):
            """
            Get the list of cell types available in the database
            """
            args = pagination_parser.parse_args()
            filters = (self.BaseEntity.aggregated_snps_count > 0, )
            if 'filter' in args:
                filter_exp = args.pop('filter')
                filters += (self.BaseEntity.name.like(filter_exp))
            result = self.paginate(args, extra_filters=filters, default_order_clauses=(desc(self.BaseEntity.aggregated_snps_count), ))
            total = self.items_count(extra_filters=filters)
            return {'results': result, 'total': total}


    @browse_nsp.route('/total')
    @api.hide
    @set_release_service(release_service)
    class FrontPageStatistics(Resource):
        @api.marshal_with(release_serializers.frontpage_statistics_model)
        def get(self):
            return self.release_service.get_overall_statistics()


    if release.name != 'dnase':
        @search_nsp.route('/tf/hint')
        @api.hide
        @set_release_service(release_service)
        class TranscriptionFactorHint(Resource):
            @api.expect(used_hints_parser)
            @api.marshal_list_with(release_serializers.transcription_factor_model)
            def get(self):
                args = used_hints_parser.parse_args()
                return self.release_service.get_hints('TF', args.get('search', ''), args.get('options', []))


    @search_nsp.route('/cl/hint')
    @api.hide
    @set_release_service(release_service)
    class CellLineHint(Resource):
        @api.expect(used_hints_parser)
        @api.marshal_list_with(release_serializers.cell_line_model)
        def get(self):
            args = used_hints_parser.parse_args()
            return self.release_service.get_hints('CL', args.get('search', ''), args.get('options', []))


    @search_nsp.route('/gene_name/hint')
    @api.hide
    @set_release_service(release_service)
    class GeneNameHint(Resource):
        @api.expect(used_hints_parser)
        @api.marshal_list_with(release_serializers.gene_model)
        def get(self):
            args = used_hints_parser.parse_args()
            return self.release_service.get_hints_for_gene_name(args.get('search', ''))


    @search_nsp.route('/eqtl_gene_name/hint')
    @api.hide
    @set_release_service(release_service)
    class GeneNameHint(Resource):
        @api.expect(used_hints_parser)
        @api.marshal_list_with(release_serializers.gene_model)
        def get(self):
            args = used_hints_parser.parse_args()
            return self.release_service.get_hints_for_eqtl_gene_name(args.get('search', ''))


    @snp_nsp.route('/<int:rs_id>/<string:alt>/<string:what_for>/tsv')
    @set_release_service(release_service)
    class SNPItemCSV(Resource):
        @api.expect(csv_columns_parser)
        def get(self, rs_id, alt, what_for):
            """
            Get complete information about a SNP by rs-ID and Alternative allele in TSV file
            """
            args = csv_columns_parser.parse_args()
            try:
                return self.release_service.get_full_snp_tsv(what_for, rs_id, alt, args['columns'])
            except NoResultFound:
                api.abort(404)
