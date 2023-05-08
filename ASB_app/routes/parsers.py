from flask_restplus import inputs
from ASB_app.constants import chromosomes, fdr_choices, es_choices, background_choices, nucleotides
from werkzeug.datastructures import FileStorage

from ASB_app.releases import Release, current_release

for release in Release.__subclasses__():
    api = release.api

    file_parser = api.parser()
    file_parser.add_argument('file', type=FileStorage, location='files')

    pagination_parser = api.parser()
    pagination_parser.add_argument('page', type=inputs.positive, help='Page number', default=1)
    pagination_parser.add_argument('size', type=inputs.int_range(1, 1200), help='Items per page, ≤1000')
    pagination_parser.add_argument('offset', type=inputs.natural, help='Skip first N items', default=0)
    pagination_parser.add_argument('filter', help='Comma-separated filters: field1 EQ "value1", field2 GE "value2"')
    pagination_parser.add_argument('order_by', help='ORDER BY criterion: "field1", "-field2"')

    search_parser = pagination_parser.copy()
    search_parser.add_argument('fdr', help='FDR threshold', choices=fdr_choices)
    search_parser.add_argument('es', help='Effect size threshold', choices=es_choices)
    search_parser.add_argument('atac', type=lambda value: value.split('@'),
                               help='"@" separated list of cell types, search SNPs ASB for every cell type specified')
    search_parser.add_argument('dnase', action='split', help='"@" separated list of cell types, search SNPs ASB for every cell type specified')
    search_parser.add_argument('faire', action='split',
                               help='"@" separated list of cell types, search SNPs ASB for every cell type specified')
    search_parser.add_argument('chromosome', choices=chromosomes, help='Search only SNPs on the specified chromosome')
    search_parser.add_argument('start', type=inputs.positive, help='Search SNPs in interval from specified position, Requires "chromosome" and "end", 1-based')
    search_parser.add_argument('end', type=inputs.positive, help='Search SNPs in interval to specified position, Requiers "chromosome" and "start", 1-based')
    search_parser.add_argument('phenotype_databases', action='split', help='Comma-separated list of databases, possible choices {grasp, ebi, clinvar, phewas, finemapping, QTL}, earch SNPs that have phenotype associations in all specified databases')

    search_parser_tsv = search_parser.copy()
    search_parser_tsv.replace_argument('size', type=inputs.natural, help='Items per page', default=0)

    browse_parser = pagination_parser.copy()
    browse_parser.add_argument('regexp', help="Regular expression for filtering resulting cell types' names")

    aggregated_snp_parser = api.parser()
    aggregated_snp_parser.add_argument('chromosome', choices=chromosomes, required=True, help='Chromosome')
    aggregated_snp_parser.add_argument('position', type=inputs.positive, required=True, help='Genome position in chromosome, 1-based')
    aggregated_snp_parser.add_argument('alt', choices=nucleotides, required=True, help='Alternative alllele')
    aggregated_snp_parser.add_argument('aggregation_name', required=True, help='Cell type name')

    used_hints_parser = api.parser()
    used_hints_parser.add_argument('options', action='split')
    used_hints_parser.add_argument('search')

    csv_columns_parser = api.parser()
    csv_columns_parser.add_argument('columns', action='split', required=True)
    csv_columns_parser.add_argument('filter')
#
# result_param_parser = current_release.api.parser()
# result_param_parser.add_argument('result_param', choices=('tf', 'cl', 'tf_sum', 'cl_sum', 'all', 'not_found'), default='all')
# result_param_parser.add_argument('format', choices=('json', 'tsv'), default='json')
# result_param_parser.add_argument('page', type=inputs.positive, help='Page number', default=1)
# result_param_parser.add_argument('size', type=inputs.int_range(1, 1000), help='Items per page, ≤1000')
# result_param_parser.add_argument('filter', help='comma-separated tf or cell-line names, "-" for exclude', action='split')
# result_param_parser.add_argument('order_by', help='order by column, "-" for desc')
#
# thresholds_parser = current_release.api.parser()
# thresholds_parser.add_argument('fdr', help='FDR threshold', default='0.05', choices=fdr_choices)
# thresholds_parser.add_argument('background', help='Background comparison SNP set', default='WG', choices=background_choices)
# # thresholds_parser.add_argument('es', help='Effect size threshold', default='0.6', choices=es_choices)
