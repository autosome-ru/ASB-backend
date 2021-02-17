from flask_restplus import inputs
from ASB_app.constants import chromosomes, fdr_choices
from werkzeug.datastructures import FileStorage

from ASB_app.releases import Release, current_release

for release in Release.__subclasses__():
    api = release.api

    file_parser = api.parser()
    file_parser.add_argument('file', type=FileStorage, location='files')

    pagination_parser = api.parser()
    pagination_parser.add_argument('page', type=inputs.positive, help='Page number', default=1)
    pagination_parser.add_argument('size', type=inputs.natural, help='Items per page or 0 for all items', default=0)
    pagination_parser.add_argument('offset', type=inputs.natural, help='Skip first N items', default=0)
    pagination_parser.add_argument('filter', help='Comma-separated filters: field1 EQ "value1", field2 GE "value2"')
    pagination_parser.add_argument('order_by', help='ORDER BY criterion: "field1", "-field2"')

    search_parser = pagination_parser.copy()
    if int(release.version) >= 3:
        search_parser.add_argument('fdr', help='FDR threshold', default='0.05', choices=fdr_choices)
    search_parser.add_argument('cell_types', action='split', help='Comma-separated list of cell types, search SNPs ASB for every cell type scpecified')
    search_parser.add_argument('transcription_factors', action='split', help='Comma-separated list of cell types, search SNPs ASB for every cell type scpecified')
    search_parser.add_argument('chromosome', choices=chromosomes, help='Search only SNPs on the specified chromosome')
    search_parser.add_argument('start', type=inputs.positive, help='Search SNPs in interval from specified position, Requires "chromosome" and "end"')
    search_parser.add_argument('end', type=inputs.positive, help='Search SNPs in interval to specified position, Requiers "chromosome" and "start"')
    search_parser.add_argument('phenotype_databases', action='split', help='Comma-separated list of databases, possible choices {grasp, ebi, clinvar, phewas, finemapping, QTL}, earch SNPs that have phenotype associations in all specified databases')
    search_parser.add_argument('motif_concordance', action='split', help='Comma-separated list of motif concordance values, possible choices {Concordant. Discordant, Weak Concordant, Weak Discordant}, if no TF specified will search SNPs with any TF having any of the specified concordance valuese, else only SNPs ASB for the specified TFs and any of the specified concordance values for these TFs')

    used_hints_parser = api.parser()
    used_hints_parser.add_argument('options', action='split')
    used_hints_parser.add_argument('search')

    csv_columns_parser = api.parser()
    csv_columns_parser.add_argument('columns', action='split', required=True)
    csv_columns_parser.add_argument('filter')

result_param_parser = current_release.api.parser()
result_param_parser.add_argument('result_param', choices=('tf', 'cl', 'tf_sum', 'cl_sum'), default='tf')
result_param_parser.add_argument('format', choices=('json', 'tsv'), default='json')
result_param_parser.add_argument('limit', type=inputs.positive, default=0)

fdr_parser = current_release.api.parser()
fdr_parser.add_argument('fdr', help='FDR threshold', default='0.05', choices=fdr_choices)


def parse_fdr(fdr):
    try:
        fdr = float(fdr)
        assert 0 < fdr <= 0.25
    except (ValueError, AssertionError):
        api.abort(400)
    return fdr
