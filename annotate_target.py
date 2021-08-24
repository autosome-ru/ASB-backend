from sqlalchemy.orm import aliased

from ASB_app import *
from ASB_app.releases import ReleaseSusan as release
from ASB_app.service.adastra_service import ReleaseService
from ASB_app.utils.statistics import get_corresponding_fdr_classes, get_corresponding_es_classes

SNP = release.SNP

factor_name = 'FOXA2_HUMAN'
release_service = ReleaseService(release)
fdr = '0.1'
es = '0'


def get_snps(filters_object):
    # file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
    # csv_writer = csv.writer(file, dialect=TsvDialect)

    headers = ['Chromosome', 'Position', 'Ref', 'Alt', 'rsID']
    names = dict(zip(headers, ['chromosome', 'position', 'ref', 'alt', 'rs_id']))

    join_tuples = []
    additional_columns = []
    query_args = [release.SNP]
    aliases = dict()

    if not filters_object['fdr']:
        filters_object['fdr'] = 0.1
    if not filters_object['es']:
        filters_object['es'] = 0

    for what_for in ('TF', 'CL'):
        aggregation_class = {'TF': release.TranscriptionFactor, 'CL': release.CellLine}[what_for]
        aggregated_snp_class = {'TF': release.TranscriptionFactorSNP, 'CL': release.CellLineSNP}[what_for]
        id_field = {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
        agr_snp_class_alias = aliased(aggregated_snp_class)
        agr_class_alias = aliased(aggregation_class)
        query_args.append(release.db.func.group_concat(agr_class_alias.name.distinct()))
        headers.append('{}-ASBs'.format({'TF': 'TF', 'CL': 'Cell type'}[what_for]))


        join_tuples += [
            (
                agr_snp_class_alias,
                (agr_snp_class_alias.chromosome == release.SNP.chromosome) &
                (agr_snp_class_alias.position == release.SNP.position) &
                (agr_snp_class_alias.alt == release.SNP.alt) &
                (agr_snp_class_alias.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr']))) &
                (agr_snp_class_alias.es_class.in_(get_corresponding_es_classes(filters_object['es']))),
            ),
            (
                agr_class_alias,
                getattr(agr_class_alias, id_field) == getattr(agr_snp_class_alias, id_field),
            )
        ]

    for what_for in ('TF', 'CL'):
        filter_object_key = {'TF': 'transcription_factors', 'CL': 'cell_types'}[what_for]
        aggregation_class = {'TF': release.TranscriptionFactor, 'CL': release.CellLine}[what_for]
        aggregated_snp_class = {'TF': release.TranscriptionFactorSNP, 'CL': release.CellLineSNP}[what_for]
        id_field = {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
        if filters_object[filter_object_key]:
            for name in filters_object[filter_object_key]:
                aliases[name] = aliased(aggregated_snp_class)
                aggregation_entity = aggregation_class.query.filter_by(name=name).one_or_none()
                if not aggregation_entity:
                    continue
                aggregation_id = getattr(aggregation_entity, id_field)
                join_tuples.append(
                    (aliases[name],
                     (getattr(
                         aliases[name],
                         {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
                     ) == aggregation_id) &
                     (aliases[name].chromosome == release.SNP.chromosome) &
                     (aliases[name].position == release.SNP.position) &
                     (aliases[name].alt == release.SNP.alt))
                )
                for field, label in [
                    ('log_p_value_ref', '{}_FDR_Ref'.format(name)),
                    ('log_p_value_alt', '{}_FDR_Alt'.format(name)),
                    ('es_ref', '{}_Effect_Size_Ref'.format(name)),
                    ('es_alt', '{}_Effect_Size_Alt'.format(name)),
                ]:
                    headers.append(label)
                    additional_columns.append(
                        release.db.func.coalesce(getattr(aliases[name], field)).label(label))

    Target = aliased(release.Gene)
    join_tuples.append((Target, release.SNP.target_genes))
    additional_columns.append(release.db.func.group_concat(Target.gene_name.distinct()).label('targets'))
    headers.append('Targets')

    found_snps = release.session.query(*query_args)
    found_snps = found_snps.filter(*release_service.construct_advanced_filters(filters_object))
    for cls, condition in join_tuples:
        found_snps = found_snps.join(cls, condition, isouter=True)
    found_snps = found_snps.add_columns(*additional_columns)
    found_snps = found_snps.group_by(release.SNP)

    with open('D:\Sashok/FOXA2_targets.tsv', 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for tup in found_snps.join():
            snp = tup[0]
            columns = tup[1:]
            line = \
                [getattr(snp, names[header]) if names[header] != 'rs_id' else 'rs' + str(getattr(snp, names[header]))
                 for header in headers[:len(names.keys())]] + list(columns)
            f.write('\t'.join(map(str, line)) + '\n')

get_snps({'transcription_factors': [factor_name], 'cell_types': None, 'fdr': fdr, 'es': es, 'motif_concordance': None, 'chromosome': None, 'position': None, 'alt': None, 'phenotype_databases': None})