from sqlalchemy.orm import aliased

from ASB_app.models import abstract_models, abstract_models_dnase
from ASB_app.utils.aggregates import db_name_property_dict, TsvDialect
from sqlalchemy import not_, or_
import csv
import tempfile
from flask import send_file
import numpy as np

from ASB_app.constants import stats_dict

from math import ceil

from ASB_app.utils.statistics import get_corresponding_fdr_classes


class ReleaseService:
    def __init__(self, release):
        self.release = release
        if release.name != 'dnase':
            for model in abstract_models:
                setattr(self, model.__name__, getattr(release, model.__name__))
        else:
            for model in abstract_models_dnase:
                setattr(self, model.__name__, getattr(release, model.__name__))

    def get_filters_by_fdr(self, fdr):
        if int(self.release.version) >= 3:
            return (self.SNP.fdr_class.in_(get_corresponding_fdr_classes(fdr)),)
        else:
            return tuple()

    def generate_tf_link(self, tf_name):
        return 'https://adastra.autosome.ru/{}/search/advanced?tf={}'.format(self.release.name, tf_name)

    def generate_snp_name(self, snp):
        return 'rs{0.rs_id}:{0.ref}>{0.alt}'.format(snp)

    def generate_snp_link(self, snp):
        return 'https://adastra.autosome.ru/{0.name}/snps/rs{1.rs_id}/{1.alt}'.format(self.release, snp)

    def get_tf_links(self):
        return [{'name': tf.name, 'link': self.generate_tf_link(tf.name)} for tf in
                self.TranscriptionFactor.query.filter(self.TranscriptionFactor.aggregated_snps_count > 0).order_by(
                    self.TranscriptionFactor.name)]

    def get_snp_links(self, page, on_page=20000):
        return (
            [{'name': self.generate_snp_name(snp), 'link': self.generate_snp_link(snp)} for snp in
             self.SNP.query.order_by(self.SNP.rs_id).offset(on_page * page).limit(on_page)],
            ceil(self.SNP.query.count() / on_page)
        )

    def get_filters_by_rs_id(self, rs_id):
        return (self.SNP.rs_id == rs_id,)

    def get_full_snp(self, rs_id, alt):
        return self.SNP.query.filter(
            (self.SNP.rs_id == rs_id) &
            (self.SNP.alt == alt)
        ).one()

    @staticmethod
    def format_header(header):
        format_dict = {
            'transcription_factor.name': 'transcription_factor',
            'cell_line.name': 'cell_type'
        }
        if header in format_dict:
            return format_dict[header]
        return header

    def get_full_snp_tsv(self, what_for, rs_id, alt, headers):
        agg_snps = getattr(self.get_full_snp(rs_id, alt), what_for + '_aggregated_snps')

        file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
        csv_writer = csv.writer(file, dialect=TsvDialect)

        csv_writer.writerow([self.format_header(h) for h in headers])

        getter = lambda snp, header: getter(getattr(snp, header[:header.find('.')]),
                                            header[header.find('.') + 1:]) if '.' in header else getattr(snp, header)

        for snp in agg_snps:
            csv_writer.writerow([getter(snp, header) for header in headers])

        file.flush()

        return send_file(
            file.name,
            cache_timeout=0,
            mimetype="text/tsv",
            as_attachment=True
        )

    def get_gene_by_name(self, gene_name):
        return self.Gene.query.filter_by(gene_name=gene_name).one_or_none()

    def get_gene_by_id(self, gene_id):
        return self.Gene.query.get(gene_id)

    def get_filters_by_gene(self, gene):
        if gene.orientation:
            return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos - 500, 1),
                                                                                     gene.end_pos)
        else:
            return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos, 1),
                                                                                     gene.end_pos + 500)

    def construct_advanced_filters(self, filters_object):
        filters = []
        if int(self.release.version) >= 3:
            if filters_object['transcription_factors']:
                filters += [self.SNP.tf_aggregated_snps.any(
                    (self.TranscriptionFactorSNP.tf_id == getattr(self.TranscriptionFactor.query.filter(
                        self.TranscriptionFactor.name == tf_name
                    ).one_or_none(), 'tf_id', None)) &
                    (self.TranscriptionFactorSNP.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr']))))
                    for tf_name in filters_object['transcription_factors']]

            if filters_object['cell_types']:
                filters += [self.SNP.cl_aggregated_snps.any(
                    (self.CellLineSNP.cl_id == getattr(self.CellLine.query.filter(
                        self.CellLine.name == cl_name
                    ).one_or_none(), 'cl_id', None)) &
                    (self.CellLineSNP.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr']))))
                    for cl_name in filters_object['cell_types']]

            if not filters_object['transcription_factors'] and not filters_object['cell_types']:
                filters += self.get_filters_by_fdr(filters_object['fdr'])

            if filters_object['motif_concordance']:
                search_null = False
                if 'None' in filters_object['motif_concordance']:
                    search_null = True
                    filters_object['motif_concordance'] = [x for x in filters_object['motif_concordance'] if
                                                           x != 'None']
                if filters_object['transcription_factors']:
                    filters += [self.SNP.tf_aggregated_snps.any(
                        (self.TranscriptionFactorSNP.motif_concordance.in_(filters_object['motif_concordance']) |
                         (self.TranscriptionFactorSNP.motif_concordance.is_(None) if search_null else False)) &
                        self.TranscriptionFactorSNP.transcription_factor.has(
                            self.TranscriptionFactor.name.in_(filters_object['transcription_factors'])
                        ) & (self.TranscriptionFactorSNP.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr'])))
                    )]
                else:
                    filters += [self.SNP.tf_aggregated_snps.any(
                        (self.TranscriptionFactorSNP.motif_concordance.in_(filters_object['motif_concordance']) |
                         (self.TranscriptionFactorSNP.motif_concordance.is_(None) if search_null else False)) &
                        (self.TranscriptionFactorSNP.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr'])))
                    ) | (~self.SNP.tf_aggregated_snps.any() if search_null else False)]
        else:
            if filters_object['transcription_factors']:
                filters += [self.SNP.tf_aggregated_snps.any(
                    self.TranscriptionFactorSNP.tf_id == getattr(self.TranscriptionFactor.query.filter(
                        self.TranscriptionFactor.name == tf_name
                    ).one_or_none(), 'tf_id', None))
                    for tf_name in filters_object['transcription_factors']]

            if filters_object['cell_types']:
                filters += [self.SNP.cl_aggregated_snps.any(
                    self.CellLineSNP.cl_id == getattr(self.CellLine.query.filter(
                        self.CellLine.name == cl_name
                    ).one_or_none(), 'cl_id', None))
                    for cl_name in filters_object['cell_types']]

            if filters_object['motif_concordance']:
                search_null = False
                if 'None' in filters_object['motif_concordance']:
                    search_null = True
                    filters_object['motif_concordance'] = [x for x in filters_object['motif_concordance'] if
                                                           x != 'None']
                if filters_object['transcription_factors']:
                    filters += [self.SNP.tf_aggregated_snps.any(
                        (self.TranscriptionFactorSNP.motif_concordance.in_(filters_object['motif_concordance']) |
                         (self.TranscriptionFactorSNP.motif_concordance.is_(None) if search_null else False)) &
                        self.TranscriptionFactorSNP.transcription_factor.has(
                            self.TranscriptionFactor.name.in_(filters_object['transcription_factors'])
                        )
                    )]
                else:
                    filters += [self.SNP.tf_aggregated_snps.any(
                        self.TranscriptionFactorSNP.motif_concordance.in_(filters_object['motif_concordance']) |
                        (self.TranscriptionFactorSNP.motif_concordance.is_(None) if search_null else False)
                    ) | (~self.SNP.tf_aggregated_snps.any() if search_null else False)]

        if filters_object['chromosome']:
            if not filters_object['start'] or not filters_object['end']:
                filters += [self.SNP.chromosome == filters_object['chromosome']]
            else:
                filters += [self.SNP.chromosome == filters_object['chromosome'],
                            self.SNP.position.between(filters_object['start'], filters_object['end'])]

        if filters_object['phenotype_databases']:
            filters += [or_(*(getattr(self.SNP, db_name_property_dict[phenotype_db])
                              for phenotype_db in filters_object['phenotype_databases']))]

        return filters

    def get_snps_by_advanced_filters_tsv(self, filters_object):
        file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
        csv_writer = csv.writer(file, dialect=TsvDialect)

        headers = ['Chromosome', 'Position', 'Ref', 'Alt', 'rsID']
        names = dict(zip(headers, ['chromosome', 'position', 'ref', 'alt', 'rs_id']))

        join_tuples = []
        additional_columns = []
        query_args = [self.SNP]
        aliases = dict()

        for what_for in ('TF', 'CL'):
            aggregation_class = {'TF': self.TranscriptionFactor, 'CL': self.CellLine}[what_for]
            aggregated_snp_class = {'TF': self.TranscriptionFactorSNP, 'CL': self.CellLineSNP}[what_for]
            id_field = {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
            query_args.append(self.release.db.func.group_concat(aggregation_class.name.distinct()))
            headers.append('{}-ASBs'.format({'TF': 'TF', 'CL': 'Cell type'}[what_for]))
            join_tuples += [
                (
                    aggregated_snp_class,
                    (aggregated_snp_class.chromosome == self.SNP.chromosome) &
                    (aggregated_snp_class.position == self.SNP.position) &
                    (aggregated_snp_class.alt == self.SNP.alt),
                ),
                (
                    aggregation_class,
                    getattr(aggregation_class, id_field) == getattr(aggregated_snp_class, id_field),
                )
            ]

        for what_for in ('TF', 'CL'):
            filter_object_key = {'TF': 'transcription_factors', 'CL': 'cell_types'}[what_for]
            aggregation_class = {'TF': self.TranscriptionFactor, 'CL': self.CellLine}[what_for]
            aggregated_snp_class = {'TF': self.TranscriptionFactorSNP, 'CL': self.CellLineSNP}[what_for]
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
                         (aliases[name].chromosome == self.SNP.chromosome) &
                         (aliases[name].position == self.SNP.position) &
                         (aliases[name].alt == self.SNP.alt))
                    )
                    for field, label in [
                        ('log_p_value_ref', '{}_FDR_Ref'.format(name)),
                        ('log_p_value_alt', '{}_FDR_Alt'.format(name)),
                        ('es_ref', '{}_Effect_Size_Ref'.format(name)),
                        ('es_alt', '{}_Effect_Size_Alt'.format(name)),
                    ]:
                        headers.append(label)
                        additional_columns.append(getattr(aggregated_snp_class, field).label(label))

        found_snps = self.release.session.query(*query_args)
        found_snps = found_snps.filter(*self.construct_advanced_filters(filters_object))
        for cls, condition in join_tuples:
            found_snps = found_snps.join(cls, condition, isouter=True)
        for column in additional_columns:
            found_snps = found_snps.add_column(column)
        found_snps = found_snps.group_by(self.SNP)

        csv_writer.writerow(headers)
        for tup in found_snps:
            snp = tup[0]
            columns = tup[1:]
            csv_writer.writerow(
                [getattr(snp, names[header]) if names[header] != 'rs_id' else 'rs' + str(getattr(snp, names[header]))
                 for header in headers[:len(names.keys())]] + list(columns))

        file.flush()
        return send_file(
            file.name,
            cache_timeout=0,
            mimetype="text/tsv",
            as_attachment=True
        )

    def get_hints(self, what_for, in_str, used_options):
        cls = {'TF': self.TranscriptionFactor, 'CL': self.CellLine}[what_for]
        if int(self.release.version) >= 3:
            filters = (((cls.name.like(in_str),) if in_str else ()) +
                       ((not_(cls.name.in_(used_options)),) if used_options else ()) +
                       (cls.aggregated_snps_count005,))
        else:
            filters = (((cls.name.like(in_str),) if in_str else ()) +
                       ((not_(cls.name.in_(used_options)),) if used_options else ()) +
                       (cls.aggregated_snps_count,))
        return cls.query.filter(*filters).order_by(cls.aggregated_snps_count.desc()).limit(3).all()

    def get_hints_for_gene_name(self, in_str):
        filters = (self.Gene.gene_name.like(in_str),) if in_str else ()
        return self.Gene.query.filter(*filters).order_by(self.Gene.snps_count.desc()).limit(3).all()

    def get_overall_statistics(self):
        if int(self.release.version) >= 3:
            return {
                'transcription_factors_count': self.TranscriptionFactor.query.filter(
                    self.TranscriptionFactor.aggregated_snps_count005 > 0).count(),
                'cell_types_count': self.CellLine.query.filter(self.CellLine.aggregated_snps_count005 > 0).count(),
                'snps_count': stats_dict['0.05']['possible_all_asbs_rs'],
            }
        else:
            return {
                'transcription_factors_count': self.TranscriptionFactor.query.filter(
                    self.TranscriptionFactor.aggregated_snps_count > 0).count(),
                'cell_types_count': self.CellLine.query.filter(self.CellLine.aggregated_snps_count > 0).count(),
                'snps_count': self.release.session.query(self.SNP.rs_id).distinct().count(),
            }
