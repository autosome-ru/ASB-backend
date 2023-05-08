from sqlalchemy.orm import aliased

from ASB_app.models import abstract_models
from ASB_app.utils.aggregates import db_name_property_dict, TsvDialect
from sqlalchemy import not_, or_
import csv
import tempfile
from flask import send_file

from ASB_app.constants import stats_dict, default_fdr_tr, default_es_tr

from math import ceil

from ASB_app.utils.statistics import get_corresponding_fdr_classes, get_corresponding_es_classes


class ReleaseService:
    def __init__(self, release):
        self.release = release
        self.what_for = ('dnase', 'atac', 'faire')
        self.filter_object_keys = {x: x for x in self.what_for}

        for model in abstract_models:
            setattr(self, model.__name__, getattr(release, model.__name__))
        self.aggregation_classes = {'dnase': self.Dnase, 'faire': self.Faire, 'atac': self.Atac}
        self.aggregation_snp_classes = {'dnase': self.DnaseSNP, 'faire': self.FaireSNP, 'atac': self.AtacSNP}

    def get_filters_by_fdr(self, fdr):
        if int(self.release.version) >= 3:
            return (self.SNP.fdr_class.in_(get_corresponding_fdr_classes(fdr)),)
        else:
            return tuple()

    def get_filters_by_es(self, es):
        if int(self.release.version) >= 3:
            return (self.SNP.es_class.in_(get_corresponding_es_classes(es)),)
        else:
            return tuple()

    def get_release_name(self):
        rname = self.release.name if self.release.name != 'billcipher' else 'bill-cipher'
        return f'https://adastra.autosome.org/{rname}'


    def generate_tf_link(self, tf_name):
        return '{}/search/advanced?tf={}'.format(self.get_release_name(), tf_name)

    def generate_snp_name(self, snp):
        return 'rs{0.rs_id}:{0.ref}>{0.alt}'.format(snp)

    def generate_snp_link(self, snp):
        rname = self.get_release_name()
        return '{0}/snps/rs{1.rs_id}/{1.alt}'.format(rname, snp)

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

    def get_aggregated_snp(self, chromosome, position, alt, ag_level, ag_name):
        if ag_level not in self.what_for:
            raise ValueError(ag_level)
        AgSNPClass = self.aggregation_snp_classes[ag_level]
        AgClass = self.aggregation_classes[ag_level]
        ag_attr = 'cl_id'
        ag_snp = AgSNPClass.query.filter(
            (AgSNPClass.chromosome == chromosome) &
            (AgSNPClass.position == position) &
            (AgSNPClass.alt == alt)
        ).join(
            AgClass,
            getattr(AgSNPClass, ag_attr) == getattr(AgClass, ag_attr)
        ).filter(
            (AgClass.name == ag_name)
        ).one()
        return ag_snp

    @staticmethod
    def format_header(header):
        format_dict = {
            'dnase.name': 'dnase.cell_type',
            'atac.name': 'atac.cell_type',
            'faire.name': 'faire.cell_type'
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
        return self.Gene.query.filter((self.Gene.gene_id == gene_id) | (self.Gene.gene_id.like('{}.%'.format(gene_id)))
                                      ).one_or_none()

    def get_filters_by_gene(self, gene):
        if gene.orientation:
            return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos - 5000, 1),
                                                                                     gene.end_pos)
        else:
            return self.SNP.chromosome == gene.chromosome, self.SNP.position.between(max(gene.start_pos, 1),
                                                                                     gene.end_pos + 5000)

    def get_filters_by_eqtl_gene(self, gene):
        return (self.SNP.target_genes.any(self.Gene.gene_id == gene.gene_id),)

    def construct_advanced_filters(self, filters_object):
        filters = []
        if not filters_object['fdr']:
            filters_object['fdr'] = default_fdr_tr(int(self.release.version))
        if not filters_object['es']:
            filters_object['es'] = default_es_tr(int(self.release.version))
        if filters_object['transcription_factors']:
            filters += [self.SNP.tf_aggregated_snps.any(
                (self.TranscriptionFactorSNP.tf_id == getattr(self.TranscriptionFactor.query.filter(
                    self.TranscriptionFactor.name == tf_name
                ).one_or_none(), 'tf_id', None)) &
                (self.TranscriptionFactorSNP.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr']))) &
                (self.TranscriptionFactorSNP.es_class.in_(get_corresponding_es_classes(filters_object['es']))))
                for tf_name in filters_object['transcription_factors']]

        if filters_object['cell_types']:
            filters += [self.SNP.cl_aggregated_snps.any(
                (self.CellLineSNP.cl_id == getattr(self.CellLine.query.filter(
                    self.CellLine.name == cl_name
                ).one_or_none(), 'cl_id', None)) &
                (self.CellLineSNP.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr']))) &
                (self.CellLineSNP.es_class.in_(get_corresponding_es_classes(filters_object['es']))))
                for cl_name in filters_object['cell_types']]

        if not filters_object['transcription_factors'] and not filters_object['cell_types']:
            filters += self.get_filters_by_fdr(filters_object['fdr']) + \
                       self.get_filters_by_es(filters_object['es'])

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

        if int(self.release.version) >= 3:
            if not filters_object['fdr']:
                filters_object['fdr'] = default_fdr_tr(int(self.release.version))
            if not filters_object['es']:
                filters_object['es'] = default_es_tr(int(self.release.version))

        for what_for in self.what_for:
            aggregation_class = self.aggregation_classes[what_for]
            aggregated_snp_class = self.aggregation_snp_classes[what_for]
            id_field = 'cl_id'
            agr_snp_class_alias = aliased(aggregated_snp_class)
            agr_class_alias = aliased(aggregation_class)
            query_args.append(self.release.db.func.group_concat(agr_class_alias.name.distinct()))
            headers.append('{}-ASBs'.format({'TF': 'TF', 'CL': 'Cell type'}[what_for]))

            if int(self.release.version) >= 3:
                join_tuples += [
                    (
                        agr_snp_class_alias,
                        (agr_snp_class_alias.chromosome == self.SNP.chromosome) &
                        (agr_snp_class_alias.position == self.SNP.position) &
                        (agr_snp_class_alias.alt == self.SNP.alt) &
                        (agr_snp_class_alias.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr']))) &
                        (agr_snp_class_alias.es_class.in_(get_corresponding_es_classes(filters_object['es']))),
                    ),
                    (
                        agr_class_alias,
                        getattr(agr_class_alias, id_field) == getattr(agr_snp_class_alias, id_field),
                    )
                ]
            else:
                join_tuples += [
                    (
                        agr_snp_class_alias,
                        (agr_snp_class_alias.chromosome == self.SNP.chromosome) &
                        (agr_snp_class_alias.position == self.SNP.position) &
                        (agr_snp_class_alias.alt == self.SNP.alt),
                    ),
                    (
                        agr_class_alias,
                        getattr(agr_class_alias, id_field) == getattr(agr_snp_class_alias, id_field),
                    )
                ]

        for what_for in self.what_for:
            filter_object_key = self.filter_object_keys[what_for]
            aggregation_class = self.aggregation_classes[what_for]
            aggregated_snp_class = self.aggregation_snp_classes[what_for]
            id_field = 'cl_id'
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
                             'cl_id'
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
                        additional_columns.append(
                            self.release.db.func.coalesce(getattr(aliases[name], field)).label(label))

        found_snps = self.release.session.query(*query_args)
        found_snps = found_snps.filter(*self.construct_advanced_filters(filters_object))
        for cls, condition in join_tuples:
            found_snps = found_snps.join(cls, condition, isouter=True)
        found_snps = found_snps.add_columns(*additional_columns)
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

    def get_snps_by_advanced_filters_tsv_with_targets(self, filters_object):
        file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
        csv_writer = csv.writer(file, dialect=TsvDialect)

        headers = ['Chromosome', 'Position', 'Ref', 'Alt', 'rsID']
        names = dict(zip(headers, ['chromosome', 'position', 'ref', 'alt', 'rs_id']))

        join_tuples = []
        additional_columns = []
        query_args = [self.SNP]
        aliases = dict()

        if int(self.release.version) >= 3:
            if not filters_object['fdr']:
                filters_object['fdr'] = default_fdr_tr(int(self.release.version))
            if not filters_object['es']:
                filters_object['es'] = default_es_tr(int(self.release.version))

        for what_for in self.what_for:
            aggregation_class = self.aggregation_classes[what_for]
            aggregated_snp_class = self.aggregation_snp_classes[what_for]
            id_field = 'cl_id'
            agr_snp_class_alias = aliased(aggregated_snp_class)
            agr_class_alias = aliased(aggregation_class)
            query_args.append(self.release.db.func.group_concat(agr_class_alias.name.distinct()))
            headers.append('{}-ASBs'.format(self.filter_object_keys[what_for]))

            if int(self.release.version) >= 3:
                join_tuples += [
                    (
                        agr_snp_class_alias,
                        (agr_snp_class_alias.chromosome == self.SNP.chromosome) &
                        (agr_snp_class_alias.position == self.SNP.position) &
                        (agr_snp_class_alias.alt == self.SNP.alt) &
                        (agr_snp_class_alias.fdr_class.in_(get_corresponding_fdr_classes(filters_object['fdr']))) &
                        (agr_snp_class_alias.es_class.in_(get_corresponding_es_classes(filters_object['es']))),
                    ),
                    (
                        agr_class_alias,
                        getattr(agr_class_alias, id_field) == getattr(agr_snp_class_alias, id_field),
                    )
                ]
            else:
                join_tuples += [
                    (
                        agr_snp_class_alias,
                        (agr_snp_class_alias.chromosome == self.SNP.chromosome) &
                        (agr_snp_class_alias.position == self.SNP.position) &
                        (agr_snp_class_alias.alt == self.SNP.alt),
                    ),
                    (
                        agr_class_alias,
                        getattr(agr_class_alias, id_field) == getattr(agr_snp_class_alias, id_field),
                    )
                ]

        for what_for in self.what_for:
            filter_object_key = self.filter_object_keys[what_for]
            aggregation_class = self.aggregation_classes[what_for]
            aggregated_snp_class = self.aggregation_snp_classes[what_for]
            id_field = 'cl_id'
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
                             'cl_id'
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
                        additional_columns.append(
                            self.release.db.func.coalesce(getattr(aliases[name], field)).label(label))

        Target = aliased(self.Gene)
        join_tuples.append((Target, self.SNP.target_genes))
        additional_columns.append(self.release.db.func.group_concat(Target.gene_name.distinct()).label('targets'))
        headers.append('Targets')

        found_snps = self.release.session.query(*query_args)
        found_snps = found_snps.filter(*self.construct_advanced_filters(filters_object))
        for cls, condition in join_tuples:
            found_snps = found_snps.join(cls, condition, isouter=True)
        found_snps = found_snps.add_columns(*additional_columns)
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
        cls = self.aggregation_classes[what_for]
        if int(self.release.version) >= 3:
            filters = (((cls.name.like(in_str),) if in_str else ()) +
                       ((not_(cls.name.in_(used_options)),) if used_options else ()) +
                       (cls.aggregated_snps_count010,))
            return cls.query.filter(*filters).order_by(cls.aggregated_snps_count010.desc()).limit(3).all()
        else:
            filters = (((cls.name.like(in_str),) if in_str else ()) +
                       ((not_(cls.name.in_(used_options)),) if used_options else ()) +
                       (cls.aggregated_snps_count,))
            return cls.query.filter(*filters).order_by(cls.aggregated_snps_count.desc()).limit(3).all()

    def get_gene_locus(self, gene, offset=5000):
        if int(self.release.version) >= 3:
            if gene.orientation:
                return max(gene.start_pos - offset, 1), gene.end_pos
            else:
                return gene.start_pos, gene.end_pos + offset
        else:
            return gene.start_pos, gene.end_pos

    def get_hints_for_gene_name(self, in_str):
        filters = (self.Gene.gene_name.like(in_str), self.Gene.gene_name != self.Gene.gene_id) if in_str else (
            self.Gene.gene_name != self.Gene.gene_id,)
        genes = self.Gene.query.filter(*filters).order_by(self.Gene.snps_count.desc()).limit(3).all()
        for g in genes:
            g.locus_start, g.locus_end = self.get_gene_locus(g)
        return genes

    def get_hints_for_eqtl_gene_name(self, in_str):
        filters = (self.Gene.gene_name.like(in_str),) if in_str else ()
        genes = self.Gene.query.filter(*filters).order_by(self.Gene.eqtl_snps_count.desc()).limit(3).all()
        for g in genes:
            g.locus_start, g.locus_end = self.get_gene_locus(g)
        return genes

    def get_overall_statistics(self):
        if int(self.release.version) >= 3:
            if int(self.release.version) == 3:
                stats_dict_to_use = {'0.01': {'expected_tf_asbs': 100178,
                                              'expected_cl_asbs': 156005,
                                              'expected_all_asbs': 256183,
                                              'expected_tf_asbs_rs': 79713,
                                              'expected_cl_asbs_rs': 115757,
                                              'expected_all_asbs_rs': 114925},
                                     '0.05': {'expected_tf_asbs': 183756,
                                              'expected_cl_asbs': 262013,
                                              'expected_all_asbs': 445769,
                                              'expected_tf_asbs_rs': 136045,
                                              'expected_cl_asbs_rs': 180299,
                                              'expected_all_asbs_rs': 169816},
                                     '0.1': {'expected_tf_asbs': 255812,
                                             'expected_cl_asbs': 349490,
                                             'expected_all_asbs': 605302,
                                             'expected_tf_asbs_rs': 180149,
                                             'expected_cl_asbs_rs': 229260,
                                             'expected_all_asbs_rs': 204619},
                                     '0.15': {'expected_tf_asbs': 323556,
                                              'expected_cl_asbs': 427750,
                                              'expected_all_asbs': 751306,
                                              'expected_tf_asbs_rs': 219099,
                                              'expected_cl_asbs_rs': 270903,
                                              'expected_all_asbs_rs': 228154},
                                     '0.25': {'expected_tf_asbs': 468555,
                                              'expected_cl_asbs': 585750,
                                              'expected_all_asbs': 1054305,
                                              'expected_tf_asbs_rs': 297146,
                                              'expected_cl_asbs_rs': 349733,
                                              'expected_all_asbs_rs': 253630},
                                     }
            else:
                stats_dict_to_use = stats_dict
            return {
                'transcription_factors_count': self.TranscriptionFactor.query.filter(
                    self.TranscriptionFactor.aggregated_snps_count > 0).count(),
                'cell_types_count': self.CellLine.query.filter(self.CellLine.aggregated_snps_count > 0).count(),
                'snps_count': stats_dict_to_use['0.25']['expected_all_asbs_rs'],
                'asbs_count': stats_dict_to_use['0.25']['expected_all_asbs'],
                'snps_count010': stats_dict_to_use['0.1']['expected_all_asbs_rs'],
                'asbs_count010': stats_dict_to_use['0.1']['expected_all_asbs'],
                'snps_count005': stats_dict_to_use['0.05']['expected_all_asbs_rs'],
                'asbs_count005': stats_dict_to_use['0.05']['expected_all_asbs'],
            }
        else:
            return {
                'transcription_factors_count': self.TranscriptionFactor.query.filter(
                    self.TranscriptionFactor.aggregated_snps_count > 0).count(),
                'cell_types_count': self.CellLine.query.filter(self.CellLine.aggregated_snps_count > 0).count(),
                'snps_count': self.release.session.query(self.SNP.rs_id).distinct().count(),
            }
