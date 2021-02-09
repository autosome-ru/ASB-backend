import math

from ASB_app import *
from ASB_app import constants
from ASB_app.models import *
import os
import json
import numpy as np
import pandas as pd

from sqlalchemy.sql import case

current_release = releases.ReleaseFord
session = current_release.session

TranscriptionFactor, \
CellLine, \
Experiment, \
ExpSNP, \
SNP, \
TranscriptionFactorSNP, \
CellLineSNP, \
Phenotype, \
PhenotypeSNPCorrespondence, \
BADGroup, \
GeneSNPCorrespondence, \
Gene = \
current_release.TranscriptionFactor, \
current_release.CellLine, \
current_release.Experiment, \
current_release.ExpSNP, \
current_release.SNP, \
current_release.TranscriptionFactorSNP, \
current_release.CellLineSNP, \
current_release.Phenotype, \
current_release.PhenotypeSNPCorrespondence, \
current_release.BADGroup, \
current_release.GeneSNPCorrespondence, \
current_release.Gene


TF = 0
CL = 0
tr = 0.25
EXP = 0
TF_DICT = 0
CL_DICT = 1
PHEN = 0
CONTEXT = 0
CONTROLS = 0
BAD_GROUP = 0
GENES = 0
TARGET_GENES = 0
PROMOTER_GENES = 0


release_path = os.path.expanduser('~/RESULTS/DataChIP/')
parameters_path = os.path.expanduser('~/RESULTS/Configs/')

conv_bad = dict(zip(
    (1, 4 / 3, 3 / 2, 2, 5 / 2, 3, 4, 5, 6),
    ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')
))

if __name__ == '__main__':
    with open(os.path.join(release_path, 'release_stats', 'convert_cl_names.json')) as file:
        cl_dict = json.loads(file.readline())

    cl_dict_reverse = {}
    for key, value in cl_dict.items():
        cl_dict_reverse[value] = key

    if EXP:
        table = pd.read_table(parameters_path + 'master-chip.txt')
        counter = 1
        exps = []
        tfs = []
        cls = []
        used_tf_names = {}
        used_cl_ids = set()
        for index, row in table.iterrows():
            if (index + 1) % 1000 == 0:
                print(index + 1)

            if row['TF_UNIPROT_NAME'] is None or pd.isna(row['TF_UNIPROT_NAME']):
                assert row['EXP_TYPE'] in ('chip_control', 'chipexo_control')
                continue

            if row['TF_UNIPROT_NAME'] not in used_tf_names:
                tfs.append(TranscriptionFactor(tf_id=counter, uniprot_ac=row['TF_UNIPROT_ID'], name=row['TF_UNIPROT_NAME']))
                used_tf_names[row['TF_UNIPROT_NAME']] = counter
                counter += 1
            if row['CELL_ID'] not in used_cl_ids:
                cls.append(CellLine(cl_id=int(row['CELL_ID']), name=row['CELLS']))
                used_cl_ids.add(row['CELL_ID'])

            exp = Experiment(exp_id=row['#EXP'],
                             align=row['ALIGNS'],
                             geo_gse=row['GEO'] if row['GEO'] != '' and not pd.isna(row['GEO']) else None,
                             encode=row['ENCODE'] if row['ENCODE'] != '' and not pd.isna(row['ENCODE']) else None,
                             tf_id=used_tf_names[row['TF_UNIPROT_NAME']],
                             cl_id=int(row['CELL_ID']))

            exps.append(exp)

        session.add_all(tfs + cls + exps)
        session.commit()
        session.close()

    for param in ['TF'] * TF + ['CL'] * CL:
        pv_path = os.path.join(release_path, '{}_P-values/'.format(param))
        for file in sorted(os.listdir(pv_path)):
            with open(pv_path + file, 'r') as table:
                name = file.replace('.tsv', '')
                if param == 'CL':
                    name = cl_dict_reverse[name]
                # elif param == 'TF':
                #     name = name.replace('_HUMAN', '')
                print(name)

                AgrClass = {'TF': TranscriptionFactor, 'CL': CellLine}[param]
                SNPClass = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[param]

                ag = AgrClass.query.filter(AgrClass.name == name).first()
                if not ag:
                    print('There is no {} {}'.format(param, name))
                    ag = AgrClass(name=name)
                    session.add(ag)
                    session.commit()
                if param == 'CL':
                    ag_id = ag.cl_id
                else:
                    ag_id = ag.tf_id

                ag_snps = []
                snps = []
                header = []
                for index, row in enumerate(table):
                    print(index + 1) if (index + 1) % 50000 == 0 else ...

                    if row[0] == '#':
                        header = row.strip('\n').split('\t')
                        continue
                    else:
                        row = dict(zip(header, row.strip('\n').split('\t')))
                    float_fields = ['fdrp_bh_ref', 'fdrp_bh_alt',
                                    'es_mean_ref', 'es_mean_alt', 'mean_BAD']
                    int_fields = ['pos', 'n_peak_calls', 'n_peak_callers']
                    if param == "TF":
                        float_fields += ['motif_log_pref', 'motif_log_palt', 'motif_fc']
                        int_fields += ['motif_pos']
                        row['motif_orient'] = {'+': True, '-': False, '': None}[row['motif_orient']]
                        row['motif_conc'] = None if row['motif_conc'] in ('None', '') else row['motif_conc']

                    for field in float_fields:
                        if row[field] == '' or row[field] == '.':
                            row[field] = None
                        else:
                            row[field] = float(row[field])

                    min_pv = min(
                        row['fdrp_bh_ref'] if row['fdrp_bh_ref'] else 1,
                        row['fdrp_bh_alt'] if row['fdrp_bh_alt'] else 1,
                    )

                    if min_pv > tr:
                        continue

                    for field in int_fields:
                        if row[field] == '' or row[field] == '.':
                            row[field] = None
                        else:
                            row[field] = int(row[field])

                    row['ID'] = int(row['ID'][row['ID'].rfind('rs') + 2:])
                    mutation = SNP.query.filter((SNP.rs_id == row['ID']) &
                                                (SNP.alt == row['alt'])).first()
                    if not mutation:
                        mutation = SNP(
                            rs_id=row['ID'],
                            chromosome=row['#chr'],
                            position=row['pos'],
                            ref=row['ref'],
                            alt=row['alt'],
                        )

                    snps.append(mutation)

                    ag_data = {
                        'chromosome': row['#chr'],
                        'position': int(row['pos']),
                        'alt': row['alt'],
                        ({'TF': 'tf_id', 'CL': 'cl_id'}[param]): ag_id,
                        'log_p_value_ref': -np.log10(row['fdrp_bh_ref']),
                        'log_p_value_alt': -np.log10(row['fdrp_bh_alt']),
                        'es_ref': row['es_mean_ref'],
                        'es_alt': row['es_mean_alt'],
                        'is_asb': min_pv <= 0.05,
                        'mean_bad': row['mean_BAD'],
                        'peak_calls': row['n_peak_calls'],
                        'peak_callers': row['n_peak_callers'],
                    }
                    if param == 'TF':
                        ag_data.update({'motif_log_p_ref': row['motif_log_pref'],
                                        'motif_log_p_alt': row['motif_log_palt'],
                                        'motif_log_2_fc': row['motif_fc'],
                                        'motif_orientation': row['motif_orient'],
                                        'motif_position': row['motif_pos'],
                                        'motif_concordance': row['motif_conc'],
                                        })
                    ag_snps.append(SNPClass(**ag_data))

            session.add_all(snps)
            session.commit()

            session.add_all(ag_snps)
            session.commit()

            session.close()

    if PHEN:
        table = pd.read_table(os.path.join(release_path, 'release_stats', 'phenotypes_stats.tsv'))
        for index, row in table.iterrows():
            if (index + 1) % 1000 == 0:
                print(index + 1)
            mutations = SNP.query.filter(SNP.rs_id == int(row['RSID'][row['RSID'].rfind('rs') + 2:])).all()
            if not mutations:
                print('No snps for ', int(row['RSID'][2:]))
            for database in ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping', 'QTL']:
                if str(row[database]) == 'nan':
                    continue
                ph_names = row[database].strip('\n').split(';')
                for mutation in mutations:
                    mutation.phenotypes += [
                        Phenotype(**{
                            'db_name': database,
                            'phenotype_name': name
                        }) for name in ph_names
                    ]
        session.commit()

    session.close()

    for param in ['TF'] * TF_DICT + ['CL'] * CL_DICT:
        pv_path = release_path + '{}_DICTS/'.format(param)
        for file in sorted(os.listdir(pv_path)):

            name = file.replace('.json', '')
            if param == 'CL':
                name = cl_dict_reverse[name]
            print(name)

            with open(pv_path + file, 'r') as info:
                content = json.loads(info.readline())

            AgrClass = {'TF': TranscriptionFactor, 'CL': CellLine}[param]
            SNPClass = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[param]

            ag = AgrClass.query.filter(AgrClass.name == name).one()
            if param == 'CL':
                ag_id = ag.cl_id
            else:
                ag_id = ag.tf_id

            # exp_snp = ExpSNP.query.filter(
            #     getattr(ExpSNP, {'TF': 'tf_aggregated_snp', 'CL': 'cl_aggregated_snp'}[param]).has(
            #         getattr(SNPClass, {'TF': 'tf_id', 'CL': 'cl_id'}[param]) == ag_id,
            #     ),
            # ).first()
            # if exp_snp:
            #     continue

            items_length = len(content)

            items = list(content.items())

            processed = 0
            chunk_size = 100000
            while items_length - processed > 0:
                exp_snps = []
                for index, (key, value) in enumerate(items[processed: min(items_length, processed + chunk_size)], 1):
                    if index % 10000 == 0:
                        print(index)
                    chromosome, position, rs_id, ref, alt = key.strip().split('\t')[:5]
                    position = int(position)
                    rs_id = int(rs_id[rs_id.rfind('rs') + 2:])
                    ag_snp = SNPClass.query.filter(
                        SNPClass.chromosome == chromosome,
                        SNPClass.position == position,
                        SNPClass.alt == alt,
                        getattr(SNPClass, {'TF': 'tf_id', 'CL': 'cl_id'}[param]) == ag_id,
                    ).first()

                    if not ag_snp:
                        continue

                    ag_snp_id = getattr(ag_snp, {'TF': 'tf_snp_id', 'CL': 'cl_snp_id'}[param])

                    AnotherAgrClass = {'CL': TranscriptionFactor, 'TF': CellLine}[param]
                    AnotherSNPClass = {'CL': TranscriptionFactorSNP, 'TF': CellLineSNP}[param]

                    another_ag_snps = AnotherSNPClass.query.filter(
                        AnotherSNPClass.chromosome == chromosome,
                        AnotherSNPClass.position == position,
                        AnotherSNPClass.alt == alt,
                    ).all()

                    another_dict = {}
                    another_id = {'CL': 'tf_snp_id', 'TF': 'cl_snp_id'}[param]
                    another_class = {'CL': 'tf_id', 'TF': 'cl_id'}[param]
                    for snp in another_ag_snps:
                        another_dict[AnotherAgrClass.query.get(getattr(snp, another_class)).name] = getattr(snp,
                                                                                                            another_id)

                    del value['ref_ef']
                    del value['alt_ef']
                    if 'logitp_ref' in value:
                        del value['logitp_ref']
                        del value['logitp_alt']

                    parameters_list = [dict(zip(
                        value.keys(),
                        [val[i] for val in value.values()],
                    ))
                        for i in range(len(value['aligns']))]

                    for parameter in parameters_list:
                        # FIXME TEMPORARY
                        exp_id = Experiment.query.filter(Experiment.align == parameter['aligns'][0]).one().exp_id

                        exp_snp = ExpSNP.query.filter(
                            ExpSNP.exp_id == exp_id,
                            getattr(ExpSNP, {'TF': 'tf_snp_id', 'CL': 'cl_snp_id'}[param]) == ag_snp_id,
                        ).first()

                        if not exp_snp:
                            exp_snp = ExpSNP(**{
                                'ref_readcount': parameter['ref_counts'],
                                'alt_readcount': parameter['alt_counts'],
                                'p_value_ref': parameter['ref_pvalues'],
                                'p_value_alt': parameter['alt_pvalues'],
                                'bad': conv_bad[parameter['BAD']],
                                'tf_snp_id': {'TF': ag_snp_id, 'CL': another_dict.get(parameter.get('TF'))}[param],
                                'cl_snp_id': {'TF': another_dict.get(cl_dict_reverse.get(parameter.get('CL'))),
                                              'CL': ag_snp_id}[param],
                                'exp_id': exp_id,
                            })
                        else:
                            other_id = getattr(exp_snp, {'TF': 'cl_snp_id', 'CL': 'tf_snp_id'}[param])
                            assert other_id == {'TF': another_dict.get(cl_dict_reverse.get(parameter.get('CL'))),
                                                'CL': another_dict.get(parameter.get('TF'))}[param]
                            assert exp_snp.ref_readcount == parameter['ref_counts']
                            assert exp_snp.p_value_alt == parameter['alt_pvalues']
                            assert exp_snp.bad == conv_bad[parameter['BAD']]

                        exp_snps.append(exp_snp)

                session.add_all(exp_snps)
                session.commit()
                session.close()
                processed += chunk_size

    if CONTEXT:
        used = set()
        with open(os.path.join(release_path, 'Sarus', 'all_tfs.fasta')) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if line.startswith('>') and line[-3:] == 'ref' and line not in used:
                    used.add(line)
                    alt = line.split(';')[-1].split('_')[0]
                    rs = int(line.split(';')[0][3:])
                    snp = SNP.query.filter(SNP.rs_id == rs, SNP.alt == alt).one_or_none()
                    context = file.readline().strip('\n')
                    if snp:
                        snp.context = context
                line = file.readline()
        session.commit()

    if CONTROLS:
        table = pd.read_table(parameters_path + 'master-chip.txt')
        exps = []
        cls = []
        used_exp_ids = set()
        used_cl_ids = set([x[0] for x in session.query(CellLine.cl_id.distinct())])
        for index, row in table.iterrows():
            if (index + 1) % 1000 == 0:
                print(index + 1)

            if len(exps) >= 990:
                session.add_all(cls + exps)
                session.commit()
                exps = []
                cls = []
                session.close()

            if not (row['TF_UNIPROT_NAME'] is None or pd.isna(row['TF_UNIPROT_NAME'])):
                continue
            assert row['EXP_TYPE'] in ('chip_control', 'chipexo_control')

            if row['#EXP'] in used_exp_ids:
                continue
            used_exp_ids.add(row['#EXP'])
            if row['CELL_ID'] not in used_cl_ids:
                cls.append(CellLine(cl_id=int(row['CELL_ID']), name=row['CELLS']))
                used_cl_ids.add(row['CELL_ID'])

            exp = Experiment(exp_id=row['#EXP'],
                             align=row['ALIGNS'],
                             geo_gse=row['GEO'] if row['GEO'] != '' and not pd.isna(row['GEO']) else None,
                             encode=row['ENCODE'] if row['ENCODE'] != '' and not pd.isna(row['ENCODE']) else None,
                             tf_id=None,
                             cl_id=int(row['CELL_ID']),
                             is_control=True)

            exps.append(exp)

        session.add_all(cls + exps)
        session.commit()
        session.close()

    if BAD_GROUP:
        with open(os.path.join(release_path, 'release_stats', 'badmaps_dict.json')) as f:
            cell_lines_dict = json.loads(f.readline())
        exps = []
        bad_groups = []
        for key, value in cell_lines_dict.items():
            print(key)
            name = key
            bad_group = BADGroup.query.filter(BADGroup.bad_group_name == name).one_or_none()
            if not bad_group:
                bad_group = BADGroup(
                    bad_group_name=name
                )
            bad_groups.append(bad_group)
            for path in value:
                if len(exps) >= 300:
                    session.add_all(exps)
                    session.commit()
                    exps = []
                    session.close()
                exp_id = path.split('/')[-2]
                exp = Experiment.query.get(exp_id)
                if not exp:
                    continue
                exp.bad_group = bad_group
                exps.append(exp)
                print(exp)
        session.add_all(exps + bad_groups)
        session.commit()
        session.close()

    if GENES:
        genes = []
        genes_ids = set()
        with open(os.path.expanduser('~/REFERENCE/gencode.v35.annotation.gtf')) as inp:
            for index, line in enumerate(inp):
                if line.startswith('#'):
                    continue
                line = line.strip('\n').split('\t')
                chrom, start_pos, end_pos, orient = line[0], int(line[3]), int(line[4]), line[6]
                if chrom not in constants.chromosomes or line[2] != 'gene':
                    continue
                if index % 1000 == 0:
                    print(index, len(genes))
                params_dict = dict(map(lambda x: tuple(x.split(' ')), line[8].split('; ')))
                gene_name = params_dict['gene_name'].strip('"')
                gene_id = params_dict['gene_id'].strip('"')
                if orient == '+':
                    start_pos = max(start_pos - 5000, 1)
                elif orient == '-':
                    end_pos = end_pos + 5000
                else:
                    raise ValueError

                snps = SNP.query.filter(SNP.chromosome == chrom,
                                        SNP.position.between(start_pos, end_pos)).count()

                gene = Gene(gene_id=gene_id, gene_name=gene_name, start_pos=start_pos, end_pos=end_pos, chromosome=chrom,
                            orientation=True if orient == '+' else False if orient == '-' else None, snps_count=snps)
                if gene_id in genes_ids:
                    print(gene_id, chrom, start_pos, end_pos)
                    continue
                genes.append(gene)
                genes_ids.add(gene_id)

        gene_names = [g.gene_name for g in genes]
        repeating_gene_names = set()
        used_names = set()
        for name in gene_names:
            if name in used_names:
                repeating_gene_names.add(name)
            else:
                used_names.add(name)

        genes = [g for g in genes if g.gene_name not in repeating_gene_names]

        session.add_all(genes)
        session.commit()

    if TARGET_GENES:
        # table = pd.read_table(os.path.join(release_path, 'release_stats', 'phenotypes_stats.tsv'))
        table = pd.read_table(os.path.join(release_path, 'release_stats', 'phenotypes_stats.tsv'))
        print(len(table.index))
        genes = []
        for index, row in table.iterrows():
            if (index + 1) % 1000 == 0:
                print(index + 1)
            if str(row['QTLg']) in ('nan', '', 'None'):
                continue

            all_target_genes = []
            for id in row['QTLg'].strip('\n').split(';'):
                target_genes = Gene.query.filter(Gene.gene_id.like(id.split('.')[0] + '%')).all()
                if target_genes:
                    if len(set(g.gene_name for g in target_genes)) != 1:
                        print('Bad genes: {}'.format(target_genes))
                    gene = target_genes[0]
                    all_target_genes.append(gene)
                else:
                    gene = Gene(gene_id=id, gene_name=id, chromosome='chr1', start_pos=1, end_pos=1, orientation=True)
                    genes.append(gene)
                    all_target_genes.append(gene)

            mutations = SNP.query.filter(SNP.rs_id == int(row['RSID'][row['RSID'].rfind('rs') + 2:])).all()
            if not mutations:
                print('No snps for ', int(row['RSID'][2:]))

            for mutation in mutations:
                mutation.target_genes = all_target_genes
        session.add_all(genes)
        session.commit()

    if PROMOTER_GENES:
        genes = []
        for index, gene in enumerate(Gene.query.filter(~((Gene.start_pos == 1) & (Gene.end_pos == 1)))):
            if (index + 1) % 1000 == 0:
                print(index + 1)

            gene.proximal_promoter_snps = SNP.query.filter(
                    SNP.chromosome == gene.chromosome,
                    SNP.position.between(gene.start_pos - 500, gene.end_pos) if gene.orientation
                    else SNP.position.between(gene.start_pos, gene.end_pos + 500)
            ).all()
            genes.append(gene)
        session.commit()
