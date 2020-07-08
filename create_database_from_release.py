from ASB_app import *
from ASB_app.models import *
import os
import json
import numpy as np
import pandas as pd

TF = 0
CL = 0
tr = 0.05
EXP = 0
UNIPROT = 0
TF_DICT = 0
CL_DICT = 0
PHEN = 0
CONTEXT = 1

release_path = os.path.expanduser('~/RESULTS/release-220620_Waddles/')
parameters_path = os.path.expanduser('~/PARAMETERS/')

conv_bad = dict(zip(
    (1, 4 / 3, 3 / 2, 2, 5 / 2, 3, 4, 5, 6),
    ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')
))

if __name__ == '__main__':
    with open(parameters_path + 'CONVERT_CL_NAMES.json') as file:
        cl_dict = json.loads(file.readline())

    cl_dict_reverse = {}
    for key, value in cl_dict.items():
        cl_dict_reverse[value] = key

    if EXP:
        table = pd.read_table(parameters_path + 'Master-lines.tsv')
        counter = 1
        exps = []
        tfs = []
        cls = []
        used_tf_names = {}
        used_cl_ids = set()
        for index, row in table.iterrows():
            if (index + 1) % 1000 == 0:
                print(index + 1)

            if row['tf_uniprot_ac'] not in used_tf_names:
                tfs.append(TranscriptionFactor(tf_id=counter, name=row['tf_uniprot_ac']))
                used_tf_names[row['tf_uniprot_ac']] = counter
                counter += 1
            if row['cell_id'] not in used_cl_ids:
                cls.append(CellLine(cl_id=int(row['cell_id']), name=row['cell_title']))
                used_cl_ids.add(row['cell_id'])

            exp = Experiment(exp_id=int(row['#id'][3:]),
                             align=int(row['ALIGNS'][6:]),
                             geo_gse=row['GEO_GSE'] if row['GEO_GSE'] != 'None' else None,
                             encode=row['ENCODE'] if row['ENCODE'] != 'None' else None,
                             tf_id=used_tf_names[row['tf_uniprot_ac']],
                             cl_id=int(row['cell_id']))

            exps.append(exp)

        session.add_all(tfs + cls + exps)
        session.commit()
        session.close()

    if UNIPROT:
        table = pd.read_table(parameters_path + 'slice(GTRD).csv')
        counter = 1
        tfs = []
        used_tf_names = {}
        for index, row in table.iterrows():
            if (index + 1) % 1000 == 0:
                print(index + 1)

            if row['tf_uniprot_ac'] not in used_tf_names:
                tf = TranscriptionFactor.query.filter(TranscriptionFactor.name == row['tf_uniprot_ac']).first()
                if not tf:
                    continue
                tf.uniprot_ac = row['tf_uniprot_id']
                tfs.append(tf)
                used_tf_names[row['tf_uniprot_ac']] = counter
                counter += 1

        session.add_all(tfs)
        session.commit()
        session.close()

    for param in ['TF'] * TF + ['CL'] * CL:
        pv_path = release_path + '{}_P-values/'.format(param)
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
                    int_fields = ['pos']
                    if param == "TF":
                        float_fields += ['motif_log_pref', 'motif_log_palt', 'motif_fc']
                        int_fields += ['motif_pos']
                        row['motif_orient'] = {'+': True, '-': False, '': None}[row['motif_orient']]
                        row['motif_conc'] = {'concordant': True,
                                             'discordant': False,
                                             '': None,
                                             'None': None}[row['motif_conc']]

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

                    row['ID'] = int(row['ID'][2:])
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
        table = pd.read_table(release_path + '00_QTL_TFCL_fdrp_bh_0.05snpphtfASB_220620_Waddles.tsv')
        for index, row in table.iterrows():
            if (index + 1) % 1000 == 0:
                print(index + 1)
            mutations = SNP.query.filter(SNP.rs_id == int(row['RSID'][2:])).all()
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

            name = file.replace('_DICT.json', '')
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
                    rs_id = int(rs_id[2:])
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
                        exp_id = Experiment.query.filter(Experiment.align == int(parameter['aligns'][6:])).one().exp_id

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
        for file in os.path.expanduser('~/SARUS_ANNOTATION/'):
            with open(os.path.expanduser('~/SARUS_ANNOTATION/') + file):
                line = file.readline()
                while line:
                    if line.starts_with('>') and line[-3:] == 'ref' and line not in used:
                        used.add(line)
                        line = line.strip()
                        alt = line.split(';')[-1].split('_')[0]
                        rs = int(line.split(';')[0][3:])
                        snp = SNP.query.filter(SNP.rs_id == rs, SNP.alt == alt).one_or_none()
                        context = file.readline()
                        if snp:
                            snp.context = context
                    line = file.readline()
        session.commit()
