from ASB_app import *
from ASB_app.models import *
import os
import json
import numpy as np
import pandas as pd

TF = False
CL = False
tr = 0.05
EXP = False
TF_DICT = False
CL_DICT = False
PHEN = True

release_path = os.path.expanduser('~/Releases/latest/')

conv_bad = dict(zip(
    (1, 4 / 3, 3 / 2, 2, 5 / 2, 3, 4, 5, 6),
    ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')
))

with open(release_path + 'CONVERT_CL_NAMES.json') as file:
    cl_dict = json.loads(file.readline())

cl_dict_reverse = {}
for key, value in cl_dict.items():
    cl_dict_reverse[value] = key

for param in ['TF'] * TF + ['CL'] * CL:
    pv_path = release_path + '{}_P-values/'.format(param)
    for file in sorted(os.listdir(pv_path)):
        with open(pv_path + file, 'r') as table:
            name = file.replace('.tsv', '')
            if param == 'CL':
                name = cl_dict_reverse[name]
            elif param == 'TF':
                name = name.replace('_HUMAN', '')
            print(name)

            AgrClass = {'TF': TranscriptionFactor, 'CL': CellLine}[param]
            SNPClass = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[param]

            ag = AgrClass.query.filter(AgrClass.name == name).first()
            if not ag:
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

if PHEN:
    table = pd.read_table(release_path + '00_QTL_TFCL_fdrp_bh_0.05snpphtfASB_230120_Durland.tsv')
    for index, row in table.iterrows():
        if (index + 1) % 1000 == 0:
            print(index + 1)
        mutations = SNP.query.filter(SNP.rs_id == int(row['RSID'][2:])).all()
        # if not mutations:
        #     print('ЗАЧЕМ ТЫ УДАЛИЛ НАМ БАЗУ?', int(row['RSID'][2:]))
        # print('GOOD BOI', int(row['RSID'][2:]))
        for database in ['grasp', 'ebi', 'clinvar', 'phewas', 'finemapping', 'QTL']:
            if str(row[database]) == 'nan':
                continue
            ph_names = row[database].strip('\n').split(';')
            for mutation in mutations:
                mutation.phenotypes = [
                    Phenotype(**{
                        'db_name': database,
                        'phenotype_name': name
                    }) for name in ph_names
                ]
    session.commit()

# for param in ['TF'] * TF_DICT, ['CL'] * CL_DICT:
#     pv_path = release_path + '{}_DICTS/'.format(param)
#     for file in sorted(os.listdir(pv_path)):
#         with open(pv_path + file, 'r') as info:
#             content = json.loads(info.readline())
#             name = file.replace('.tsv', '')
#             if param == 'CL':
#                 name = cl_dict_reverse[name]
#             elif param == 'TF':
#                 name = name.replace('_HUMAN', '')
#             print(name)
#
#             AgrClass = {'TF': TranscriptionFactor, 'CL': CellLine}[param]
#             SNPClass = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[param]
#
#             ag = AgrClass.query.filter(AgrClass.name == name).first()
#             if not ag:
#                 ag = AgrClass(name=name)
#                 session.add(ag)
#                 session.commit()
#             if param == 'CL':
#                 ag_id = ag.cl_id
#             else:
#                 ag_id = ag.tf_id
