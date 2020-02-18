from ASB_app import *
from ASB_app.models import *
import gzip
import os
import json

TF = True
CL = False
tr = 0.2
TF_DICT = False
CL_DICT = False

release_path = os.path.expanduser('~/Releases/latest/')

conv_bad = dict(zip(
                    (1, 4/3, 3/2, 2, 5/2, 3, 4, 5, 6),
                    ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')
                ))

with open(os.path.expanduser('~/Releases/CONVERT_CL_NAMES.json')) as file:
    cl_dict = json.loads(file.readline())

cl_dict_reverse = {}
for key, value in cl_dict.items():
    cl_dict_reverse[value] = key

for param in ['TF' * TF, 'CL' * CL]:
    pv_path = release_path + '{}_P-values/'.format(param)
    for file in sorted(os.listdir(pv_path)):
        with gzip.open(pv_path + file, 'rt') as table:
            name = file.replace('.tsv.gz', '')
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

                for field in 'fdrp_bh_ref', 'fdrp_bh_alt', 'es_mean_ref', 'es_mean_alt', 'mean_BAD':
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

                for field in ['pos']:
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
                    'ref': row['ref'],
                    'alt': row['alt'],
                    ({'TF': 'tf_id', 'CL': 'cl_id'}[param]): ag_id,
                    'p_value_ref': row['fdrp_bh_ref'],
                    'p_value_alt': row['fdrp_bh_alt'],
                    'es_ref': row['es_mean_ref'],
                    'es_alt': row['es_mean_alt'],
                    'is_asb': min_pv <= 0.05
                }

                ag_snps.append(SNPClass(*ag_data))

        session.add_all(snps)
        session.commit()

        session.add_all(ag_snps)
        session.commit()
