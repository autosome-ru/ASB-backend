from ASB_app import *
from ASB_app.models import *
import gzip
import os
import json

TF = True
CL = False
tr = 0.05  # 0.2
TF_DICT = False
CL_DICT = False

conv_bad = dict(zip(
                    (1, 4/3, 3/2, 2, 5/2, 3, 4, 5, 6),
                    ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')
                ))

with open(os.path.expanduser('~/Releases/CONVERT_CL_NAMES.json')) as file:
    cl_dict = json.loads(file.readline())

cl_dict_reverse = {}
for key, value in cl_dict.items():
    cl_dict_reverse[value] = key

if TF:
    for file in sorted(os.listdir(os.path.expanduser('~/Releases/TF'))):
        with gzip.open(os.path.expanduser('~/Releases/TF/' + file), 'rt') as table:

            name = file.replace('.tsv.gz', '')
            print(name)
            tf = TranscriptionFactor.query.filter(TranscriptionFactor.name == name).first()
            if not tf:
                tf = TranscriptionFactor(name=name)
                session.add(tf)
                session.commit()
            tf_id = tf.tf_id

            tf_snps = []
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

                tf_snps.append(TranscriptionFactorSNP(
                    chromosome=row['#chr'],
                    position=int(row['pos']),
                    ref=row['ref'],
                    alt=row['alt'],
                    tf_id=tf_id,
                    p_value_ref=row['fdrp_bh_ref'],
                    p_value_alt=row['fdrp_bh_alt'],
                    es_ref=row['es_mean_ref'],
                    es_alt=row['es_mean_alt'],
                    is_asb=min_pv <= 0.05
                ))

        session.add_all(snps)
        session.commit()

        session.add_all(tf_snps)
        session.commit()

if CL:
    for file in sorted(os.listdir(os.path.expanduser('~/Releases/CL'))):
        with gzip.open(os.path.expanduser('~/Releases/CL/' + file), 'rt') as table:

            name = cl_dict_reverse[file.replace('.tsv.gz', '')]
            print(name)
            cl = CellLine.query.filter(CellLine.name == name).first()
            if not cl:
                cl = CellLine(name=name)
                session.add(cl)
                session.commit()
            cl_id = cl.cl_id

            cl_snps = []
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

                cl_snps.append(CellLineSNP(
                    chromosome=row['#chr'],
                    position=int(row['pos']),
                    ref=row['ref'],
                    alt=row['alt'],
                    cl_id=cl_id,
                    p_value_ref=row['fdrp_bh_ref'],
                    p_value_alt=row['fdrp_bh_alt'],
                    es_ref=row['es_mean_ref'],
                    es_alt=row['es_mean_alt'],
                    is_asb=min_pv <= 0.05
                ))

        session.add_all(snps)
        session.commit()

        session.add_all(cl_snps)
        session.commit()
