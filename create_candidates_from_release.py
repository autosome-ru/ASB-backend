from ASB_app import *
import os
import json
import numpy as np

from ASB_app.models import CandidateSNP, CandidateRS, CandidateTFRS, CandidateCLRS
from ASB_app.utils.statistics import get_fdr_class

current_release = releases.ReleaseSusan
session = current_release.session
db = current_release.db

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
current_release.Gene


TF = 1
CL = 1
SNP_RS = 1
TF_SNP = 1
CL_SNP = 1


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

    for param in ['TF'] * TF + ['CL'] * CL:
        pv_path = release_path + '{}_P-values/'.format(param)
        snps = []
        for file in sorted(os.listdir(pv_path)):
            with open(pv_path + file, 'r') as table:
                name = file.replace('.tsv', '')
                if param == 'CL':
                    name = cl_dict_reverse[name]
                print(name)

                AgrClass = {'TF': TranscriptionFactor, 'CL': CellLine}[param]
                SNPClass = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[param]

                ag = AgrClass.query.filter(AgrClass.name == name).first()
                if not ag:
                    print('There is no {} {}'.format(param, name))
                    exit(1)
                if param == 'CL':
                    ag_id = ag.cl_id
                else:
                    ag_id = ag.tf_id

                header = []
                for index, row in enumerate(table):
                    print(index + 1) if (index + 1) % 50000 == 0 else ...

                    if row[0] == '#':
                        header = row.strip('\n').split('\t')
                        continue
                    else:
                        row = dict(zip(header, row.strip('\n').split('\t')))
                    if row['fdrp_bh_ref'] == '' or row['fdrp_bh_ref'] == '.':
                        continue
                    int_fields = ['pos']
                    for field in int_fields:
                        if row[field] == '' or row[field] == '.':
                            row[field] = None
                        else:
                            row[field] = int(row[field])

                    float_fields = ['fdrp_bh_ref', 'fdrp_bh_alt']
                    for field in float_fields:
                        if row[field] == '' or row[field] == '.':
                            row[field] = None
                        else:
                            row[field] = float(row[field])

                    row['ID'] = int(row['ID'][row['ID'].rfind('rs') + 2:])

                    mutation = CandidateSNP(
                        rs_id=row['ID'],
                        chromosome=row['#chr'],
                        position=row['pos'],
                        ref=row['ref'],
                        alt=row['alt'],
                        ag_level=param,
                        ag_id=ag_id,
                        best_p_value=-np.log10(min(row['fdrp_bh_ref'], row['fdrp_bh_alt']))
                    )

                    snps.append(mutation)

        session.add_all(snps)
        session.commit()

        session.close()

    if SNP_RS:
        q = session.query(CandidateSNP.rs_id, db.func.max(CandidateSNP.best_p_value)).group_by(CandidateSNP.rs_id)
        snps = []
        for i, (rs_id, fdr) in enumerate(q, 1):
            if i == 1:
                print('start')
            if i % 100000 == 0:
                print(i)
            snp = CandidateRS(
                rs_id=rs_id,
                best_p_value=fdr,
                fdr_class=get_fdr_class(fdr)
            )
            snps.append(snp)
        session.add_all(snps)
        session.commit()

    for param in ['TF'] * TF_SNP + ['CL'] * CL_SNP:
        q = session.query(CandidateSNP.rs_id, db.func.max(CandidateSNP.best_p_value)).filter(CandidateSNP.ag_level==param).group_by(CandidateSNP.rs_id)
        snps = []
        SNPClass = {'TF': CandidateTFRS, 'CL': CandidateCLRS}[param]
        for i, (rs_id, fdr) in enumerate(q, 1):
            if i == 1:
                print('start {}'.format(param))
            if i % 100000 == 0:
                print(i)
            snp = SNPClass(
                rs_id=rs_id,
                best_p_value=fdr,
                fdr_class=get_fdr_class(fdr)
            )
            snps.append(snp)
        session.add_all(snps)
        session.commit()
