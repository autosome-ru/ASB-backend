from ASB_app import *
import os
import json
import numpy as np

from ASB_app.constants import chromosomes
from ASB_app.models import CandidateSNP, CandidateRS, CandidateTFRS, CandidateCLRS, PositionHash, LDIslandsInfo
from ASB_app.utils.statistics import get_fdr_class, get_es_class

current_release = releases.ReleaseBillCipher
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


TF = 0
CL = 0
SNP_RS = 0
TF_SNP = 0
CL_SNP = 0
HASH = 0
LD = 0

release_path = os.path.expanduser('~/adastra/DataChipBillCipher0501/')
parameters_path = os.path.expanduser('~/Configs/')

conv_bad = dict(zip(
    (1, 4 / 3, 3 / 2, 2, 5 / 2, 3, 4, 5, 6),
    ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6')
))

if __name__ == '__main__':
    with open(os.path.join(release_path, 'release_stats', 'convert_cell_lines.json')) as file:
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

                    float_fields = ['fdrp_bh_ref', 'fdrp_bh_alt', 'es_mean_ref', 'es_mean_alt']
                    for field in float_fields:
                        if row[field] == '' or row[field] == '.':
                            row[field] = None
                        else:
                            row[field] = float(row[field])

                    row['ID'] = int(row['ID'][row['ID'].rfind('rs') + 2:])

                    min_pv = -np.log10(min(row['fdrp_bh_ref'], row['fdrp_bh_alt'])) \
                        if (row['fdrp_bh_ref'] != 0 and row['fdrp_bh_alt'] != 0) else 310
                    max_es = max([x for x in (row['es_mean_ref'],
                                             row['es_mean_alt']) if x is not None], default=None)

                    mutation = CandidateSNP(
                        rs_id=row['ID'],
                        chromosome=row['#chr'],
                        position=row['pos'],
                        ref=row['ref'],
                        alt=row['alt'],
                        ag_level=param,
                        ag_id=ag_id,
                        best_p_value=min_pv,
                        fdr_class=get_fdr_class(min_pv),
                        best_es=max_es,
                        es_class=get_es_class(max_es),
                    )

                    snps.append(mutation)

        session.add_all(snps)
        session.commit()

        session.close()

    if SNP_RS:
        q = session.query(CandidateSNP.rs_id, db.func.max(CandidateSNP.best_p_value), db.func.max(CandidateSNP.best_es)).group_by(CandidateSNP.rs_id)
        snps = []
        for i, (rs_id, fdr, es) in enumerate(q, 1):
            if i == 1:
                print('start')
            if i % 100000 == 0:
                print(i)
            snp = CandidateRS(
                rs_id=rs_id,
                best_p_value=fdr,
                fdr_class=get_fdr_class(fdr),
                best_es=es,
                es_class=get_es_class(es),
            )
            snps.append(snp)
        session.add_all(snps)
        session.commit()

    for param in ['TF'] * TF_SNP + ['CL'] * CL_SNP:
        q = session.query(CandidateSNP.rs_id, db.func.max(CandidateSNP.best_p_value), db.func.max(CandidateSNP.best_es)).filter(CandidateSNP.ag_level==param).group_by(CandidateSNP.rs_id)
        snps = []
        SNPClass = {'TF': CandidateTFRS, 'CL': CandidateCLRS}[param]
        for i, (rs_id, fdr, es) in enumerate(q, 1):
            if i == 1:
                print('start {}'.format(param))
            if i % 100000 == 0:
                print(i)
            snp = SNPClass(
                rs_id=rs_id,
                best_p_value=fdr,
                fdr_class=get_fdr_class(fdr),
                best_es=es,
                es_class=get_es_class(es),
            )
            snps.append(snp)
        session.add_all(snps)
        session.commit()

    if HASH:
        phs = []
        rs_ids = set()
        for cand in CandidateSNP.query:
            if cand.rs_id in rs_ids:
                continue
            else:
                phs.append(
                    PositionHash(
                        rs_id=cand.rs_id,
                        position_hash=chromosomes.index(cand.chromosome) * 10 ** 9 + cand.position
                    )
                )
                rs_ids.add(cand.rs_id)
        session.add_all(phs)
        session.commit()
        session.close()

    if LD:
        lds = []
        for cand in CandidateRS.query:
            lds.append(
                LDIslandsInfo(
                    rs_id=cand.rs_id,
                )
            )
        session.add_all(lds)
        session.commit()
        session.close()

        print('created')

        adastra_rs = set(*zip(*session.query(CandidateRS.rs_id)))
        for file, column in zip(('asn_rs.tsv', 'afr_rs.tsv', 'eur_rs.tsv'), ('ld_asn', 'ld_afr', 'ld_eur')):
            print(column)
            island_num = {}
            num = 0
            num_rs = {}
            with open(os.path.join(os.path.expanduser('~/PARAMETERS/ananastra_LD'), file)) as f:
                for line in f:
                    line = line.strip('\n').split('\t')
                    rs = int(line[0][2:])
                    if rs not in adastra_rs:
                        continue
                    island = (line[1], line[2], line[3])
                    if island not in island_num:
                        print(island)
                        island_num[island] = num
                        num += 1
                    num_rs.setdefault(num, set()).add(rs)
            for num, rs_set in num_rs.items():
                print(num)
                session.execute(db.update(LDIslandsInfo).where(LDIslandsInfo.rs_id.in_(rs_set)).values(**{column: num}))
            session.commit()
            session.close()
