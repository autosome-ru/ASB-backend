from ASB_app import *
from ASB_app import constants
from ASB_app.models import *
from ASB_app.utils.aggregates import update_phenotype_associations, \
    update_best_p_value, update_best_es
import os
import json
import numpy as np
import pandas as pd

from tqdm import tqdm

from ASB_app.utils.statistics import get_fdr_class, get_es_class

current_release = releases.ReleaseBmo
session = current_release.session

Dnase, \
Faire, \
Atac, \
Experiment, \
ExpSNP, \
SNP, \
DnaseSNP, \
FaireSNP, \
AtacSNP, \
Phenotype, \
PhenotypeSNPCorrespondence, \
BADGroup, \
GeneSNPCorrespondence, \
Gene = \
current_release.Dnase, \
current_release.Faire, \
current_release.Atac, \
current_release.Experiment, \
current_release.ExpSNP, \
current_release.SNP, \
current_release.DnaseSNP, \
current_release.FaireSNP, \
current_release.AtacSNP, \
current_release.Phenotype, \
current_release.PhenotypeSNPCorrespondence, \
current_release.BADGroup, \
current_release.GeneSNPCorrespondence, \
current_release.Gene

tr = 0.25

EXP = 0
FAIRE = 0
DNASE = 0
ATAC = 0
PHEN = 0
FAIRE_DICT = 0
ATAC_DICT = 0
DNASE_DICT = 0

CONTEXT = 0
BAD_GROUP = 0
GENES = 0
TARGET_GENES = 0
PROMOTER_GENES = 0  # not needed at first time
TARGET_GENE_SNP_COUNT = 0

# Don't forget to change current_release in releases.py
UPDATE_PHEN_COUNT = 0
UPDATE_BEST_P_VALUE = 1
UPDATE_BEST_ES = 1
PROMOTER_GENE_COUNT = 1
TARGET_GENE_COUNT_010 = 0
PROMOTER_GENE_COUNT_010 = 1
SET_NONE_TO_ZERO = 1
CHECK_NONE = 1

# Gene name in tfs is not updated


release_path = os.path.expanduser('~/adastra/DataChipBillCipher0626')
parameters_path = os.path.expanduser('~/Configs/')

conv_bad = dict(zip(
    (1, 4 / 3, 3 / 2, 2, 5 / 2, 3, 4, 5, 6, 1.33333),
    ('1', '4/3', '3/2', '2', '5/2', '3', '4', '5', '6', '4/3')
))
agr_class_dict = {'dnase': Dnase, 'faire': Faire, 'atac': Atac}
snp_class_dict = {'dnase': DnaseSNP, 'faire': FaireSNP, 'atac': AtacSNP}
if __name__ == '__main__':
    cl_dict = {}
    cl_dict_reverse = {}
    with open('/home/abramov/udacha_v1_id2name.tsv') as file:
        for line in file:
            if not line.strip():
                continue
            try:
                key, value = line.strip().split('\t')
            except:
                print(line)
            cl_dict_reverse[value] = key
            cl_dict[key] = value

    if EXP:
        for param in ['faire'] * FAIRE + ['dnase'] * DNASE + ['atac'] * ATAC:
            print('Loading experiments')
            table = pd.read_table(f'/home/abramov/metadata_{param}.tsv')

            table2 = pd.read_table('/home/abramov/Configs/' + {'faire': 'faire+header.txt', 'atac': 'master-atac.txt', 'dnase': 'master-dnase.txt'}[param])
            print(param, len(table2.index), len(table.index), len(table.merge(table2).index))
            table = table.merge(table2)
            counter = 1
            exps = []
            cls = []
            used_tf_names = {}
            used_cl_ids = set()

            for index, row in tqdm(table.drop_duplicates(subset='#EXP').iterrows(), total=len(table.index)):
                ### CELL ID IS MISSING
                if row['CELL_ID'] not in used_cl_ids:
                    cls.append(agr_class_dict[param](cl_id=int(row['CELL_ID']), name=row['CELLS']))
                    used_cl_ids.add(row['CELL_ID'])

                exp = Experiment(exp_id=row['#EXP'],
                                 align=row['ALIGNS'],
                                 geo_gse=row['GEO'] if row['GEO'] != '' and not pd.isna(row['GEO']) else None,
                                 encode=row['ENCODE'] if row['ENCODE'] != '' and not pd.isna(row['ENCODE']) else None,
                                 **{f'{param}_id': int(row['CELL_ID'])})

                exps.append(exp)

            session.add_all(cls + exps)
            session.commit()
        session.close()

    for param in ['faire'] * FAIRE + ['dnase'] * DNASE + ['atac'] * ATAC:
        print('Loading {} ASBs'.format(param))
        pv_path = f'/home/safronov/Projects/UDACHA/release_BMO/{param}/'
        for file in tqdm(sorted(os.listdir(pv_path))):
            df = pd.read_table(pv_path + file)
            name = file.replace('.tsv', '')
            name = cl_dict_reverse[name]
            AgrClass = agr_class_dict[param]
            SNPClass = snp_class_dict[param]
            ag = AgrClass.query.filter(AgrClass.name == name).first()
            if not ag:
                # print('There is no {} {}'.format(param, name))
                ag = AgrClass(name=name)
                session.add(ag)
                session.commit()
            ag_id = ag.cl_id
            ag_snps = []
            snps = []
            header = []
            for index, row in tqdm(df.iterrows(), leave=False, total=len(df.index)):
                row['fdrp_bh_ref'] = row['ref_fdr_comb_pval']
                row['fdrp_bh_alt'] = row['alt_fdr_comb_pval']
                row['es_mean_ref'] = row['ref_comb_es']
                row['es_mean_alt'] = row['alt_comb_es']

                min_pv = min(
                    row['ref_fdr_comb_pval'] if row['ref_fdr_comb_pval'] else 1,
                    row['alt_fdr_comb_pval'] if row['alt_fdr_comb_pval'] else 1,
                )

                if min_pv > tr:
                    continue

                max_es = max(x for x in (row['ref_comb_es'],
                                         row['alt_comb_es']) if x is not None)

                row['ID'] = int(row['id'][row['id'].rfind('rs') + 2:])
                mutation = SNP.query.filter((SNP.rs_id == row['ID']) &
                                            (SNP.alt == row['alt'])).first()
                if not mutation:
                    mutation = SNP(
                        rs_id=row['ID'],
                        chromosome=row['#chr'],
                        position=row['end'],
                        ref=row['ref'],
                        alt=row['alt'],
                    )

                snps.append(mutation)

                ag_data = {
                    'chromosome': row['#chr'],
                    'position': int(row['end']),
                    'alt': row['alt'],
                    'cl_id': ag_id,
                    'log_p_value_ref': -np.log10(row['fdrp_bh_ref']) if row['fdrp_bh_ref'] != 0 else 310,
                    'log_p_value_alt': -np.log10(row['fdrp_bh_alt']) if row['fdrp_bh_alt'] != 0 else 310,
                    'best_p_value': -np.log10(min_pv),
                    'best_es': max_es,
                    'fdr_class': get_fdr_class(-np.log10(min_pv)),
                    'es_class': get_es_class(max_es),
                    'es_ref': row['es_mean_ref'],
                    'es_alt': row['es_mean_alt'],
                    'is_asb': min_pv <= 0.1,
                    'mean_bad': row['mean_bad'],
                }

                ag_snps.append(SNPClass(**ag_data))

            session.add_all(snps)
            session.commit()

            session.add_all(ag_snps)
            session.commit()

            session.close()

    if PHEN:
        print('Loading phenotypes')
        table = pd.read_table(os.path.join(release_path, 'release_stats', 'phenotypes_stats.tsv'))
        for index, row in tqdm(table.iterrows(), total=len(table.index)):
            mutations = SNP.query.filter(SNP.rs_id == int(row['RSID'][row['RSID'].rfind('rs') + 2:])).all()
            # if not mutations:
                # print('No snps for ', int(row['RSID'][2:]))
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

    for param in ['faire'] * FAIRE_DICT + ['dnase'] * DNASE_DICT + ['atac'] * ATAC_DICT:
        print('Loading {} experiment snps'.format(param))
        pv_path = f'/home/safronov/Projects/UDACHA/release_BMO/raw_{param}/'
        for file in tqdm(sorted(os.listdir(pv_path))):
            if file.startswith('.'):
                continue
            name = cl_dict_reverse[file]
            wd = f'{pv_path}/{file}/'
            AgrClass = agr_class_dict[param]
            SNPClass = snp_class_dict[param]
            cdf_list = []
            for f in os.listdir(wd):
                a = pd.read_table(f'{wd}/{f}')
                a['aligns'] = f.split('.')[0]
                cdf_list.append(a)
            df_l = pd.concat(cdf_list)
            content = {}
            # Ugly hack
            vals = ['aligns', 'ref_count', 'alt_count', 'ref_pval', 'alt_pval', 'bad']
            for index, row in df_l.iterrows():
                key = '\t'.join(map(str, [row['#chr'], row['end'], row['id'], row['ref'], row['alt']]))
                for x in vals:
                    content.setdefault(key, {}).setdefault(x, []).append(row[x])
            ag = AgrClass.query.filter(AgrClass.name == name).one()
            ag_id = ag.cl_id

            items_length = len(content)
            total_length = items_length

            items = list(content.items())

            processed = 0
            chunk_size = 100000
            while items_length - processed > 0:
                exp_snps = []
                for key, value in tqdm(items[processed: min(items_length, processed + chunk_size)],
                                       leave=False,
                                       initial=processed,
                                       total=total_length):
                    chromosome, position, rs_id, ref, alt = key.strip().split('\t')[:5]
                    position = int(position)
                    rs_id = int(rs_id[rs_id.rfind('rs') + 2:])
                    ag_snp = SNPClass.query.filter(
                        SNPClass.chromosome == chromosome,
                        SNPClass.position == position,
                        SNPClass.alt == alt,
                        getattr(SNPClass, 'cl_id') == ag_id,
                    ).first()

                    if not ag_snp:
                        continue

                    ag_snp_id = ag_snp.cl_snp_id

                    parameters_list = [dict(zip(
                        value.keys(),
                        [val[i] for val in value.values()],
                    ))
                        for i in range(len(value['aligns']))]

                    for parameter in parameters_list:
                        exp_id = Experiment.query.filter(Experiment.align == parameter['aligns']).one().exp_id

                        exp_snp = ExpSNP.query.filter(
                            ExpSNP.exp_id == exp_id,
                            getattr(ExpSNP, f'{param}_snp_id') == ag_snp_id,
                        ).first()

                        if not exp_snp:
                            exp_snp = ExpSNP(**{
                                'ref_readcount': parameter['ref_count'],
                                'alt_readcount': parameter['alt_count'],
                                'p_value_ref': parameter['ref_pval'],
                                'p_value_alt': parameter['alt_pval'],
                                'bad': conv_bad[parameter['bad']],
                                f'{param}_snp_id': ag_snp_id,
                                'exp_id': exp_id,
                            })
                        else:
                            assert exp_snp.ref_readcount == parameter['ref_count']
                            try:
                                assert round(exp_snp.p_value_alt, 3) == round(parameter['alt_pval'], 3)
                            except AssertionError:
                                print(exp_snp.p_value_alt, parameter['alt_pval'])
                            assert exp_snp.bad == conv_bad[parameter['bad']]

                        exp_snps.append(exp_snp)

                session.add_all(exp_snps)
                session.commit()
                session.close()
                processed += chunk_size

    # if CONTEXT:
    #     print('Loading SNP context')
    #     used = set()
    #     with open(os.path.join(release_path, 'Sarus', 'all_tfs.fasta')) as file:
    #         line = file.readline()
    #         while line:
    #             line = line.strip('\n')
    #             if line.startswith('>') and line not in used:
    #                 # '> rs000000@A@ref'
    #                 used.add(line)
    #                 rs, alt, allele = line[2:].strip('\n').split('@')
    #                 context = file.readline().strip('\n')
    #                 if allele == 'alt':
    #                     # '> rs000000@G@alt'
    #                     rs = int(rs[2:])
    #                     snp = SNP.query.filter(SNP.rs_id == rs, SNP.alt == alt).one_or_none()
    #                     if snp:
    #                         snp.context = context
    #             line = file.readline()
    #     session.commit()

    if BAD_GROUP:
        for param in ['faire'] * FAIRE + ['dnase'] * DNASE + ['atac'] * ATAC:
            print('Loading BAD groups')
            with open(f'/home/safronov/Projects/UDACHA/meta_info/badmaps_dicts/badmaps_dict_{param}.json') as f:
                cell_lines_dict = json.load(f)
            exps = []
            bad_groups = []
            for key, value in tqdm(cell_lines_dict.items()):
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
            session.add_all(exps + bad_groups)
            session.commit()
            session.close()

    if GENES:
        print('Loading genes')
        genes = []
        genes_ids = set()
        with open(os.path.expanduser('~/REFERENCE/gencode.v35.annotation.gtf')) as inp:
            num_items = len(inp.readlines())
        with open(os.path.expanduser('~/REFERENCE/gencode.v35.annotation.gtf')) as inp:
            for line in tqdm(inp, total=num_items):
                if line.startswith('#'):
                    continue
                line = line.strip('\n').split('\t')
                chrom, start_pos, end_pos, orient = line[0], int(line[3]), int(line[4]), line[6]
                if chrom not in constants.chromosomes or line[2] != 'gene':
                    continue
                params_dict = dict(map(lambda x: tuple(x.split(' ')), line[8].split('; ')))
                gene_name = params_dict['gene_name'].strip('"')
                gene_id = params_dict['gene_id'].strip('"')
                if orient == '+':
                    start_pos_ext = max(start_pos - 5000, 1)
                    end_pos_ext = end_pos
                elif orient == '-':
                    start_pos_ext = start_pos
                    end_pos_ext = end_pos + 5000
                else:
                    raise ValueError

                snps = SNP.query.filter(SNP.chromosome == chrom,
                                        SNP.position.between(start_pos_ext, end_pos_ext)).all()

                gene = Gene(gene_id=gene_id, gene_name=gene_name, start_pos=start_pos, end_pos=end_pos, chromosome=chrom,
                            orientation=True if orient == '+' else False if orient == '-' else None, snps_count=len(snps),
                            proximal_promoter_snps=snps)
                if gene_id in genes_ids:
                    # print(gene_id, chrom, start_pos, end_pos)
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
        print('Loading target genes')
        # table = pd.read_table(os.path.join(release_path, 'release_stats', 'phenotypes_stats.tsv'))
        table = pd.read_table(os.path.join(release_path, 'release_stats', 'phenotypes_stats.tsv'))
        genes = []
        for index, row in tqdm(table.iterrows(), total=len(table.index)):
            if str(row['QTLg']) in ('nan', '', 'None'):
                continue

            all_target_genes = []
            for id in row['QTLg'].strip('\n').split(';'):
                target_genes = Gene.query.filter(Gene.gene_id.like(id.split('.')[0] + '%')).all()
                if target_genes:
                    # if len(set(g.gene_name for g in target_genes)) != 1:
                    #     print('Bad genes: {}'.format(target_genes))
                    gene = target_genes[0]
                    all_target_genes.append(gene)
                else:
                    gene = Gene(gene_id=id, gene_name=id, chromosome='chr1', start_pos=1, end_pos=1, orientation=True)
                    genes.append(gene)
                    all_target_genes.append(gene)

            mutations = SNP.query.filter(SNP.rs_id == int(row['RSID'][row['RSID'].rfind('rs') + 2:])).all()
            # if not mutations:
            #     print('No snps for ', int(row['RSID'][2:]))

            for mutation in mutations:
                mutation.target_genes = all_target_genes
        session.add_all(genes)
        session.commit()

    if PROMOTER_GENES:
        print('Updating promoter snps')
        genes = []
        q = Gene.query.filter(~((Gene.start_pos == 1) & (Gene.end_pos == 1)))
        for gene in tqdm(Gene.query.filter(~((Gene.start_pos == 1) & (Gene.end_pos == 1))), total=q.count()):

            gene.proximal_promoter_snps = SNP.query.filter(
                    SNP.chromosome == gene.chromosome,
                    SNP.position.between(gene.start_pos - 5000, gene.end_pos) if gene.orientation
                    else SNP.position.between(gene.start_pos, gene.end_pos + 5000)
            ).all()
            genes.append(gene)
        session.commit()

    if TARGET_GENE_SNP_COUNT:
        print('Updating target snp count')
        q = session.query(Gene, db.func.count('*')).join(SNP, Gene.snps_by_target).group_by(Gene)
        for gene, count in tqdm(
            session.query(Gene, db.func.count('*')).join(SNP, Gene.snps_by_target).group_by(Gene),
            total=q.count()
        ):
            gene.eqtl_snps_count = count
        session.commit()
        for gene in Gene.query.filter(Gene.eqtl_snps_count.is_(None)):
            gene.eqtl_snps_count = 0
        session.commit()
        session.close()

    # if REDO_CONCORDANCE:
    #     print('Rereading motif columns')
    #     def to_type(val, typ):
    #         if val == '' or val == '.' or pd.isna(val):
    #             return None
    #         else:
    #             return typ(val)
    #     float_field = ['motif_log_pref', 'motif_log_palt', 'motif_fc']
    #     int_field = ['motif_pos']
    #     for tf in tqdm(TranscriptionFactor.query.all(), position=0):
    #         edited_snps = []
    #         path = os.path.join(release_path, 'TF_P-values', tf.name + '.tsv')
    #         if not os.path.exists(path):
    #             continue
    #         tf_pval_df = pd.read_table(path)
    #         tf_pval_df['key'] = tf_pval_df.apply(lambda x:
    #                                              '@'.join(map(str,
    #                                                           [x['#chr'],
    #                                                            x['pos'],
    #                                                            x['alt']])),
    #                                              axis=1)
    #         tf_pval_df = tf_pval_df.set_index('key')
    #         for tf_snp, snp in tqdm(session.query(
    #                 TranscriptionFactorSNP, SNP
    #         ).join(
    #             SNP,
    #             TranscriptionFactorSNP.snp
    #         ).filter(TranscriptionFactorSNP.tf_id == tf.tf_id).all(), position=1, leave=False):
    #             key = '@'.join(map(str, [snp.chromosome, snp.position, snp.alt]))
    #             snp_df = tf_pval_df.loc[key]
    #             tf_snp.motif_log_p_ref = to_type(snp_df['motif_log_pref'], float)
    #             tf_snp.motif_log_p_alt = to_type(snp_df['motif_log_palt'], float)
    #             tf_snp.motif_log_2_fc = to_type(snp_df['motif_fc'], float)
    #             tf_snp.motif_position = to_type(snp_df['motif_pos'], int)
    #             tf_snp.motif_orientation = {'+': True, '-': False}.get(snp_df['motif_orient'])
    #             conc = snp_df['motif_conc']
    #             tf_snp.motif_concordance = None if conc in ('None', '') or pd.isna(conc) else conc
    #             edited_snps.append(tf_snp)
    #         session.commit()
    #

    # if UPDATE_CONCORDANCE:
    #     print('Updating motif concordance')
    #     update_motif_concordance()

    if UPDATE_PHEN_COUNT:
        print('Updating phenotype associations counts')
        update_phenotype_associations()

    if UPDATE_BEST_P_VALUE:
        print('Updating best p-value')
        update_best_p_value()

    if UPDATE_BEST_ES:
        print('Updating best effect_size')
        update_best_es()

    if TARGET_GENE_COUNT_010:
        print('Updating target snp count 010')
        q = session.query(Gene, db.func.count('*')).join(SNP, Gene.snps_by_target).filter(SNP.fdr_class.in_(['0.01', '0.05', '0.1'])).group_by(Gene)
        for gene, count in tqdm(
            session.query(Gene, db.func.count('*')).join(SNP, Gene.snps_by_target).filter(
                SNP.fdr_class.in_(['0.01', '0.05', '0.1'])).group_by(Gene),
            total=q.count()
        ):
            gene.eqtl_snps_count010 = count
        session.commit()
        for gene in Gene.query.filter(Gene.eqtl_snps_count010.is_(None)):
            gene.eqtl_snps_count010 = 0
        session.commit()
        session.close()

    if PROMOTER_GENE_COUNT:
        print('Updating promoter snp count')
        q = session.query(Gene, db.func.count('*')).join(SNP, Gene.proximal_promoter_snps).group_by(Gene)
        for gene, count in tqdm(
            session.query(Gene, db.func.count('*')).join(SNP, Gene.proximal_promoter_snps).group_by(Gene),
            total=q.count()
        ):
            gene.snps_count = count
        session.commit()
        for gene in Gene.query.filter(Gene.snps_count.is_(None)):
            gene.snps_count = 0
        session.commit()
        session.close()

    if PROMOTER_GENE_COUNT_010:
        print('Updating promoter snp count 010')
        q = session.query(Gene, db.func.count('*')).join(SNP, Gene.proximal_promoter_snps).filter(SNP.fdr_class.in_(['0.01', '0.05', '0.1'])).group_by(Gene)
        for gene, count in tqdm(
            session.query(Gene, db.func.count('*')).join(SNP, Gene.proximal_promoter_snps).filter(
                SNP.fdr_class.in_(['0.01', '0.05', '0.1'])).group_by(Gene),
            total=q.count()
        ):
            gene.snps_count010 = count
        session.commit()
        for gene in Gene.query.filter(Gene.snps_count010.is_(None)):
            gene.snps_count010 = 0
        session.commit()
        session.close()

    if SET_NONE_TO_ZERO:
        print('Setting 0 when NULL')
        items_dict = {
            Gene: ['snps_count', 'snps_count010', 'eqtl_snps_count', 'eqtl_snps_count010'],
            # SNP: ['has_clinvar_associations', 'has_phewas_associations', 'has_ebi_associations',
            #       'has_qtl_associations', 'has_grasp_associations', 'has_finemapping_associations',
            #       'has_concordance'],
            Dnase: ['aggregated_snps_count', 'aggregated_snps_count005', 'aggregated_snps_count010'],
            Atac: ['aggregated_snps_count', 'aggregated_snps_count005', 'aggregated_snps_count010'],
            Faire: ['aggregated_snps_count', 'aggregated_snps_count005', 'aggregated_snps_count010']
        }
        for cls, lst in items_dict.items():
            for field in lst:
                print(cls.__name__, field)
                for item in cls.query.filter(getattr(cls, field).is_(None)):
                    setattr(item, field, 0)
                session.commit()

    if CHECK_NONE:
        print('Performing NULL checks')
        items_dict = {
            #Gene: ['snps_count', 'snps_count010', 'eqtl_snps_count', 'eqtl_snps_count010'],
            SNP: ['best_p_value', 'fdr_class', 'best_es', 'es_class', 'context', 'has_clinvar_associations', 'has_phewas_associations',
                  'has_ebi_associations', 'has_qtl_associations', 'has_grasp_associations',
                  'has_finemapping_associations', 'has_concordance'],
            #TranscriptionFactor: ['aggregated_snps_count', 'aggregated_snps_count005', 'aggregated_snps_count010'],
            Atac: ['aggregated_snps_count', 'aggregated_snps_count005', 'aggregated_snps_count010'],
            #TranscriptionFactorSNP: ['best_p_value', 'fdr_class', 'best_es', 'es_class'],
            AtacSNP: ['best_p_value', 'fdr_class', 'best_es', 'es_class'],
            Dnase: ['aggregated_snps_count', 'aggregated_snps_count005', 'aggregated_snps_count010'],
            # TranscriptionFactorSNP: ['best_p_value', 'fdr_class', 'best_es', 'es_class'],
            DnaseSNP: ['best_p_value', 'fdr_class', 'best_es', 'es_class'],
            Faire: ['aggregated_snps_count', 'aggregated_snps_count005', 'aggregated_snps_count010'],
            # TranscriptionFactorSNP: ['best_p_value', 'fdr_class', 'best_es', 'es_class'],
            FaireSNP: ['best_p_value', 'fdr_class', 'best_es', 'es_class'],
        }
        for cls, lst in items_dict.items():
            for field in lst:
                print(cls.__name__, field)
                ctr = cls.query.filter(getattr(cls, field).is_(None)).count()
                if ctr > 0:
                    print('WARN: {} occasions of {}.{} = NULL'.format(ctr, cls, field))

