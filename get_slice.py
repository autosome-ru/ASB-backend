from ASB_app import *
from ASB_app.releases import ReleaseFord
from sqlalchemy.orm import aliased
import os

current_release = ReleaseFord

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

session = current_release.session
db =current_release.db

grasp = aliased(Phenotype, name='grasp')
ebi = aliased(Phenotype, name='ebi')
clinvar = aliased(Phenotype, name='clinvar')
finemapping = aliased(Phenotype, name='finemapping')
qtl = aliased(Phenotype, name='qtl')
phewas = aliased(Phenotype, name='phewas')


def pack(line):
    return '\t'.join(map(str, line)) + '\n'


out_path = os.path.expanduser('~/Desktop/ADASTRA_ford_slice')
if not os.path.isdir(out_path):
    os.mkdir(out_path)

tf_path = os.path.join(out_path, 'TF')
if not os.path.isdir(tf_path):
    os.mkdir(tf_path)

cl_path = os.path.join(out_path, 'CL')
if not os.path.isdir(cl_path):
    os.mkdir(cl_path)

common_header = ['CHROMOSOME', 'POSITION', 'RS_ID', 'REF', 'ALT',
             'PEAK_CALLS', 'MEAN_BAD', 'LOG10_FDR_REF', 'LOG10_FDR_ALT',
             'EFFECT_SIZE_REF', 'EFFECT_SIZE_ALT']
cl_header = common_header + ['AGGREGATED_TFS']
tf_header = common_header + ['MOTIF_LOG_P_REF', 'MOTIF_LOG_P_ALT', 'MOTIF_LOG2_FC', 'MOTIF_POSITION',
             'MOTIF_ORIENTATION', 'MOTIF_CONCORDANCE', 'AGGREGATED_CELL_TYPES']


def process_row(row, param):
    header = {'TF': tf_header, 'CL': cl_header}[param]
    row_dict = dict(zip(header, row))
    row_dict['PEAK_CALLS'] = 'None' if row_dict['PEAK_CALLS'] == 0 else 'Single' if row_dict['PEAK_CALLS'] == 1 else 'Multiple'
    row_dict['RS_ID'] = 'rs' + str(row_dict['RS_ID'])
    if param == 'TF':
        row_dict['MOTIF_ORIENTATION'] = '+' if row_dict['MOTIF_ORIENTATION'] else '-' if row_dict['MOTIF_ORIENTATION'] == 0 else 'None'
        row_dict['MOTIF_CONCORDANCE'] = {
            'Concordant': 'Concordant',
            'Discordant': 'Discordant',
            'Weak Concordant': 'Weak concordant',
            'Weak Discordant': 'Weak discordant',
            'No Hit': 'No hit',
            None: 'None',
        }[row_dict['MOTIF_CONCORDANCE']]
    return tuple(row_dict[name] for name in header)


for tf in TranscriptionFactor.query:
    if tf.aggregated_snps_count == 0:
        continue
    id = tf.tf_id
    q_tf = session.query(
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.chromosome)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.position)),
        db.func.group_concat(db.func.distinct(SNP.rs_id)),
        db.func.group_concat(db.func.distinct(SNP.ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.peak_calls)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.mean_bad)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.log_p_value_ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.log_p_value_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.es_ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.es_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_log_p_ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_log_p_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_log_2_fc)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_position)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_orientation)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.motif_concordance)),
        db.func.group_concat(db.func.distinct(CellLine.name)),
        # db.func.group_concat(db.func.distinct(qtl.phenotype_name)),
        # db.func.group_concat(db.func.distinct(ebi.phenotype_name)),
        # db.func.group_concat(db.func.distinct(phewas.phenotype_name)),
        # db.func.group_concat(db.func.distinct(finemapping.phenotype_name)),
        # db.func.group_concat(db.func.distinct(grasp.phenotype_name)),
        # db.func.group_concat(db.func.distinct(clinvar.phenotype_name)),
    ).filter(
        TranscriptionFactorSNP.tf_id == id
    ).join(
        SNP,
        TranscriptionFactorSNP.snp
    # ).join(
    #     PhenotypeSNPCorrespondence,
    #     (SNP.chromosome == PhenotypeSNPCorrespondence.chromosome) &
    #     (SNP.position == PhenotypeSNPCorrespondence.position) &
    #     (SNP.alt == PhenotypeSNPCorrespondence.alt)
    ).join(
        ExpSNP,
        TranscriptionFactorSNP.exp_snps
    ).join(
        CellLine,
        ExpSNP.cell_line
    # ).join(
    #     qtl,
    #     (PhenotypeSNPCorrespondence.phenotype_id == qtl.phenotype_id) &
    #     (qtl.db_name == 'QTL'),
    #     isouter=True
    # ).join(
    #     ebi,
    #     (PhenotypeSNPCorrespondence.phenotype_id == ebi.phenotype_id) &
    #     (ebi.db_name == 'ebi'),
    #     isouter=True
    # ).join(
    #     phewas,
    #     (PhenotypeSNPCorrespondence.phenotype_id == phewas.phenotype_id) &
    #     (phewas.db_name == 'phewas'),
    #     isouter=True
    # ).join(
    #     finemapping,
    #     (PhenotypeSNPCorrespondence.phenotype_id == finemapping.phenotype_id) &
    #     (finemapping.db_name == 'finemapping'),
    #     isouter=True
    # ).join(
    #     grasp,
    #     (PhenotypeSNPCorrespondence.phenotype_id == grasp.phenotype_id) &
    #     (grasp.db_name == 'grasp'),
    #     isouter=True
    # ).join(
    #     clinvar,
    #     (PhenotypeSNPCorrespondence.phenotype_id == clinvar.phenotype_id) &
    #     (clinvar.db_name == 'clinvar'),
    #     isouter=True
    ).group_by(TranscriptionFactorSNP.tf_snp_id)

    if q_tf.count() == 0:
        continue

    with open(os.path.join(tf_path, tf.name + '.tsv'), 'w') as out:
        out.write(pack(tf_header))
        for tup in q_tf:
            out.write(pack(process_row(tup, 'TF')))

for cl in CellLine.query:
    if cl.aggregated_snps_count == 0:
        continue
    id = cl.cl_id

    q_cl = session.query(
        db.func.group_concat(db.func.distinct(CellLineSNP.chromosome)),
        db.func.group_concat(db.func.distinct(CellLineSNP.position)),
        db.func.group_concat(db.func.distinct(SNP.rs_id)),
        db.func.group_concat(db.func.distinct(SNP.ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.alt)),
        db.func.group_concat(db.func.distinct(CellLineSNP.peak_calls)),
        db.func.group_concat(db.func.distinct(CellLineSNP.mean_bad)),
        db.func.group_concat(db.func.distinct(CellLineSNP.log_p_value_ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.log_p_value_alt)),
        db.func.group_concat(db.func.distinct(CellLineSNP.es_ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.es_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactor.name)),
        # db.func.group_concat(db.func.distinct(qtl.phenotype_name)),
        # db.func.group_concat(db.func.distinct(ebi.phenotype_name)),
        # db.func.group_concat(db.func.distinct(phewas.phenotype_name)),
        # db.func.group_concat(db.func.distinct(finemapping.phenotype_name)),
        # db.func.group_concat(db.func.distinct(grasp.phenotype_name)),
        # db.func.group_concat(db.func.distinct(clinvar.phenotype_name)),
    ).filter(
        CellLineSNP.cl_id == id
    ).join(
        SNP,
        CellLineSNP.snp
    # ).join(
    #     PhenotypeSNPCorrespondence,
    #     (SNP.chromosome == PhenotypeSNPCorrespondence.chromosome) &
    #     (SNP.position == PhenotypeSNPCorrespondence.position) &
    #     (SNP.alt == PhenotypeSNPCorrespondence.alt)
    ).join(
        ExpSNP,
        CellLineSNP.exp_snps
    ).join(
        TranscriptionFactor,
        ExpSNP.transcription_factor
    # ).join(
    #     qtl,
    #     (PhenotypeSNPCorrespondence.phenotype_id == qtl.phenotype_id) &
    #     (qtl.db_name == 'QTL'),
    #     isouter=True
    # ).join(
    #     ebi,
    #     (PhenotypeSNPCorrespondence.phenotype_id == ebi.phenotype_id) &
    #     (ebi.db_name == 'ebi'),
    #     isouter=True
    # ).join(
    #     phewas,
    #     (PhenotypeSNPCorrespondence.phenotype_id == phewas.phenotype_id) &
    #     (phewas.db_name == 'phewas'),
    #     isouter=True
    # ).join(
    #     finemapping,
    #     (PhenotypeSNPCorrespondence.phenotype_id == finemapping.phenotype_id) &
    #     (finemapping.db_name == 'finemapping'),
    #     isouter=True
    # ).join(
    #     grasp,
    #     (PhenotypeSNPCorrespondence.phenotype_id == grasp.phenotype_id) &
    #     (grasp.db_name == 'grasp'),
    #     isouter=True
    # ).join(
    #     clinvar,
    #     (PhenotypeSNPCorrespondence.phenotype_id == clinvar.phenotype_id) &
    #     (clinvar.db_name == 'clinvar'),
    #     isouter=True
    ).group_by(CellLineSNP.cl_snp_id)

    if q_cl.count() == 0:
        continue

    if cl.name == 'SLK (Clear cell renal cell carcinoma. Derived from metastatic site: Skin.Clear cell renal cell carcinoma (NCIt: C4033) Derived from metastatic site: Skin)':
        clname = 'SLK (Clear cell renal cell carcinoma. Derived from metastatic site: Skin)'
    else:
        clname = cl.name

    with open(os.path.join(cl_path, clname.replace('/', '-') + '.tsv'), 'w') as out:
        out.write(pack(cl_header))
        for tup in q_cl:
            out.write(pack(process_row(tup, 'CL')))
