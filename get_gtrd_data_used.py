from ASB_app import *
from ASB_app.releases import current_release

T = current_release.TranscriptionFactor
E = current_release.Experiment
C = current_release.CellLine
B = current_release.BADGroup
session = current_release.session


with open('D:\\Sashok\\exps_{}.tsv'.format(current_release.name), 'w', encoding='utf-8') as f:
    f.write('EXP	ALIGN	TF_UNIPROT_AC	TF_UNIPROT_ID	GTRD_CELL_TYPE_ID	GTRD_CELL_TYPE_NAME	BAD_GROUP_ID	BAD_GROUP_NAME	ENCODE	GEO_GSE	IS_CONTROL\n')
    for e, t, c, b in session.query(E, T, C, B).join(T, E.transcription_factor).join(C, E.cell_line).join(B, E.bad_group):
        f.write('\t'.join(map(str, [e.exp_id, e.align, t.uniprot_ac, t.name, c.cl_id, c.name, b.bad_group_id, b.bad_group_name, e.encode, e.geo_gse, e.is_control])) + '\n')