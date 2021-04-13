from ASB_app import releases
from ASB_app.constants import fdr_classes, es_classes
from ASB_app.models import CandidateSNP
import numpy as np

current_release = releases.current_release

default_fdr_tr = 0.05

if current_release.name != 'dnase':
    TranscriptionFactor, TranscriptionFactorSNP, CellLine, CellLineSNP, \
    SNP, ExpSNP, Phenotype, PhenotypeSNPCorrespondence, Gene, Experiment = \
        current_release.TranscriptionFactor, current_release.TranscriptionFactorSNP, current_release.CellLine, current_release.CellLineSNP, \
        current_release.SNP, current_release.ExpSNP, current_release.Phenotype, current_release.PhenotypeSNPCorrespondence, current_release.Gene, current_release.Experiment
else:
    CellLine, CellLineSNP, \
    SNP, ExpSNP, Phenotype, PhenotypeSNPCorrespondence, Gene, Experiment = \
        current_release.CellLine, current_release.CellLineSNP, \
        current_release.SNP, current_release.ExpSNP, current_release.Phenotype, current_release.PhenotypeSNPCorrespondence, current_release.Gene, current_release.Experiment


def get_possible_tf_asbs(fdr, mode='all'):
    q = TranscriptionFactorSNP.query.filter(
        TranscriptionFactorSNP.best_p_value >= fdr,
    )
    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_possible_cl_asbs(fdr, mode='all'):
    q = CellLineSNP.query.filter(
        CellLineSNP.best_p_value >= fdr,
    )
    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_possible_all_asbs(fdr, mode='all'):
    q = SNP.query.filter(
        SNP.best_p_value >= fdr
    )
    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_possible_tf_candidates(fdr_low, mode='all'):
    q = CandidateSNP.query.filter(
        CandidateSNP.ag_level == 'TF',
        CandidateSNP.best_p_value < fdr_low,
    )
    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_possible_cl_candidates(fdr_low, mode='all'):
    q = CandidateSNP.query.filter(
        CandidateSNP.ag_level == 'CL',
        CandidateSNP.best_p_value < fdr_low,
    )
    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_possible_all_candidates(fdr_low, mode='all'):
    q = CandidateSNP.query.filter(
        CandidateSNP.best_p_value < fdr_low
    )
    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_fdr_class(logp):
    for fdr in fdr_classes:
        if logp >= -np.log10(float(fdr)):
            return fdr


def get_es_class(effect_size):
    for es in es_classes[:-1]:
        if effect_size >= float(es):
            return es
    return es_classes[-1]


def get_corresponding_fdr_classes(fdr_class, low=False):
    index = fdr_classes.index(fdr_class)
    if low:
        return fdr_classes[index + 1:]
    else:
        return fdr_classes[:index + 1]


def get_corresponding_es_classes(es_class, low=False):
    index = es_classes.index(es_class)
    if low:
        return es_classes[index + 1:]
    else:
        return es_classes[:index + 1]


def get_stats_dict(fdrs):
    stats_dict = {}
    for fdr_class in fdrs:
        print(fdr_class)
        fdr = -np.log10(float(fdr_class))
        possible_tf_asbs_lsit = get_possible_tf_asbs(fdr)
        print('tf_asb')
        possible_cl_asbs_list = get_possible_cl_asbs(fdr)
        print('cl_asb')
        possible_all_asbs_list = get_possible_all_asbs(fdr)
        print('all_asb')
        stats_dict[fdr_class] = {
            'possible_tf_asbs': len(possible_tf_asbs_lsit),
            'possible_cl_asbs': len(possible_cl_asbs_list),
            'possible_all_asbs': len(possible_tf_asbs_lsit) + len(possible_cl_asbs_list),
            'possible_tf_asbs_rs': len(set(x.snp.rs_id for x in possible_tf_asbs_lsit)),
            'possible_cl_asbs_rs': len(set(x.snp.rs_id for x in possible_cl_asbs_list)),
            'possible_all_asbs_rs': len(set(x.rs_id for x in possible_all_asbs_list)),
        }
    return stats_dict
