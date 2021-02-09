from ASB_app import releases
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


def get_stats_dict(fdrs):
    stats_dict = {}
    for fdr_raw in fdrs:
        print(fdr_raw)
        fdr = -np.log10(fdr_raw)
        possible_tf_asbs_lsit = get_possible_tf_asbs(fdr)
        print('tf_asb')
        possible_cl_asbs_list = get_possible_cl_asbs(fdr)
        print('cl_asb')
        possible_all_asbs_list = get_possible_all_asbs(fdr)
        print('all_asb')
        stats_dict[fdr_raw] = {
            'possible_tf_asbs': len(possible_tf_asbs_lsit),
            'possible_cl_asbs': len(possible_cl_asbs_list),
            'possible_all_asbs': len(possible_all_asbs_list),
            'possible_tf_asbs_rs': len(set(x.snp.rs_id for x in possible_tf_asbs_lsit)),
            'possible_cl_asbs_rs': len(set(x.snp.rs_id for x in possible_cl_asbs_list)),
            'possible_all_asbs_rs': len(set(x.rs_id for x in possible_all_asbs_list)),
        }
    return stats_dict
