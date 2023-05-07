import json

from ASB_app import releases
from ASB_app.constants import fdr_classes, es_classes, fdr_choices, ananastra_stats_file, chromosomes
from ASB_app.models import CandidateSNP, CandidateRS, CandidateCLRS, CandidateTFRS
import numpy as np

current_release = releases.current_release

default_fdr_tr = 0.05


# CellLine, CellLineSNP, \
# SNP, ExpSNP, Phenotype, PhenotypeSNPCorrespondence, Gene, Experiment = \
#     current_release.CellLine, current_release.CellLineSNP, \
#     current_release.SNP, current_release.ExpSNP, current_release.Phenotype, current_release.PhenotypeSNPCorrespondence, current_release.Gene, current_release.Experiment


def get_fdr_filters(alternative, fdr, snp_class):
    if alternative == 'less':
        return [snp_class.fdr_class.in_(get_corresponding_fdr_classes(fdr, low=True))]
    elif alternative == 'greater':
        return [snp_class.fdr_class.in_(get_corresponding_fdr_classes(fdr, low=False))]
    elif alternative == 'all':
        return []
    else:
        raise ValueError


def get_expected_asbs(fdr_class, level='TF', ag_id=None, mode='all'):
    return
    if level == 'ALL':
        SNPClass = SNP
    elif level == 'CL':
        SNPClass = CellLineSNP
    else:
        raise ValueError

    filters = []
    if ag_id is not None:
        if ag_id in chromosomes:
            filters += [SNPClass.chromosome == ag_id]
        else:
            filters += [getattr(SNPClass, {'TF': 'tf_id', 'CL': 'cl_id'}[level]) == ag_id]


    filters += get_fdr_filters('greater', fdr_class, SNPClass)
    q = SNPClass.query.filter(*filters)

    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_fdr_class(logp):
    for fdr in fdr_classes:
        if logp >= -np.log10(float(fdr)):
            return fdr


def get_es_class(effect_size):
    if effect_size is not None:
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


def get_stats_dict(fdrs, level='ALL'):
    return
    if level == 'ALL':
        ag_ids = [None]
    elif level == 'CHR':
        ag_ids = list(chromosomes)
    elif level == 'CL':
        ag_ids = [x.cl_id for x in CellLine.query]
    else:
        raise ValueError
    ret_val = {}
    print(level)
    for ag_id in ag_ids:
        print(ag_id)
        stats_dict = {}
        for fdr_class in fdrs:
            stats_dict[fdr_class] = {}
            print('Collecting statistics for: {}'.format(fdr_class))
            if level in ('ALL', 'CHR', 'TF'):
                expected_tf_asbs_list = get_expected_asbs(fdr_class, ag_id=ag_id, level='TF', mode='all')
                stats_dict[fdr_class]['expected_tf_asbs'] = len(expected_tf_asbs_list)
                stats_dict[fdr_class]['expected_tf_asbs_rs'] = len(set(x.snp.rs_id for x in expected_tf_asbs_list))
            if level in ('ALL', 'CHR', 'CL'):
                expected_cl_asbs_list = get_expected_asbs(fdr_class, ag_id=ag_id, level='CL', mode='all')
                stats_dict[fdr_class]['expected_cl_asbs'] = len(expected_cl_asbs_list)
                stats_dict[fdr_class]['expected_cl_asbs_rs'] = len(set(x.snp.rs_id for x in expected_cl_asbs_list))
            if level in ('ALL', 'CHR'):
                expected_all_asbs_list = get_expected_asbs(fdr_class, ag_id=ag_id, level='ALL', mode='all')
                stats_dict[fdr_class]['expected_all_asbs'] = len(expected_tf_asbs_list) + len(expected_cl_asbs_list)
                stats_dict[fdr_class]['expected_all_asbs_rs'] = len(set(x.rs_id for x in expected_all_asbs_list))
        print('Collecting statistics for candidate data')
        if level in ('ALL', 'CHR'):
            if level == 'CHR':
                add_filters = [CandidateSNP.chromosome == ag_id]
            else:
                add_filters = []
            tf_candidates = CandidateSNP.query.filter(
                CandidateSNP.ag_level == 'TF',
                *add_filters
            ).count()
            cl_candidates = CandidateSNP.query.filter(
                CandidateSNP.ag_level == 'CL',
                *add_filters
            ).count()
            if level == 'CHR':
                candidates_rs = current_release.session.query(CandidateSNP.rs_id)\
                    .filter(*add_filters)\
                    .group_by(CandidateSNP.rs_id).count()
            else:
                candidates_rs = CandidateRS.query.count()
                try:
                    assert CandidateCLRS.query.count() == candidates_rs
                    assert CandidateTFRS.query.count() == candidates_rs
                except AssertionError:
                    print(CandidateCLRS.query.filter(*add_filters).count(),
                          CandidateTFRS.query.filter(*add_filters).count(),
                          candidates_rs)
            stats_dict['1'] = {
                'total_tf_candidates': tf_candidates,
                'total_cl_candidates': cl_candidates,
                'total_all_candidates': tf_candidates + cl_candidates,
                'total_tf_candidates_rs': candidates_rs,
                'total_cl_candidates_rs': candidates_rs,
                'total_all_candidates_rs': candidates_rs,
            }
        else:
            stats_dict['1'] = {
                'total_{}_candidates'.format(level.lower()): CandidateSNP.query.filter(
                    CandidateSNP.ag_level == level,
                    CandidateSNP.ag_id == ag_id,
                ).count(),
                'total_{}_candidates_rs'.format(level.lower()): current_release.session.query(
                    CandidateSNP.rs_id.distinct()).filter(
                    CandidateSNP.ag_level == level,
                    CandidateSNP.ag_id == ag_id,
                ).count()
            }
        ret_val[ag_id] = stats_dict
    if level == 'ALL':
        ret_val = ret_val[None]
    return ret_val


def collect_ananastra_stats(fdrs=fdr_choices, level=None):
    if level is None:
        with open(ananastra_stats_file, 'w') as f:
            json.dump({
                'chr_stats_dict': get_stats_dict(fdrs, level='CHR'),
                'cl_stats_dict': get_stats_dict(fdrs, level='CL'),
                'stats_dict': get_stats_dict(fdrs),
            }, f, indent=2)
    else:
        with open(ananastra_stats_file.replace('ananastra_constants', 'ananastra_constants_{}'.format(level)), 'w') as f:
            json.dump(get_stats_dict(fdrs, level=level), f, indent=2)
