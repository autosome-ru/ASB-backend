import gzip
import os
from sqlalchemy import tuple_, or_
import re
from datetime import datetime
from ASB_app import logger, executor
from ASB_app.constants import stats_dict, tf_stats_dict, cl_stats_dict, chr_stats_dict, chromosomes, max_nrows, \
    max_comments
from ASB_app.service import ananastra_service
from ASB_app.utils import pack, process_row, group_concat_distinct_sep
from ASB_app.utils.statistics import get_stats_dict, get_corresponding_fdr_classes
from sqlalchemy.orm import aliased
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, distributions
from statsmodels.stats.multitest import multipletests

from ASB_app.models import CandidateSNP, CandidateRS, CandidateTFRS, CandidateCLRS, PositionHash, LDIslandsInfo

from ASB_app.releases import current_release

session = current_release.session
db = current_release.db

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


class ConvError(ValueError):
    pass


class TooBigError(ValueError):
    pass


def logit_combine_p_values(pvalues):
    pvalues = np.array([pvalue for pvalue in pvalues if 1 > pvalue > 0])
    if len(pvalues) == 0:
        return 1
    elif len(pvalues) == 1:
        return pvalues[0]

    statistic = -np.sum(np.log(pvalues)) + np.sum(np.log1p(-pvalues))
    k = len(pvalues)
    nu = np.int_(5 * k + 4)
    approx_factor = np.sqrt(np.int_(3) * nu / (np.int_(k) * np.square(np.pi) * (nu - np.int_(2))))
    pval = distributions.t.sf(statistic * approx_factor, nu)
    return pval


def convert_rs_to_int(rs_str):
    if pd.isna(rs_str):
        raise ConvError(rs_str)
    rs_str = rs_str.strip()
    if not re.match(r'^rs\d+$', rs_str):
        raise ConvError(rs_str)
    return int(rs_str[2:])


def get_tf_query(rs_ids, fdr):
    grasp = aliased(Phenotype, name='grasp')
    ebi = aliased(Phenotype, name='ebi')
    clinvar = aliased(Phenotype, name='clinvar')
    finemapping = aliased(Phenotype, name='finemapping')
    qtl = aliased(Phenotype, name='qtl')
    phewas = aliased(Phenotype, name='phewas')

    return session.query(
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.chromosome)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.position)),
        db.func.group_concat(db.func.distinct(SNP.rs_id)),
        db.func.group_concat(db.func.distinct(SNP.ref)),
        db.func.group_concat(db.func.distinct(TranscriptionFactorSNP.alt)),
        db.func.group_concat(db.func.distinct(SNP.context)),
        db.func.group_concat(db.func.distinct(TranscriptionFactor.name)),
        db.func.group_concat(db.func.distinct(TranscriptionFactor.uniprot_ac)),
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
        db.func.group_concat(db.func.distinct(CellLine.name), separator=', '),
        db.func.group_concat(db.func.distinct(qtl.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(ebi.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(phewas.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(finemapping.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(grasp.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(clinvar.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(Gene.gene_name), separator=', '),
    ).join(
        SNP,
        TranscriptionFactorSNP.snp
    ).filter(
        SNP.rs_id.in_(rs_ids),
        TranscriptionFactorSNP.fdr_class.in_(get_corresponding_fdr_classes(fdr))
    ).join(
        TranscriptionFactor,
        TranscriptionFactorSNP.transcription_factor
    ).join(
        PhenotypeSNPCorrespondence,
        (SNP.chromosome == PhenotypeSNPCorrespondence.chromosome) &
        (SNP.position == PhenotypeSNPCorrespondence.position) &
        (SNP.alt == PhenotypeSNPCorrespondence.alt),
        isouter=True
    ).join(
        ExpSNP,
        TranscriptionFactorSNP.exp_snps,
        isouter=True
    ).filter(
        ((ExpSNP.p_value_ref - ExpSNP.p_value_alt) * (
                TranscriptionFactorSNP.log_p_value_alt - TranscriptionFactorSNP.log_p_value_ref) > 0) |
        (ExpSNP.p_value_ref == None)
    ).join(
        Experiment,
        ExpSNP.experiment,
        isouter=True
    ).join(
        CellLine,
        Experiment.cell_line,
        isouter=True
    ).join(
        Gene,
        SNP.target_genes,
        isouter=True
    ).join(
        qtl,
        (PhenotypeSNPCorrespondence.phenotype_id == qtl.phenotype_id) &
        (qtl.db_name == 'QTL'),
        isouter=True
    ).join(
        ebi,
        (PhenotypeSNPCorrespondence.phenotype_id == ebi.phenotype_id) &
        (ebi.db_name == 'ebi'),
        isouter=True
    ).join(
        phewas,
        (PhenotypeSNPCorrespondence.phenotype_id == phewas.phenotype_id) &
        (phewas.db_name == 'phewas'),
        isouter=True
    ).join(
        finemapping,
        (PhenotypeSNPCorrespondence.phenotype_id == finemapping.phenotype_id) &
        (finemapping.db_name == 'finemapping'),
        isouter=True
    ).join(
        grasp,
        (PhenotypeSNPCorrespondence.phenotype_id == grasp.phenotype_id) &
        (grasp.db_name == 'grasp'),
        isouter=True
    ).join(
        clinvar,
        (PhenotypeSNPCorrespondence.phenotype_id == clinvar.phenotype_id) &
        (clinvar.db_name == 'clinvar'),
        isouter=True
    ).group_by(TranscriptionFactorSNP.tf_snp_id)


def get_cl_query(rs_ids, fdr):
    grasp = aliased(Phenotype, name='grasp')
    ebi = aliased(Phenotype, name='ebi')
    clinvar = aliased(Phenotype, name='clinvar')
    finemapping = aliased(Phenotype, name='finemapping')
    qtl = aliased(Phenotype, name='qtl')
    phewas = aliased(Phenotype, name='phewas')

    return session.query(
        db.func.group_concat(db.func.distinct(CellLineSNP.chromosome)),
        db.func.group_concat(db.func.distinct(CellLineSNP.position)),
        db.func.group_concat(db.func.distinct(SNP.rs_id)),
        db.func.group_concat(db.func.distinct(SNP.ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.alt)),
        db.func.group_concat(db.func.distinct(SNP.context)),
        db.func.group_concat(db.func.distinct(CellLine.name)),
        db.func.group_concat(db.func.distinct(CellLine.cl_id)),
        db.func.group_concat(db.func.distinct(CellLineSNP.peak_calls)),
        db.func.group_concat(db.func.distinct(CellLineSNP.mean_bad)),
        db.func.group_concat(db.func.distinct(CellLineSNP.log_p_value_ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.log_p_value_alt)),
        db.func.group_concat(db.func.distinct(CellLineSNP.es_ref)),
        db.func.group_concat(db.func.distinct(CellLineSNP.es_alt)),
        db.func.group_concat(db.func.distinct(TranscriptionFactor.name), separator=', '),
        db.func.group_concat(db.func.distinct(qtl.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(ebi.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(phewas.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(finemapping.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(grasp.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(clinvar.phenotype_name), separator=', '),
        db.func.group_concat(db.func.distinct(Gene.gene_name), separator=', '),
    ).join(
        SNP,
        CellLineSNP.snp
    ).filter(
        SNP.rs_id.in_(rs_ids),
        CellLineSNP.fdr_class.in_(get_corresponding_fdr_classes(fdr))
    ).join(
        CellLine,
        CellLineSNP.cell_line
    ).join(
        PhenotypeSNPCorrespondence,
        (SNP.chromosome == PhenotypeSNPCorrespondence.chromosome) &
        (SNP.position == PhenotypeSNPCorrespondence.position) &
        (SNP.alt == PhenotypeSNPCorrespondence.alt),
        isouter=True
    ).join(
        ExpSNP,
        CellLineSNP.exp_snps,
        isouter=True
    ).filter(
        ((ExpSNP.p_value_ref - ExpSNP.p_value_alt) * (CellLineSNP.log_p_value_alt - CellLineSNP.log_p_value_ref) > 0) |
        (ExpSNP.p_value_ref == None)
    ).join(
        Experiment,
        ExpSNP.experiment,
        isouter=True
    ).join(
        TranscriptionFactor,
        Experiment.transcription_factor,
        isouter=True
    ).join(
        Gene,
        SNP.target_genes,
        isouter=True
    ).join(
        qtl,
        (PhenotypeSNPCorrespondence.phenotype_id == qtl.phenotype_id) &
        (qtl.db_name == 'QTL'),
        isouter=True
    ).join(
        ebi,
        (PhenotypeSNPCorrespondence.phenotype_id == ebi.phenotype_id) &
        (ebi.db_name == 'ebi'),
        isouter=True
    ).join(
        phewas,
        (PhenotypeSNPCorrespondence.phenotype_id == phewas.phenotype_id) &
        (phewas.db_name == 'phewas'),
        isouter=True
    ).join(
        finemapping,
        (PhenotypeSNPCorrespondence.phenotype_id == finemapping.phenotype_id) &
        (finemapping.db_name == 'finemapping'),
        isouter=True
    ).join(
        grasp,
        (PhenotypeSNPCorrespondence.phenotype_id == grasp.phenotype_id) &
        (grasp.db_name == 'grasp'),
        isouter=True
    ).join(
        clinvar,
        (PhenotypeSNPCorrespondence.phenotype_id == clinvar.phenotype_id) &
        (clinvar.db_name == 'clinvar'),
        isouter=True
    ).group_by(CellLineSNP.cl_snp_id)


def get_asbs_by_level(rs_ids=None, fdr='0.05', level='TF', mode='all'):
    if level == 'ALL':
        SNPClass = SNP
    elif level == 'TF':
        SNPClass = TranscriptionFactorSNP
    elif level == 'CL':
        SNPClass = CellLineSNP
    else:
        raise ValueError

    fdr_filters = get_fdr_filters('greater', fdr, SNPClass)
    if level == 'ALL':
        rs_filters = get_rs_filters(rs_ids, SNPClass)
        q = SNPClass.query.filter(*(fdr_filters + rs_filters))
    else:
        rs_filters = get_rs_filters(rs_ids, SNP)
        q = SNPClass.query.filter(*fdr_filters).join(
            SNP, SNPClass.snp
        ).filter(*rs_filters)

    if mode == 'count':
        return q.count()
    elif mode == 'all':
        return q.all()


def get_candidates_by_level(rs_ids=None, fdr='0.05', intervals=None, lds=None, ld_type='ld_eur', level='TF',
                            alternative='less', mode='all', rs=False):
    if rs:
        if level == 'TF':
            SNPClass = CandidateTFRS
        elif level == 'CL':
            SNPClass = CandidateCLRS
        elif level == 'ALL':
            SNPClass = CandidateRS
        else:
            raise ValueError
    else:
        SNPClass = CandidateSNP

    filters = get_rs_filters(rs_ids, SNPClass)
    filters += get_fdr_filters(alternative, fdr, SNPClass)

    if not rs and level != 'ALL':
        filters += [CandidateSNP.ag_level == level]

    q = SNPClass.query.filter(*filters)

    return_empty = False

    if lds is not None:
        assert ld_type in ('ld_eur', 'ld_asn', 'ld_afr')
        assert intervals is None
        if len(lds) == 0:
            return_empty = True
        else:
            q = q.join(LDIslandsInfo, SNPClass.rs_id == LDIslandsInfo.rs_id).filter(
                *get_ld_filters(lds, ld_type)
            )
    elif intervals is not None:
        if len(intervals) == 0:
            return_empty = True
        else:
            q = q.join(PositionHash, SNPClass.rs_id == PositionHash.rs_id).filter(
                *get_intervals_filters(intervals)
            )

    if return_empty:
        if mode == 'count':
            return 0
        elif mode == 'all':
            return []
    else:
        if mode == 'count':
            return q.count()
        elif mode == 'all':
            return q.all()


def get_rs_filters(rs_ids, snp_class):
    if rs_ids is not None:
        return [snp_class.rs_id.in_(rs_ids)]
    else:
        return []


def get_intervals_filters(intervals):
    # print(*['OR position_hash.position_hash BETWEEN {} AND {}\n'.format(st, ed) for st, ed in intervals])
    return [or_(*[PositionHash.position_hash.between(st, ed) for st, ed in intervals])]


def get_ld_filters(lds, ld_type):
    return [getattr(LDIslandsInfo, ld_type).in_(lds)]


def get_fdr_filters(alternative, fdr, snp_class):
    if alternative == 'less':
        return [snp_class.fdr_class.in_(get_corresponding_fdr_classes(fdr, low=True))]
    elif alternative == 'greater':
        return [snp_class.fdr_class.in_(get_corresponding_fdr_classes(fdr, low=False))]
    elif alternative == 'all':
        return []
    else:
        raise ValueError


def divide_chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def divide_query(get_query, values, chunk_size=100000):
    for chunk in divide_chunks(values, chunk_size):
        yield get_query(chunk)


def get_alleles(df):
    assert len(set(df['REF'].tolist())) == 1
    ref = df['REF'].tolist()[0]
    alts = list(set(df['ALT'].tolist()))
    return '/'.join([ref] + alts)


def get_preferences(df):
    ref = len(df[df['LOG10_FDR_REF'] > df['LOG10_FDR_ALT']].index) > 0
    alt = len(df[df['LOG10_FDR_ALT'] > df['LOG10_FDR_REF']].index) > 0
    assert ref or alt
    if ref and alt:
        return 'Both'
    elif ref:
        return 'Ref'
    else:
        return 'Alt'


def modify_counts(asb_data, counts=None, top=False):
    if top:
        if not counts:
            return []
        tuples = []
        data_by_name = {data['name']: data for data in asb_data}
        for count in counts:
            data = data_by_name[count['name']]
            tuples.append((count, data))
        counts_list, _ = (list(a) for a in
                          zip(*sorted(tuples, key=lambda x: (x[0]['count'], x[1]['odds']), reverse=True)))
        for item in counts_list:
            # In background it'd be more precise to group expected asbs by rs+alt, but expected_asbs_rs is a good approx.
            item['background_count'] = data_by_name[item['name']]['expected_asbs_rs'] + item['count']
        if len(counts_list) > 7:
            counts_list = counts_list[:6] + [{'name': 'Other', 'count': sum(x['count'] for x in counts_list[6:]),
                                              'background_count': sum(x['background_count'] for x in counts_list[6:])}]
    else:
        counts_list = sorted(asb_data, key=lambda x: (x['asbs'], x['odds']), reverse=True)
        if len(counts_list) > 7:
            counts_list = [{'name': x['name'], 'count': x['asbs'],
                            'background_count': x['expected_asbs'] + x['asbs']} for x in counts_list[:6]] + [
                              {'name': 'Other', 'count': sum(x['asbs'] for x in counts_list[6:]),
                               'background_count': sum(x['asbs'] + x['expected_asbs'] for x in counts_list[6:])}]
        else:
            counts_list = [{'name': x['name'], 'count': x['asbs'],
                            'background_count': x['expected_asbs'] + x['asbs']} for x in counts_list]
    return counts_list


def update_ticket_status(ticket, status):
    session.close()
    session.add(ticket)
    meta_info = dict(ticket.meta_info)
    meta_info.update({
        'status_details': status,
        'last_status_update_at': str(datetime.now()),
    })
    ticket.meta_info = meta_info
    session.commit()


def get_rs_ids_by_chr_pos_query(chromosome, tuples, candidates=False):
    if candidates:
        return CandidateSNP.query.filter(CandidateSNP.chromosome == chromosome,
                                         tuple_(CandidateSNP.position, CandidateSNP.ref, CandidateSNP.alt).in_(
                                             tuples)).all()
    else:
        return SNP.query.filter(SNP.chromosome == chromosome, tuple_(SNP.position, SNP.ref, SNP.alt).in_(tuples)).all()


def get_rs_ids_from_vcf(data):
    if len(data.columns) < 5:
        raise ConvError('number of columns in a VCF file.')
    snps = []
    filtered_data = data[[0, 1, 3, 4]].drop_duplicates()
    unique_submitted_snps_count = len(filtered_data.index)
    all_snps = filtered_data.agg(lambda x: '_'.join(map(str, x)), axis=1)
    for chromosome in data[0].unique():
        if chromosome not in chromosomes:
            fixed_chromosome = f'chr{chromosome}'
            if fixed_chromosome not in chromosomes:
                continue
                # raise ConvError('chromosome: {}'.format(chr))
        try:
            tuples = [(int(position), ref.upper(), alt.upper())
                      for index, (position, ref, alt)
                      in data.loc[data[0] == chromosome, [1, 3, 4]].iterrows()]
        except ValueError as e:
            raise ConvError('position: {}'.format(e.args[0]))
        for snps_chunk in divide_query(lambda poss: get_rs_ids_by_chr_pos_query(fixed_chromosome, poss), tuples,
                                       chunk_size=3000):
            snps += snps_chunk
        for snps_chunk in divide_query(
                lambda poss: get_rs_ids_by_chr_pos_query(fixed_chromosome, poss, candidates=True),
                tuples, chunk_size=3000):
            snps += snps_chunk
    found_snps = set((x.chromosome, x.position, x.ref, x.alt) for x in snps)
    return list(set(x.rs_id for x in snps)), \
           all_snps[~all_snps.isin({'_'.join(map(str, x)) for x in found_snps})].tolist(), \
           unique_submitted_snps_count


def get_snps_from_interval(interval_str):
    interval_str = interval_str.strip()
    match = re.match(r'^(chr)?(\d|1\d|2[0-2]|X|Y):([1-9]\d*)-([1-9]\d*)$', interval_str)
    if match:
        chr, start, end = match.groups()[1:]
        chr = 'chr' + chr
        start = int(start)
        end = int(end)
        if end - start > 10000000:
            raise TooBigError(interval_str)
        return list(set(x for (x,) in session.query(SNP.rs_id).filter(
            SNP.chromosome == chr,
            SNP.position.between(start, end)
        )) | set(x.rs_id for x in CandidateSNP.query.filter(
            CandidateSNP.chromosome == chr,
            CandidateSNP.position.between(start, end)
        )))
    else:
        raise ConvError(interval_str)


def transform_list(obj):
    if len(obj) > 1:
        new_list = []
        merged_elem = None
        for k in range(1, len(obj)):
            if merged_elem is None:
                merged_elem = obj[k - 1]
            if merged_elem[1] >= obj[k][0]:
                merged_elem = (merged_elem[0], obj[k][1])
            else:
                new_list.append(merged_elem)
                merged_elem = None
        if merged_elem is not None:
            new_list.append(merged_elem)
    else:
        return obj
    return new_list


def decode_snv_from_hash(hash):
    chr_index, position = divmod(hash, 10 ** 9)
    return chromosomes[chr_index], position


def get_local_snv_hash(hash, w_size):
    chromosome, position = decode_snv_from_hash(hash)
    start_hash = max(1, position - w_size)
    end_hash = min(position + w_size, 10 ** 9)
    chr_hash = chromosomes.index(chromosome) * 10 ** 9
    return chr_hash + start_hash, chr_hash + end_hash


def snvs_to_hash(hash_list=(), w_size=10 ** 6):
    snvs_list = [get_local_snv_hash(hash, w_size=w_size) for hash in hash_list]
    return transform_list(sorted(snvs_list, key=lambda x: x[0]))


def marshal_inf(odds):
    if np.isinf(odds):
        return 'infinity'
    elif np.isnan(odds):
        return None
    else:
        return str(odds)


def marshal_logp(p):
    if p == 0:
        return 'infinity'
    if np.isnan(p):
        return None
    if p is None:
        return p
    return str(-np.log10(p))


def marshal_data(asb_data):
    return [{k: marshal_inf(v) if k in ('log10_p_value', 'log10_fdr', 'odds') else v for (k, v) in elem.items()} for
            elem in asb_data]


def marshal_chr_data(asb_data):
    return [{k: marshal_inf(v) if k in ('log10_p_value', 'log10_fdr', 'odds',
                                        'tf_log10_p_value', 'tf_log10_fdr', 'tf_odds',
                                        'cl_log10_p_value', 'cl_log10_fdr', 'cl_odds') else v for (k, v) in
             elem.items()} for
            elem in asb_data]


def get_asb_counts_by_chromosome(asbs_list, negatives_list, chromosome, ag_snp=False):
    asbs = len([x for x in asbs_list if x.chromosome == chromosome])
    if ag_snp:
        asbs_rs = len(set(x.snp.rs_id for x in asbs_list if x.chromosome == chromosome))
    else:
        asbs_rs = len(set(x.rs_id for x in asbs_list if x.chromosome == chromosome))
    negatives = len([cand for cand in negatives_list if cand.chromosome == chromosome])
    negatives_rs = len(set(cand.rs_id for cand in negatives_list if cand.chromosome == chromosome))
    return {
        'asbs': asbs,
        'asbs_rs': asbs_rs,
        'negatives': negatives,
        'negatives_rs': negatives_rs,
    }


def get_p_and_odds_by_chromosome(counts_data, expected_asbs_rs, expected_negatives_rs):
    return fisher_exact(
        ((counts_data['asbs_rs'], counts_data['negatives_rs']), (expected_asbs_rs, expected_negatives_rs)),
        alternative='greater')


def read_table_as_df(f):
    not_null_rows = 0
    for index, line in enumerate(f, 1):
        if not line.startswith('#'):
            not_null_rows += 1
        if index - not_null_rows > max_comments or not_null_rows > max_nrows:
            return False, None
    f.seek(0)
    return True, pd.read_table(f,
                               sep='\t',
                               header=None,
                               encoding='utf-8',
                               dtype=str,
                               comment='#',
                               keep_default_na=False)


def read_file_as_df(input_file_name, ticket):
    try:
        with gzip.open(input_file_name, 'rt') as f:
            status, data = read_table_as_df(f)

    except OSError:
        try:
            with open(input_file_name, 'r') as f:
                status, data = read_table_as_df(f)

        except:
            print(input_file_name)
            update_ticket_status(ticket,
                                 'Processing failed: the file must be a valid utf-8 text file with a single SNP rs-ID '
                                 'on each line or a valid .vcf(.gz) file')
            raise ConvError
    return status, data


def get_rs_ids_from_list(rs_list):
    not_found = []
    formatted_rs = []
    for rs_id in rs_list:
        rs_id = rs_id.strip()
        if re.match(r'^rs\d+$', rs_id):
            formatted_rs.append(int(rs_id[2:]))
        elif re.match(r'^\d+$', rs_id):
            formatted_rs.append(int(rs_id))
        else:
            not_found.append(rs_id)
    return formatted_rs, not_found


def check_in_class(rs_id, rs_set_dict, what_for='asb', ag_class='tf'):
    assert ag_class in ('tf', 'cl', 'all')
    assert what_for in ('asb', 'neg', 'non-neg')

    return rs_id in rs_set_dict[what_for][ag_class]


def get_snp_class_label(rs_id, rs_set_dict, not_found=None, ag_class='tf'):
    if not_found is None:
        not_found = set()
    if rs_id in not_found:
        return 'NOT_IN_ADASTRA'
    elif check_in_class(rs_id, rs_set_dict, what_for='asb', ag_class=ag_class):
        return 'ASB'
    elif check_in_class(rs_id, rs_set_dict, what_for='neg', ag_class=ag_class):
        return 'NON-ASB'
    elif check_in_class(rs_id, rs_set_dict, what_for='non-neg', ag_class='all'):
        return 'UNDEFINED'
    else:
        return 'NOT_IN_ADASTRA'


@executor.job
def process_snp_file(ticket_id, fdr_class='0.05', background='WG'):
    processing_start_time = datetime.now()
    input_file_name = ananastra_service.get_path_by_ticket_id(ticket_id)
    ticket = ananastra_service.get_ticket(ticket_id)
    ticket.fdr = fdr_class
    ticket.background = background
    logger.info('start')
    session.commit()
    change_status_on_fail = False

    try:
        logger.info('start parsing')
        ticket.status = 'Processing'
        ticket.meta_info = {'processing_started_at': str(datetime.now())}
        update_ticket_status(ticket, 'Processing started')
        if fdr_class not in stats_dict:
            update_ticket_status(ticket, 'Non-standard fdr threshold values are not supported: {}'.format(fdr_class))
            raise ConvError

        status, data = read_file_as_df(input_file_name, ticket)
        submitted_snps_count = 0
        if status:
            submitted_snps_count = len(data.index)
            if len(data.columns) != 1:
                try:
                    rs_ids, not_found, unique_submitted_snps_count = get_rs_ids_from_vcf(data)
                except ConvError as e:
                    update_ticket_status(ticket,
                                         'Processing failed: the file must contain a single SNP rs-ID on each line or be '
                                         'a valid vcf file, invalid {}'.format(
                                             e.args[0]))
                    raise ConvError
                except:
                    change_status_on_fail = True
                    raise
            else:
                unique_submitted_snps_count = data[0].nunique()
                rs_ids, not_found = get_rs_ids_from_list(data[0].unique())

        if not status or (submitted_snps_count > max_nrows and ticket.user_id != 'adminas'):
            update_ticket_status(ticket,
                                 'Too many SNPs found (>{}). Consider using complete ADASTRA database dump ('
                                 'https://adastra.autosome.ru/downloads) and performing stand-alone enrichment '
                                 'analysis.'.format(max_nrows))
            raise ConvError

        change_status_on_fail = True

        common_header_1 = ['CHROMOSOME', 'POSITION', 'RS_ID', 'REF', 'ALT', 'SEQUENCE']
        common_header_2 = ['PEAK_CALLS', 'MEAN_BAD', 'LOG10_FDR_REF', 'LOG10_FDR_ALT',
                           'EFFECT_SIZE_REF', 'EFFECT_SIZE_ALT']
        common_header_3 = ['GTEX_EQTL', 'EBI', 'PHEWAS', 'FINEMAPPING', 'GRASP', 'CLINVAR', 'GTEX_EQTL_TARGET_GENES']
        cl_header = common_header_1 + ['CELL_TYPE'] + ['CELL_TYPE_GTRD_ID'] + common_header_2 + [
            'SUPPORTING_TFS'] + common_header_3
        tf_header = common_header_1 + ['TRANSCRIPTION_FACTOR'] + ['TF_UNIPROT_AC'] + common_header_2 + \
                    ['MOTIF_LOG_P_REF', 'MOTIF_LOG_P_ALT', 'MOTIF_LOG2_FC', 'MOTIF_POSITION',
                     'MOTIF_ORIENTATION', 'MOTIF_CONCORDANCE', 'SUPPORTING_CELL_TYPES'] + common_header_3

        tf_asb_counts = {}
        tf_sum_counts = {}
        ananastra_service.create_processed_path(ticket_id)
        logger.info('Ticket {}: processing started'.format(ticket_id))
        update_ticket_status(ticket, 'Searching for ASBs of transcription factors (TF-ASBs)')

        tf_path = ananastra_service.get_path_by_ticket_id(ticket_id, path_type='tf', ext='.tsv')

        with open(tf_path, 'w', encoding='utf-8') as out:
            out.write(pack(tf_header))

        for q_tf in divide_query(lambda x: get_tf_query(x, fdr_class), rs_ids):
            with open(tf_path, 'a', encoding='utf-8') as out:
                for tup in q_tf:
                    tf_name = tup[6]
                    tf_asb_counts.setdefault(tf_name, {
                        'name': tf_name,
                        'count': 0
                    })['count'] += 1

                    out.write(pack(process_row(tup, 'TF', tf_header)))

        logger.info('Ticket {}: tf done'.format(ticket_id))
        update_ticket_status(ticket, 'Aggregating TF-ASBs information')

        tf_table = pd.read_table(tf_path, encoding='utf-8', na_values=['None', 'NaN', 'nan'])
        tf_table['LOG10_TOP_FDR'] = tf_table[['LOG10_FDR_REF', 'LOG10_FDR_ALT']].max(axis=1)
        tf_table['IS_EQTL'] = tf_table['GTEX_EQTL_TARGET_GENES'].apply(lambda x: False if pd.isna(x) else True)
        idx = tf_table.groupby(['RS_ID', 'ALT'])['LOG10_TOP_FDR'].transform(max) == tf_table['LOG10_TOP_FDR']
        tf_sum_table = tf_table.loc[idx].copy()
        if len(idx) > 0:
            tf_sum_table['TOP_EFFECT_SIZE'] = tf_sum_table.apply(
                lambda row: row['EFFECT_SIZE_REF'] if row['LOG10_FDR_REF'] >= row['LOG10_FDR_ALT'] else row[
                    'EFFECT_SIZE_ALT'], axis=1)
            tf_sum_table['PREFERRED_ALLELE'] = tf_sum_table.apply(
                lambda row: 'Ref ({})'.format(row['REF']) if row['LOG10_FDR_REF'] >= row[
                    'LOG10_FDR_ALT'] else 'Alt ({})'.format(row['ALT']), axis=1)
            tf_sum_table['MINOR_ALLELE'] = tf_sum_table.apply(
                lambda row: 'Alt ({})'.format(row['ALT']) if row['LOG10_FDR_REF'] >= row[
                    'LOG10_FDR_ALT'] else 'Ref ({})'.format(row['REF']), axis=1)
            tf_table.drop(columns=['LOG10_TOP_FDR'], inplace=True)
            tf_sum_table.drop(columns=['LOG10_FDR_REF', 'LOG10_FDR_ALT', 'EFFECT_SIZE_REF', 'EFFECT_SIZE_ALT'],
                              inplace=True)
            tf_sum_table['IS_EQTL'] = tf_sum_table['GTEX_EQTL_TARGET_GENES'].apply(
                lambda x: False if pd.isna(x) else True)
            tf_sum_table['ALLELES'] = tf_sum_table.apply(
                lambda row: get_alleles(tf_table.loc[tf_table['RS_ID'] == row['RS_ID'], ['REF', 'ALT']]), axis=1)
            tf_sum_table['TF_BINDING_PREFERENCES'] = tf_sum_table.apply(lambda row: get_preferences(
                tf_table.loc[tf_table['RS_ID'] == row['RS_ID'], ['LOG10_FDR_REF', 'LOG10_FDR_ALT']]), axis=1)
            tf_sum_table.drop(columns=['REF', 'ALT'])
            tf_sum_table.to_csv(ananastra_service.get_path_by_ticket_id(ticket_id, 'tf_sum'), sep='\t', index=False)
            tf_table.to_csv(tf_path, sep='\t', index=False)
            tf_sum_counts = tf_sum_table['TRANSCRIPTION_FACTOR'].value_counts().to_dict()
        else:
            tf_sum_table.to_csv(ananastra_service.get_path_by_ticket_id(ticket_id, 'tf_sum'), sep='\t', index=False)

        logger.info('Ticket {}: tf_sum done'.format(ticket_id))
        update_ticket_status(ticket, 'Searching for cell type-ASBs')

        cl_asb_counts = {}
        cl_sum_counts = {}
        cl_path = ananastra_service.get_path_by_ticket_id(ticket_id, path_type='cl', ext='.tsv')

        with open(cl_path, 'w', encoding='utf-8') as out:
            out.write(pack(cl_header))

        for q_cl in divide_query(lambda x: get_cl_query(x, fdr_class), rs_ids):
            with open(cl_path, 'a', encoding='utf-8') as out:
                for tup in q_cl:
                    cl_name = tup[6]
                    cl_asb_counts.setdefault(cl_name, {
                        'name': cl_name,
                        'count': 0
                    })['count'] += 1

                    out.write(pack(process_row(tup, 'CL', cl_header)))

        logger.info('Ticket {}: cl done'.format(ticket_id))
        update_ticket_status(ticket, 'Aggregating CL-ASBs information')

        cl_table = pd.read_table(cl_path, encoding='utf-8', na_values=['None', 'NaN', 'nan'])
        cl_table['LOG10_TOP_FDR'] = cl_table[['LOG10_FDR_REF', 'LOG10_FDR_ALT']].max(axis=1)
        cl_table['IS_EQTL'] = cl_table['GTEX_EQTL_TARGET_GENES'].apply(lambda x: False if pd.isna(x) else True)
        idx = cl_table.groupby(['RS_ID', 'ALT'])['LOG10_TOP_FDR'].transform(max) == cl_table['LOG10_TOP_FDR']
        cl_sum_table = cl_table.loc[idx].copy()
        if len(idx) > 0:
            cl_sum_table['TOP_EFFECT_SIZE'] = cl_sum_table.apply(
                lambda row: row['EFFECT_SIZE_REF'] if row['LOG10_FDR_REF'] >= row['LOG10_FDR_ALT'] else row[
                    'EFFECT_SIZE_ALT'], axis=1)
            cl_sum_table['PREFERRED_ALLELE'] = cl_sum_table.apply(
                lambda row: 'Ref ({})'.format(row['REF']) if row['LOG10_FDR_REF'] >= row[
                    'LOG10_FDR_ALT'] else 'Alt ({})'.format(row['ALT']), axis=1)
            cl_sum_table['MINOR_ALLELE'] = cl_sum_table.apply(
                lambda row: 'Alt ({})'.format(row['ALT']) if row['LOG10_FDR_REF'] >= row[
                    'LOG10_FDR_ALT'] else 'Ref ({})'.format(row['REF']), axis=1)
            cl_table.drop(columns=['LOG10_TOP_FDR'], inplace=True)
            cl_sum_table.drop(columns=['LOG10_FDR_REF', 'LOG10_FDR_ALT', 'EFFECT_SIZE_REF', 'EFFECT_SIZE_ALT'],
                              inplace=True)
            cl_sum_table['ALLELES'] = cl_sum_table.apply(
                lambda row: get_alleles(cl_table.loc[cl_table['RS_ID'] == row['RS_ID'], ['REF', 'ALT']]), axis=1)
            cl_sum_table['TF_BINDING_PREFERENCES'] = cl_sum_table.apply(lambda row: get_preferences(
                cl_table.loc[cl_table['RS_ID'] == row['RS_ID'], ['LOG10_FDR_REF', 'LOG10_FDR_ALT']]), axis=1)
            cl_sum_table.drop(columns=['REF', 'ALT'])
            cl_sum_table.to_csv(ananastra_service.get_path_by_ticket_id(ticket_id, 'cl_sum'), sep='\t', index=False)
            cl_table.to_csv(cl_path, sep='\t', index=False)
            cl_sum_counts = cl_sum_table['CELL_TYPE'].value_counts().to_dict()
        else:
            cl_sum_table.to_csv(ananastra_service.get_path_by_ticket_id(ticket_id, 'cl_sum'), sep='\t', index=False)

        logger.info('Ticket {}: cl_sum done'.format(ticket_id))
        update_ticket_status(ticket, 'Checking the control data of candidate but non-significant ASBs (non-ASBs)')

        tf_sum_counts = [{'name': key, 'count': value} for key, value in tf_sum_counts.items()]
        cl_sum_counts = [{'name': key, 'count': value} for key, value in cl_sum_counts.items()]

        all_rs = len(rs_ids)
        tf_asbs_list = [x for query in divide_query(lambda x: get_asbs_by_level(x, fdr_class, level='TF'), rs_ids) for x
                        in query]
        tf_asbs = len(tf_asbs_list)
        tf_asbs_rs_set = set(x.snp.rs_id for x in tf_asbs_list)
        tf_asbs_rs = len(tf_asbs_rs_set)
        cl_asbs_list = [x for query in divide_query(lambda x: get_asbs_by_level(x, fdr_class, level='CL'), rs_ids) for x
                        in query]
        cl_asbs = len(cl_asbs_list)
        cl_asbs_rs_set = set(x.snp.rs_id for x in cl_asbs_list)
        cl_asbs_rs = len(cl_asbs_rs_set)
        all_asbs_list = [x for query in divide_query(lambda x: get_asbs_by_level(x, fdr_class, level='ALL'), rs_ids) for
                         x in query]
        all_asbs = tf_asbs + cl_asbs
        all_asbs_rs_set = set(x.rs_id for x in all_asbs_list)
        all_asbs_rs = len(all_asbs_rs_set)

        logger.info('Ticket {}: query count asb done'.format(ticket_id))
        update_ticket_status(ticket, 'Checking the control data of candidate but non-significant ASBs (non-ASBs)')

        fdr_class_neg = '0.25'
        # tf_non_negative_candidates_rs_set = [x.rs_id for query in divide_query(lambda x: get_candidates_by_level(x, fdr_class_neg, level='TF', alternative='greater', rs=True), rs_ids) for x in query]
        # cl_non_negative_candidates_rs_set = [x.rs_id for query in divide_query(lambda x: get_candidates_by_level(x, fdr_class_neg, level='CL', alternative='greater', rs=True), rs_ids) for x in query]
        # all_non_negative_candidates_rs_set = [x.rs_id for query in divide_query(lambda x: get_all_candidates(x, fdr_class_neg, alternative='greater', rs=True), rs_ids) for x in query]
        all_non_negative_candidates_rs_set = set(x.rs_id for query in divide_query(
            lambda x: get_candidates_by_level(x, fdr_class_neg, level='ALL', alternative='greater', rs=True, ), rs_ids)
                                                 for x in query)
        all_non_negative_candidates_rs = len(all_non_negative_candidates_rs_set)
        logger.info('Ticket {}: non-negative candidates done'.format(ticket_id))

        tf_negatives_list = [x for query in
                             divide_query(lambda x: get_candidates_by_level(x, fdr_class_neg, level='TF'), rs_ids) for x
                             in query]
        tf_negatives = len(tf_negatives_list)
        tf_negatives_rs_set = set(x.rs_id for query in divide_query(
            lambda x: get_candidates_by_level(x, fdr_class_neg, level='TF', rs=True), rs_ids)
                                  for x in query)
        tf_negatives_rs = len(tf_negatives_rs_set)
        cl_negatives_list = [x for query in
                             divide_query(lambda x: get_candidates_by_level(x, fdr_class_neg, level='CL'), rs_ids) for x
                             in query]
        cl_negatives = len(cl_negatives_list)
        cl_negatives_rs_set = set(x.rs_id for query in divide_query(
            lambda x: get_candidates_by_level(x, fdr_class_neg, level='CL', rs=True), rs_ids)
                                  for x in query)
        cl_negatives_rs = len(cl_negatives_rs_set)
        all_negatives_list = [x for query in divide_query(
            lambda x: get_candidates_by_level(x, fdr_class_neg, level='ALL', rs=False, mode='all'), rs_ids) for x in
                              query]
        all_negatives = tf_negatives + cl_negatives
        all_negatives_rs_set = set(x.rs_id for query in divide_query(
            lambda x: get_candidates_by_level(x, fdr_class_neg, level='ALL', rs=True), rs_ids)
                                   for x in query)
        all_negatives_rs = len(all_negatives_rs_set)

        undefined = all_non_negative_candidates_rs - all_asbs_rs

        logger.info('Ticket {}: query count candidates done'.format(ticket_id))
        update_ticket_status(ticket, 'Aggregating SNP info by rs id')

        rs_set_dict = {
            'asb': {
                'tf': tf_asbs_rs_set,
                'cl': cl_asbs_rs_set,
                'all': all_asbs_rs_set,
            },
            'neg': {
                'tf': tf_negatives_rs_set,
                'cl': cl_negatives_rs_set,
                'all': all_negatives_rs_set,
            },
            'non-neg': {
                'all': all_non_negative_candidates_rs_set,
            }
        }

        header = ['CHROMOSOME', 'POSITION', 'RS_ID', 'REF', 'ALT',
                  'TRANSCRIPTION_FACTOR', 'TF_UNIPROT_AC', 'PEAK_CALLS (TF)', 'MEAN_BAD (TF)', 'MOTIF_LOG_P_REF (TF)',
                  'MOTIF_LOG_P_ALT (TF)', 'MOTIF_LOG2_FC (TF)', 'MOTIF_POSITION (TF)',
                  'MOTIF_ORIENTATION (TF)', 'MOTIF_CONCORDANCE (TF)', 'SUPPORTING_CELL_TYPES (TF)',
                  'LOG10_TOP_FDR (TF)', 'TOP_EFFECT_SIZE (TF)', 'PREFERRED_ALLELE (TF)',
                  'MINOR_ALLELE (TF)', 'ALLELES (TF)',
                  'TF_BINDING_PREFERENCES (TF)',
                  'CELL_TYPE', 'CELL_TYPE_GTRD_ID', 'PEAK_CALLS (CL)', 'MEAN_BAD (CL)', 'SUPPORTING_TFS (CL)',
                  'LOG10_TOP_FDR (CL)', 'TOP_EFFECT_SIZE (CL)', 'PREFERRED_ALLELE (CL)',
                  'MINOR_ALLELE (CL)', 'ALLELES (CL)',
                  'TF_BINDING_PREFERENCES (CL)',
                  'GTEX_EQTL_TARGET_GENES', 'IS_EQTL',
                  'GTEX_EQTL', 'EBI', 'PHEWAS', 'FINEMAPPING', 'GRASP', 'CLINVAR',
                  ]

        attribute_tuples = {tag: ('CL', tag[:-5]) if tag.endswith(' (CL)') else
        ('TF', tag[:-5]) if tag.endswith(' (TF)') else
        ('TF', tag) if tag in ('TRANSCRIPTION_FACTOR', 'TF_UNIPROT_AC') else
        ('CL', tag) if tag in ('CELL_TYPE', 'CELL_TYPE_GTRD_ID') else
        ('ALL', tag)
                            for tag in header
                            }

        list_of_rows = []
        for rs_id in rs_ids:
            row = {'SNP_ID': 'rs{}'.format(rs_id)}
            all_label = get_snp_class_label(rs_id, rs_set_dict, not_found=not_found, ag_class='all')
            if all_label == 'NOT_IN_ADASTRA':
                tf_label = ''
                cl_label = ''
            else:
                tf_label = get_snp_class_label(rs_id, rs_set_dict, not_found=not_found, ag_class='tf')
                cl_label = get_snp_class_label(rs_id, rs_set_dict, not_found=not_found, ag_class='cl')
            row.update({
                'ASB_STATUS': all_label,
                'TF_ASB_STATUS': tf_label,
                'CL_ASB_STATUS': cl_label,
            })
            row.update({x: '' for x in header})
            if tf_label == 'ASB':
                tf_row = tf_sum_table[tf_sum_table['RS_ID'] == 'rs{}'.format(rs_id)]
                for tag, (ag_class, table_tag) in attribute_tuples.items():
                    if ag_class in ('ALL', 'TF'):
                        row[tag] = tf_row[table_tag].values[0]
            if cl_label == 'ASB':
                cl_row = cl_sum_table[cl_sum_table['RS_ID'] == 'rs{}'.format(rs_id)]
                for tag, (ag_class, table_tag) in attribute_tuples.items():
                    if ag_class in ('ALL', 'CL'):
                        row[tag] = cl_row[table_tag].values[0]
            list_of_rows.append(row)

        for snp_id in not_found:
            row = {
                'SNP_ID': snp_id
            }
            row.update({
                'ASB_STATUS': 'NOT_IN_ADASTRA',
                'TF_ASB_STATUS': '',
                'CL_ASB_STATUS': '',
            })
            row.update({x: '' for x in header})
            list_of_rows.append(row)

        not_found_ids = [row['SNP_ID'] for row in list_of_rows if row['ASB_STATUS'] == 'NOT_IN_ADASTRA']

        all_table = pd.DataFrame(list_of_rows,
                                 columns=['SNP_ID', 'ASB_STATUS', 'TF_ASB_STATUS', 'CL_ASB_STATUS'] + header)
        all_table.to_csv(ananastra_service.get_path_by_ticket_id(ticket_id, 'all'), sep='\t', index=False)

        with open(ananastra_service.get_path_by_ticket_id(ticket_id, 'not_found'), 'w') as f:
            f.write('\n'.join(not_found_ids))

        logger.info('Ticket {}: rs report done'.format(ticket_id))
        update_ticket_status(ticket, 'Checking the control data of candidate but non-significant ASBs (non-ASBs)')

        if background == 'LOCAL' or background.startswith('LD'):
            if background == 'LOCAL':
                window_size = 10 ** 6
                local_intervals = snvs_to_hash(*zip(*(session.query(PositionHash.position_hash)
                                                      .join(CandidateRS, PositionHash.rs_id == CandidateRS.rs_id)
                                                      .filter(CandidateRS.rs_id.in_(rs_ids)))), w_size=window_size)
                args = {'intervals': local_intervals}

            elif background.startswith('LD'):
                if background == 'LD-EUR':
                    column = 'ld_eur'
                elif background == 'LD-ASN':
                    column = 'ld_asn'
                elif background == 'LD-AFR':
                    column = 'ld_afr'
                else:
                    raise ValueError

                lds = set(
                    *zip(*(session.query(getattr(LDIslandsInfo, column).distinct()).join(
                        CandidateRS, CandidateRS.rs_id == LDIslandsInfo.rs_id
                    ).filter(
                        CandidateRS.rs_id.in_(rs_ids),

                    )
                    ))) - {None}

                args = {'lds': lds, 'ld_type': column}
            else:
                raise ValueError

            expected_tf_asbs_rs = get_candidates_by_level(**args, fdr=fdr_class, alternative='greater', level='TF',
                                                          rs=True, mode='count') - tf_asbs_rs
            expected_cl_asbs_rs = get_candidates_by_level(**args, fdr=fdr_class, alternative='greater', level='CL',
                                                          rs=True, mode='count') - cl_asbs_rs
            expected_all_asbs_rs = get_candidates_by_level(**args, fdr=fdr_class, alternative='greater', level='ALL',
                                                           rs=True, mode='count') - all_asbs_rs

            expected_tf_negatives_rs = get_candidates_by_level(**args, fdr=fdr_class_neg, level='TF', rs=True,
                                                               mode='count') - tf_negatives_rs
            expected_cl_negatives_rs = get_candidates_by_level(**args, fdr=fdr_class_neg, level='CL', rs=True,
                                                               mode='count') - cl_negatives_rs
            expected_all_negatives_rs = get_candidates_by_level(**args, fdr=fdr_class_neg, level='ALL', rs=True,
                                                                mode='count') - all_negatives_rs

        elif background == 'WG':
            sd = stats_dict[fdr_class]
            cand_sd = stats_dict[fdr_class_neg]

            expected_tf_asbs_rs = sd['expected_tf_asbs_rs'] - tf_asbs_rs
            expected_cl_asbs_rs = sd['expected_cl_asbs_rs'] - cl_asbs_rs
            expected_all_asbs_rs = sd['expected_all_asbs_rs'] - all_asbs_rs
            expected_tf_negatives_rs = stats_dict['1']['total_tf_candidates_rs'] - cand_sd[
                'expected_tf_asbs_rs'] - tf_negatives_rs
            expected_cl_negatives_rs = stats_dict['1']['total_cl_candidates_rs'] - cand_sd[
                'expected_cl_asbs_rs'] - cl_negatives_rs
            expected_all_negatives_rs = stats_dict['1']['total_all_candidates_rs'] - cand_sd[
                'expected_all_asbs_rs'] - all_negatives_rs
        else:
            raise ValueError

        update_ticket_status(ticket, 'Performing statistical analysis')

        tf_odds_rs, tf_p_rs = fisher_exact(
            ((tf_asbs_rs, tf_negatives_rs), (expected_tf_asbs_rs, expected_tf_negatives_rs)), alternative='greater')
        cl_odds_rs, cl_p_rs = fisher_exact(
            ((cl_asbs_rs, cl_negatives_rs), (expected_cl_asbs_rs, expected_cl_negatives_rs)), alternative='greater')
        all_odds_rs, all_p_rs = fisher_exact(
            ((all_asbs_rs, all_negatives_rs), (expected_all_asbs_rs, expected_all_negatives_rs)), alternative='greater')

        tf_asb_data = []
        cl_asb_data = []
        for level in 'TF', 'CL':
            if level == 'TF':
                logger.info('Ticket {}: tests done'.format(ticket_id))
                update_ticket_status(ticket, 'Testing the enrichment of ASBs of individual TFs')
            else:
                logger.info('Ticket {}: tf tests done'.format(ticket_id))
                update_ticket_status(ticket, 'Testing the enrichment of ASBs of individual cell types')

            p_list = []
            asb_counts = {'TF': tf_asb_counts, 'CL': cl_asb_counts}[level]
            asbs_list = {'TF': tf_asbs_list, 'CL': cl_asbs_list}[level]
            negatives_list = {'TF': tf_negatives_list, 'CL': cl_negatives_list}[level]
            ag_stats_dict = {'TF': tf_stats_dict, 'CL': cl_stats_dict}[level]
            AgClass = {'TF': TranscriptionFactor, 'CL': CellLine}[level]
            id_attr = '{}_id'.format(level.lower())
            asb_data = {'TF': tf_asb_data, 'CL': cl_asb_data}[level]
            for name in asb_counts.keys():
                ag_id = getattr(AgClass.query.filter_by(name=name).one(), id_attr)
                asbs = asb_counts[name]['count']
                asbs_rs = len(set(x.snp.rs_id for x in asbs_list if getattr(x, id_attr) == ag_id))
                negatives = len([cand for cand in negatives_list if cand.ag_id == ag_id])
                negatives_rs = len(set(cand.rs_id for cand in negatives_list if cand.ag_id == ag_id))
                expected_asbs = ag_stats_dict[str(ag_id)][fdr_class][
                                    'expected_{}_asbs'.format(level.lower())] - asbs
                expected_asbs_rs = ag_stats_dict[str(ag_id)][fdr_class][
                                       'expected_{}_asbs_rs'.format(level.lower())] - asbs_rs
                expected_negatives_rs = ag_stats_dict[str(ag_id)]['1']['total_{}_candidates_rs'.format(level.lower())] - \
                                        ag_stats_dict[str(ag_id)][fdr_class_neg][
                                            'expected_{}_asbs_rs'.format(level.lower())] - negatives_rs
                odds, p = fisher_exact(((asbs_rs, negatives_rs), (expected_asbs_rs, expected_negatives_rs)),
                                       alternative='greater')
                p_list.append(p)
                asb_data.append({
                    'name': name,
                    'asbs': asbs,
                    'asbs_rs': asbs_rs,
                    'negatives': negatives,
                    'negatives_rs': negatives_rs,
                    'expected_asbs': expected_asbs,
                    'expected_asbs_rs': expected_asbs_rs,
                    'expected_negatives_rs': expected_negatives_rs,
                    'odds': odds,
                    'log10_p_value': -np.log10(p),
                    'log10_fdr': 0,
                })
            if len(p_list) == 0:
                pass
            else:
                _, ag_fdr, _, _ = multipletests(p_list, alpha=0.05, method='fdr_bh')
                for sig, logfdr in zip(asb_data, ag_fdr):
                    sig['log10_fdr'] = np.nan if np.isnan(logfdr) else -np.log10(logfdr)

        logger.info('Ticket {}: cl tests done'.format(ticket_id))
        update_ticket_status(ticket, 'Testing the enrichment of ASBs at individual chromosomes')

        chr_p_list = []
        chr_tf_p_list = []
        chr_cl_p_list = []
        chr_asb_data = []
        for chromosome in list(
                set(x.chromosome for x in all_asbs_list) | set(x.chromosome for x in all_negatives_list)):
            counts = get_asb_counts_by_chromosome(all_asbs_list, all_negatives_list, chromosome)
            expected_chr_asbs_rs = chr_stats_dict[chromosome][fdr_class]['expected_all_asbs_rs'] - counts['asbs_rs']
            expected_chr_negatives_rs = chr_stats_dict[chromosome]['1']['total_all_candidates_rs'] - \
                                        chr_stats_dict[chromosome][fdr_class_neg]['expected_all_asbs_rs'] - \
                                        counts['negatives_rs']
            odds, p = get_p_and_odds_by_chromosome(counts, expected_chr_asbs_rs, expected_chr_negatives_rs)
            chr_p_list.append(p)

            tf_counts = get_asb_counts_by_chromosome(tf_asbs_list, tf_negatives_list, chromosome, ag_snp=True)
            expected_chr_tf_asbs_rs = chr_stats_dict[chromosome][fdr_class]['expected_tf_asbs_rs'] - counts['asbs_rs']
            expected_chr_tf_negatives_rs = chr_stats_dict[chromosome]['1']['total_tf_candidates_rs'] - \
                                           chr_stats_dict[chromosome][fdr_class_neg]['expected_tf_asbs_rs'] - \
                                           counts['negatives_rs']
            tf_odds, tf_p = get_p_and_odds_by_chromosome(tf_counts, expected_chr_tf_asbs_rs,
                                                         expected_chr_tf_negatives_rs)
            chr_tf_p_list.append(tf_p)

            cl_counts = get_asb_counts_by_chromosome(cl_asbs_list, cl_negatives_list, chromosome, ag_snp=True)
            expected_chr_cl_asbs_rs = chr_stats_dict[chromosome][fdr_class]['expected_cl_asbs_rs'] - counts['asbs_rs']
            expected_chr_cl_negatives_rs = chr_stats_dict[chromosome]['1']['total_cl_candidates_rs'] - \
                                           chr_stats_dict[chromosome][fdr_class_neg]['expected_cl_asbs_rs'] - \
                                           counts['negatives_rs']
            cl_odds, cl_p = get_p_and_odds_by_chromosome(cl_counts, expected_chr_cl_asbs_rs,
                                                         expected_chr_cl_negatives_rs)
            chr_cl_p_list.append(cl_p)

            chr_asb_data.append({
                'name': chromosome,

                'asbs': counts['asbs'],
                'asbs_rs': counts['asbs_rs'],
                'negatives': counts['negatives'],
                'negatives_rs': counts['negatives_rs'],
                'expected_asbs_rs': expected_chr_asbs_rs,
                'expected_negatives_rs': expected_chr_negatives_rs,
                'odds': odds,
                'log10_p_value': -np.log10(p),
                'log10_fdr': 0,

                'tf_asbs': tf_counts['asbs'],
                'tf_asbs_rs': tf_counts['asbs_rs'],
                'tf_negatives': tf_counts['negatives'],
                'tf_negatives_rs': tf_counts['negatives_rs'],
                'expected_tf_asbs_rs': expected_chr_tf_asbs_rs,
                'expected_tf_negatives_rs': expected_chr_tf_negatives_rs,
                'tf_odds': tf_odds,
                'tf_log10_p_value': -np.log10(tf_p),

                'cl_asbs': cl_counts['asbs'],
                'cl_asbs_rs': cl_counts['asbs_rs'],
                'cl_negatives': cl_counts['negatives'],
                'cl_negatives_rs': cl_counts['negatives_rs'],
                'expected_cl_asbs_rs': expected_chr_cl_asbs_rs,
                'expected_cl_negatives_rs': expected_chr_cl_negatives_rs,
                'cl_odds': cl_odds,
                'cl_log10_p_value': -np.log10(cl_p),
            })
        if len(chr_p_list) == 0:
            chr_fdr = []
        else:
            _, chr_fdr, _, _ = multipletests(chr_p_list, alpha=0.05, method='fdr_bh')
        for sig, logfdr in zip(chr_asb_data, chr_fdr):
            sig['log10_fdr'] = np.nan if np.isnan(logfdr) else -np.log10(logfdr)

        chr_p_rs = logit_combine_p_values(chr_p_list)
        chr_tf_p_rs = logit_combine_p_values(chr_tf_p_list)
        chr_cl_p_rs = logit_combine_p_values(chr_cl_p_list)

        logger.info('Ticket {}: chromosome tests done'.format(ticket_id))
        update_ticket_status(ticket, 'Finalizing the report')

        tf_asb_data = sorted(tf_asb_data, key=lambda x: (x['log10_fdr'], x['log10_p_value'], x['odds']), reverse=True)
        cl_asb_data = sorted(cl_asb_data, key=lambda x: (x['log10_fdr'], x['log10_p_value'], x['odds']), reverse=True)
        chr_asb_data = sorted(chr_asb_data, key=lambda x: (x['log10_fdr'], x['log10_p_value'], x['odds']), reverse=True)

        ticket.status = 'Processed'
        meta_info = dict(ticket.meta_info)
        meta_info.update({
            'processing_time': str(datetime.now() - processing_start_time),
            'all_rs': all_rs,
            'submitted_snps_count': submitted_snps_count,
            'unique_submitted_snps_count': unique_submitted_snps_count,
            'undefined_rs': undefined,
            'tf': {
                'asbs': tf_asbs,
                'asbs_rs': tf_asbs_rs,
                'negatives': tf_negatives,
                'negatives_rs': tf_negatives_rs,
                'expected_asbs_rs': expected_tf_asbs_rs,
                'expected_negatives_rs': expected_tf_negatives_rs,
                'odds_rs': marshal_inf(tf_odds_rs),
                'log10_p_value_rs': marshal_logp(tf_p_rs),
                'asb_counts': modify_counts(tf_asb_data, top=False),
                'asb_counts_top': modify_counts(tf_asb_data, tf_sum_counts, top=True),
                'asb_data': marshal_data(tf_asb_data),
            },
            'cl': {
                'asbs': cl_asbs,
                'asbs_rs': cl_asbs_rs,
                'negatives': cl_negatives,
                'negatives_rs': cl_negatives_rs,
                'expected_asbs_rs': expected_cl_asbs_rs,
                'expected_negatives_rs': expected_cl_negatives_rs,
                'odds_rs': marshal_inf(cl_odds_rs),
                'log10_p_value_rs': marshal_logp(cl_p_rs),
                'asb_counts': modify_counts(cl_asb_data, top=False),
                'asb_counts_top': modify_counts(cl_asb_data, cl_sum_counts, top=True),
                'asb_data': marshal_data(cl_asb_data),
            },
            'all': {
                'asbs': all_asbs,
                'asbs_rs': all_asbs_rs,
                'negatives': all_negatives,
                'negatives_rs': all_negatives_rs,
                'expected_asbs_rs': expected_all_asbs_rs,
                'expected_negatives_rs': expected_all_negatives_rs,
                'odds_rs': marshal_inf(all_odds_rs),
                'log10_p_value_rs': marshal_logp(all_p_rs),
            },
            'chr': {
                'log10_p_value_rs': marshal_logp(chr_p_rs),
                'tf_log10_p_value_rs': marshal_logp(chr_tf_p_rs),
                'cl_log10_p_value_rs': marshal_logp(chr_cl_p_rs),
                'asb_data': marshal_chr_data(chr_asb_data),
            }
        })
    except Exception as e:
        if not isinstance(e, ConvError):
            logger.error(e, exc_info=True)
        ticket.status = 'Failed'
        if change_status_on_fail:
            update_ticket_status(ticket, 'Processing failed while "{}".'.format(ticket.meta_info['status_details']))
        session.commit()
        return

    ticket.meta_info = meta_info
    logger.info('Ticket {}: ticket info changed'.format(ticket_id))
    session.commit()

    logger.info('Ticket {}: session commited'.format(ticket_id))
    update_ticket_status(ticket, 'Processing finished')
    session.commit()
