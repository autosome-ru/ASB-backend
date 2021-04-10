from svgutils import transform
import numpy as np
from math import lgamma
import os

from ASB_app import *
from ASB_app.releases import current_release

session = current_release.session
SNP = current_release.SNP
TranscriptionFactor = current_release.TranscriptionFactor
TranscriptionFactorSNP = current_release.TranscriptionFactorSNP


def get_KDIC(counts, N):
    return 1 / (N * np.log(0.25)) * (
                lgamma(N + 1) + sum([x * np.log(0.25) - lgamma(x + 1) for x in counts]))


def get_heights(pcm_file, mode='freq'):
    heights = []
    lines = []
    with open(pcm_file) as f:
        for line in f:
            if line.startswith('>'):
                continue
            lines.append(list(map(float, line.strip('\n').split('\t'))))

    m = len(lines)
    for counts in lines:
        N = sum(counts)
        if mode == 'freq':
            heights.append(sorted(list(zip(letters, [x / N for x in counts])), key=lambda x: x[1], reverse=True))
        elif mode == 'KDIC':
            KDIC = get_KDIC(counts, N)
            heights.append(sorted(list(zip(letters, [(x / N) * KDIC for x in counts])), key=lambda x: x[1], reverse=True))

    return m, heights


def place_drawing_on_svg(figure, drawing_svg, x, y, h, w, self_w, self_h):
    drawing_object = transform.fromfile(drawing_svg)
    drawing_root = drawing_object.getroot()
    drawing_root.scale_xy(w/self_w, h/self_h)
    drawing_root.moveto(x, y)
    figure.append(drawing_root)


def place_letter_on_svg(figure, letter_svg, x, y, h, w):
    # 13.229 and 26.458 are letter svg view box w and h
    place_drawing_on_svg(figure, letter_svg, x, y, h, w, 13.229, 26.458)


def place_dna_on_svg(figure, dna_svg, x, y, h, w):
    place_drawing_on_svg(figure, dna_svg, x, y, h, w, 365, 100)


def place_hill_on_svg(figure, hill_svg, x, y, h, w):
    place_drawing_on_svg(figure, hill_svg, x, y, h, w, 2420.548, 479.477)


def place_concordance_on_svg(figure, concordance, x, y, h, w):
    concordance_svg = os.path.expanduser('~/letters/{}.svg'.format(concordance))
    place_drawing_on_svg(figure, concordance_svg, x, y, h, w, conc_self_width, concordance_sizes[concordance])


def renorm(position):
    letters, heights = zip(*position)
    total_height = sum(heights)
    new_total_height = 0
    new_heights = []
    for height in heights:
        if height < visible_cut_tr * total_height:
            new_heights.append(0)
        else:
            new_total_height += height
            new_heights.append(height)
    new_heights = [x * new_total_height / total_height for x in new_heights]
    return zip(letters, new_heights)


def get_scientific_text(n):
    convert_exp = {
        '0': '⁰',
        '1': '¹',
        '2': '²',
        '3': '³',
        '4': '⁴',
        '5': '⁵',
        '6': '⁶',
        '7': '⁷',
        '8': '⁸',
        '9': '⁹',
        '-': '⁻',
    }
    base, exp = '{:.2e}'.format(n).split('e')
    return base + '·10' + ''.join([convert_exp[x] for x in str(int(exp))])


if __name__ == '__main__':
    name_dict = dict()

    for file_name in os.listdir(os.path.expanduser('~/pcm/')):
        tf_name = file_name.split('.')[0]
        if tf_name in name_dict:
            print('Бросаем корабль!')
        name_dict[tf_name] = file_name

    letters = ['A', 'C', 'G', 'T']

    unit_width = 300
    unit_height = 600

    motif_gap = 0.1
    snp_gap = 0.1
    text_h = 0.5
    text_width = 0.85
    snp_text_h = 0.8
    snp_text_width = 1
    hill_sum_height = 2
    hill_width = 'to be changed'
    pseudocount = 0.15
    add_letters = 6
    p_value_text_h = 0.3
    strands_h = 0.5
    hill_gap = 0.2
    dna_h = 1/3*2/3*0.8

    concordance_indent = 1
    concordance_scale = 3000/(3351.318*(unit_height/300))
    conc_self_width = 331.204

    visible_cut_tr = 0.05

    indent = snp_text_h + snp_gap / 2

    full_gap = motif_gap + indent

    concordance_sizes = {
        'Weak Concordant': 3351.318,
        'Weak Discordant': 3226.231,
        'Concordant': 2183.019,
        'Discordant': 2044.665,
    }

    concordance_h = {concordance: concordance_sizes[concordance] * concordance_scale for concordance in concordance_sizes}
    concordance_w = conc_self_width * concordance_scale

    letter_svgs = {
        'A': os.path.expanduser('~/letters/lettarA_path.svg'),
        'C': os.path.expanduser('~/letters/lettarC_path.svg'),
        'G': os.path.expanduser('~/letters/lettarG_path.svg'),
        'T': os.path.expanduser('~/letters/lettarT_path.svg'),
    }

    black_letter_svgs = {
        'A': os.path.expanduser('~/letters/lettarA_path_black.svg'),
        'C': os.path.expanduser('~/letters/lettarC_path_black.svg'),
        'G': os.path.expanduser('~/letters/lettarG_path_black.svg'),
        'T': os.path.expanduser('~/letters/lettarT_path_black.svg'),
    }

    hill_svgs = {
        'A': os.path.expanduser('~/letters/hill_A.svg'),
        'C': os.path.expanduser('~/letters/hill_C.svg'),
        'G': os.path.expanduser('~/letters/hill_G.svg'),
        'T': os.path.expanduser('~/letters/hill_T.svg'),
    }

    get_revcomp = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
    }

    for tf_snp, snp, tf in session.query(TranscriptionFactorSNP, SNP, TranscriptionFactor).filter(TranscriptionFactorSNP.motif_concordance.isnot(None) &
                                                                                                  (TranscriptionFactorSNP.motif_concordance != 'No Hit')).join(
        SNP,
        (SNP.chromosome == TranscriptionFactorSNP.chromosome) &
        (SNP.position == TranscriptionFactorSNP.position) &
        (SNP.alt == TranscriptionFactorSNP.alt)
    ).join(
        TranscriptionFactor,
        TranscriptionFactor.tf_id == TranscriptionFactorSNP.tf_id
    ):
        for draw_revcomp in True, False:
            pcm_filename = name_dict[tf.name]

            pcm_path = os.path.expanduser('~/pcm/{}'.format(pcm_filename))

            context = ' ' * 20 + (''.join([get_revcomp[x] for x in snp.context[::-1]]) if draw_revcomp else snp.context) + ' ' * 20
            print(context)
            alt = get_revcomp[snp.alt] if draw_revcomp else snp.alt
            pos_in_motif = tf_snp.motif_position
            motif_pref = 1/np.power(10, tf_snp.motif_log_p_ref)
            motif_palt = 1/np.power(10, tf_snp.motif_log_p_alt)
            asb_is_ref = tf_snp.log_p_value_ref > tf_snp.log_p_value_alt
            revcomp = tf_snp.motif_orientation == draw_revcomp

            ef = tf_snp.es_ref if asb_is_ref else tf_snp.es_alt
            fdr = 1/np.power(10, tf_snp.log_p_value_ref) if asb_is_ref else 1/np.power(10, tf_snp.log_p_value_alt)

            m, heights = get_heights(pcm_path, mode='KDIC')
            print(full_gap, text_h, indent)
            fig_x = (m + add_letters + add_letters + concordance_indent) * unit_width
            fig_y = unit_height * (1 + full_gap + text_h + indent + hill_sum_height + strands_h + hill_gap + 0.35)
            fig = transform.SVGFigure("{}px".format(fig_x), "{}px".format(fig_y))
            txt_gen = transform.TextElement(fig_x - 0.1 * unit_width, fig_y - 0.1 * unit_height,
                                            'ADASTRA v{}'.format(current_release.full_version), size=str(p_value_text_h * unit_height) + 'px',
                                            anchor='end', color='#a8a8a8', font='PT Sans, Arial')
            txt_snp = transform.TextElement(fig_x - 0.1 * unit_width, fig_y - (2 * p_value_text_h + 0.1) * unit_height,
                                            'rs{}'.format(snp.rs_id), size=str(p_value_text_h * unit_height) + 'px',
                                            anchor='end', color='#a8a8a8', font='PT Sans, Arial')
            txt_tf = transform.TextElement(fig_x - 0.1 * unit_width, fig_y - (p_value_text_h + 0.1) * unit_height,
                                           tf.name, size=str(p_value_text_h * unit_height) + 'px',
                                           anchor='end', color='#a8a8a8', font='PT Sans, Arial')
            txt_strand = transform.TextElement(0.1 * unit_width, fig_y - (0.1) * unit_height,
                                           '({}) strand'.format('-' if draw_revcomp else '+'), size=str(p_value_text_h * unit_height) + 'px',
                                           anchor='start', color='#a8a8a8', font='PT Sans, Arial')
            fig.append([txt_gen, txt_snp, txt_tf, txt_strand])

            pos_in_motif = m - pos_in_motif - 1 if revcomp else pos_in_motif
            motif_context = context[44 - add_letters - pos_in_motif: 44 - pos_in_motif + m + add_letters]
            pos_in_motif += add_letters

            hill_height = min(min(ef, 3) / 3, 1 - pseudocount)
            hill_height = (hill_height + 1) / 2
            ref_height = hill_sum_height * hill_height
            alt_height = hill_sum_height * (1 - hill_height)

            hill_width = m

            print(motif_context, pos_in_motif)

            place_letter_on_svg(fig, os.path.expanduser('~/letters/rect.svg'), (pos_in_motif + concordance_indent) * unit_width, 0, (1 + full_gap + text_h/2 + snp_gap/2 + snp_text_h) * unit_height, unit_width)

            for pos, pack in enumerate(heights[::-1] if revcomp else heights):
                pos += add_letters + concordance_indent
                current_height = 0
                for letter, height in renorm(pack):
                    # Draw letter with offset of pos*unit_width, current_height*unit_height and height of height*unit_height
                    place_letter_on_svg(fig, letter_svgs[get_revcomp[letter] if revcomp else letter], pos*unit_width, (1-current_height - height)*unit_height, height*unit_height, unit_width)
                    current_height += height

            for pos, letter in enumerate(motif_context):
                if letter == ' ':
                    continue
                if pos != pos_in_motif:
                    place_letter_on_svg(fig, black_letter_svgs[letter.upper()], (pos + (1 - text_width) / 2 + concordance_indent) * unit_width, (1 + full_gap) * unit_height, text_h * unit_height, text_width * unit_width)
                else:
                    if not asb_is_ref:
                        ref_height, alt_height = alt_height, ref_height
                    place_hill_on_svg(fig, hill_svgs[letter], (concordance_indent + add_letters)*unit_width, (1 + full_gap + text_h + indent + hill_gap) * unit_height, ref_height * unit_height, hill_width * unit_width)
                    place_hill_on_svg(fig, hill_svgs[alt], (concordance_indent + add_letters)*unit_width, (1 + full_gap + strands_h + alt_height + ref_height + text_h + indent + hill_gap) * unit_height, -alt_height * unit_height, hill_width * unit_width)
                    place_letter_on_svg(fig, letter_svgs[letter], (pos + (1 - snp_text_width) / 2 + concordance_indent) * unit_width, (1 + full_gap + text_h/2 - snp_text_h - snp_gap/2) * unit_height, snp_text_h * unit_height, snp_text_width * unit_width)
                    place_letter_on_svg(fig, letter_svgs[alt], (pos + (1 - snp_text_width) / 2 + concordance_indent) * unit_width, (1 + full_gap + text_h/2 + snp_gap/2) * unit_height, snp_text_h * unit_height, snp_text_width * unit_width)

            for pos in range(m-1, -1, -1):
                pos += add_letters + concordance_indent
                place_dna_on_svg(fig, os.path.expanduser('~/letters/dna_grey.svg'), (pos - 13/120)*unit_width, (1 + full_gap + text_h + indent + hill_gap + ref_height) * unit_height, dna_h*unit_height, unit_width*(1 + 13/60))
                place_dna_on_svg(fig, os.path.expanduser('~/letters/dna_grey.svg'), (pos - 13/120)*unit_width, (1 + full_gap + text_h + indent + hill_gap + ref_height + strands_h - dna_h) * unit_height, dna_h*unit_height, unit_width*(1 + 13/60))

            text_x = pos_in_motif + 1 + concordance_indent
            letter_text = pos_in_motif + 1 - 0.05
            txt_ref = transform.TextElement(text_x * unit_width, (1.24 + (motif_gap + indent)/2 + p_value_text_h/2) * unit_height, 'P-value: ' + get_scientific_text(motif_pref), size=str(p_value_text_h*unit_height) + 'px', color='#000000de', font='PT Sans, Arial')
            motif_word_ref = transform.TextElement(text_x * unit_width, (1.24 + (motif_gap + indent)/2 - p_value_text_h/2) * unit_height, 'Motif', size=str(p_value_text_h*unit_height) + 'px', color='#000000de', font='PT Sans, Arial')
            txt_alt = transform.TextElement(text_x * unit_width, (1 + motif_gap + indent*3/2 + text_h + p_value_text_h/2) * unit_height, 'P-value: ' + get_scientific_text(motif_palt), size=str(p_value_text_h*unit_height) + 'px', color='#000000de', font='PT Sans, Arial')
            motif_word_alt = transform.TextElement(text_x * unit_width, (1 + motif_gap + indent*3/2 + text_h - p_value_text_h/2) * unit_height, 'Motif', size=str(p_value_text_h*unit_height) + 'px', color='#000000de', font='PT Sans, Arial')
            txt_ref_letter = transform.TextElement((letter_text) * unit_width,
                                            (1.24 + (motif_gap + indent)/2 - p_value_text_h/2) * unit_height,
                                            'Ref',
                                            size=str(p_value_text_h * unit_height) + 'px', color='#000000de',
                                            font='PT Sans, Arial', anchor='end')
            txt_alt_letter = transform.TextElement((letter_text) * unit_width, (1 + motif_gap + indent*3/2 + text_h + p_value_text_h/2) * unit_height, 'Alt',
                                                   size=str(p_value_text_h * unit_height) + 'px', color='#000000de',
                                                   font='PT Sans, Arial', anchor='end')

            fig.append([txt_ref, txt_alt, motif_word_ref, motif_word_alt, txt_ref_letter, txt_alt_letter])

            if asb_is_ref:
                ef_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height*1 - p_value_text_h/2
                fdr_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height*(1-1/3) - p_value_text_h/2
            else:
                ef_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height + strands_h + alt_height*1/3 - p_value_text_h/2
                fdr_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height + strands_h + alt_height*2/3 - p_value_text_h/2

            txt_ef = transform.TextElement(((m-1)/2 + add_letters + 0.7 + concordance_indent) * unit_width, ef_text_y * unit_height, 'ASB Effect Size, log₂ : {:.2f}'.format(ef), size=str(p_value_text_h*unit_height) + 'px', anchor='middle', color='#000000de', font='PT Sans, Arial')
            txt_fdr = transform.TextElement(((m-1)/2 + add_letters + 0.7 + concordance_indent) * unit_width, fdr_text_y * unit_height, 'ASB FDR: ' + get_scientific_text(fdr), size=str(p_value_text_h*unit_height) + 'px', anchor='middle', color='#000000de', font='PT Sans, Arial')
            fig.append([txt_ef, txt_fdr])

            bracket_thick = 10
            conc_label_space = 0.1

            place_concordance_on_svg(fig, tf_snp.motif_concordance, conc_label_space*unit_width, (1/2 +(1/2 + full_gap + text_h + indent + hill_gap + strands_h/2 + ref_height)/2)*unit_height - concordance_h[tf_snp.motif_concordance]/2, concordance_h[tf_snp.motif_concordance], concordance_w)
            place_letter_on_svg(fig, os.path.expanduser('~/letters/rect2.svg'), conc_label_space*2*unit_width+concordance_w, (1/2 + bracket_thick/600/2) * unit_height, (1/2 + full_gap + text_h + indent + hill_gap + strands_h/2 + ref_height - bracket_thick/600)*unit_height, bracket_thick/300*unit_width)
            place_letter_on_svg(fig, os.path.expanduser('~/letters/rect2.svg'), conc_label_space*2*unit_width+concordance_w, (1/2 - bracket_thick/600/2) * unit_height, bracket_thick/600*unit_height, (add_letters + concordance_indent - 3*conc_label_space)*unit_width - concordance_w)
            place_letter_on_svg(fig, os.path.expanduser('~/letters/rect2.svg'), conc_label_space*2*unit_width+concordance_w, (1 + full_gap + text_h + indent + hill_gap + strands_h/2 + ref_height - bracket_thick/600/2) * unit_height, bracket_thick/600*unit_height, (add_letters + concordance_indent - 3*conc_label_space)*unit_width - concordance_w)

            fig.save('D:\Sashok\svgs_{}/{}_{}_{}{}.svg'.format(current_release.name, tf.name, snp.rs_id, snp.alt, '_revcomp' if draw_revcomp else ''))
