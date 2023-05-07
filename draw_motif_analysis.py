from svgutils import transform
import numpy as np
from math import lgamma
import os
from tqdm import tqdm



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

    letters = ['A', 'C', 'G', 'T']

    unit_width = 30
    unit_height = 60

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
    concordance_scale = 10/(3351.318*2) * unit_width
    conc_self_width = 331.204

    visible_cut_tr = 0.05

    indent = snp_text_h + snp_gap / 2

    full_gap = motif_gap + indent


    letter_svgs = {
        'A': os.path.expanduser('~/PARAMETERS/letters/lettarA_path.svg'),
        'C': os.path.expanduser('~/PARAMETERS/letters/lettarC_path.svg'),
        'G': os.path.expanduser('~/PARAMETERS/letters/lettarG_path.svg'),
        'T': os.path.expanduser('~/PARAMETERS/letters/lettarT_path.svg'),
    }

    black_letter_svgs = {
        'A': os.path.expanduser('~/letters/lettarA_path_black.svg'),
        'C': os.path.expanduser('~/letters/lettarC_path_black.svg'),
        'G': os.path.expanduser('~/letters/lettarG_path_black.svg'),
        'T': os.path.expanduser('~/letters/lettarT_path_black.svg'),
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

    trimmed_snp_context = "CATCAGGGAAGGCTCCATAGGCTG"
    # print(full_gap, text_h, indent)
    fig = transform.SVGFigure(800, 200)

    pos_in_motif = 13
    motif_context = trimmed_snp_context
    # pos_in_motif += add_letters


    # print(motif_context, pos_in_motif)
    #
    # place_letter_on_svg(fig, os.path.expanduser('~/PARAMETERS/letters/rect.svg'), (pos_in_motif + concordance_indent) * unit_width, 0, (1 + full_gap + text_h/2 + snp_gap/2 + snp_text_h) * unit_height, unit_width)


    for pos, letter in enumerate(motif_context):
        # if letter == ' ':
        #     continue
        # if pos != pos_in_motif:
        place_letter_on_svg(fig, black_letter_svgs[letter.upper()], (pos + (1 - text_width) / 2 + concordance_indent) * unit_width, (1 + full_gap) * unit_height, text_h * unit_height, text_width * unit_width)
        # else:
        #     place_letter_on_svg(fig, letter_svgs[ref], (pos + (1 - snp_text_width) / 2 + concordance_indent) * unit_width, (1 + full_gap + text_h/2 - snp_text_h - snp_gap/2) * unit_height, snp_text_h * unit_height, snp_text_width * unit_width)
        #     place_letter_on_svg(fig, letter_svgs[alt], (pos + (1 - snp_text_width) / 2 + concordance_indent) * unit_width, (1 + full_gap + text_h/2 + snp_gap/2) * unit_height, snp_text_h * unit_height, snp_text_width * unit_width)


    # text_x = pos_in_motif + 1 + concordance_indent
    # letter_text = pos_in_motif + 1 - 0.05
    # txt_ref = transform.TextElement(text_x * unit_width, (1.24 + (motif_gap + indent)/2 + p_value_text_h/2) * unit_height, 'P-value: ' + get_scientific_text(motif_pref), size=str(p_value_text_h*unit_height), color='#000000de', font='PT Sans, Arial')
    # motif_word_ref = transform.TextElement(text_x * unit_width, (1.24 + (motif_gap + indent)/2 - p_value_text_h/2) * unit_height, 'Motif', size=str(p_value_text_h*unit_height), color='#000000de', font='PT Sans, Arial')
    # txt_alt = transform.TextElement(text_x * unit_width, (1 + motif_gap + indent*3/2 + text_h + p_value_text_h/2) * unit_height, 'P-value: ' + get_scientific_text(motif_palt), size=str(p_value_text_h*unit_height), color='#000000de', font='PT Sans, Arial')
    # motif_word_alt = transform.TextElement(text_x * unit_width, (1 + motif_gap + indent*3/2 + text_h - p_value_text_h/2) * unit_height, 'Motif', size=str(p_value_text_h*unit_height), color='#000000de', font='PT Sans, Arial')
    # txt_ref_letter = transform.TextElement((letter_text) * unit_width,
    #                                 (1.24 + (motif_gap + indent)/2 - p_value_text_h/2) * unit_height,
    #                                 'Ref',
    #                                 size=str(p_value_text_h * unit_height), color='#000000de',
    #                                 font='PT Sans, Arial', anchor='end')
    # txt_alt_letter = transform.TextElement((letter_text) * unit_width, (1 + motif_gap + indent*3/2 + text_h + p_value_text_h/2) * unit_height, 'Alt',
    #                                        size=str(p_value_text_h * unit_height), color='#000000de',
    #                                        font='PT Sans, Arial', anchor='end')
    #
    # fig.append([txt_ref, txt_alt, motif_word_ref, motif_word_alt, txt_ref_letter, txt_alt_letter])
    #
    # if asb_is_ref:
    #     ef_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height*1 - p_value_text_h/2
    #     fdr_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height*(1-1/3) - p_value_text_h/2
    # else:
    #     ef_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height + strands_h + alt_height*1/3 - p_value_text_h/2
    #     fdr_text_y = 1 + full_gap + text_h + indent + hill_gap + ref_height + strands_h + alt_height*2/3 - p_value_text_h/2
    #
    # txt_ef = transform.TextElement(((m-1)/2 + add_letters + 0.7 + concordance_indent) * unit_width, ef_text_y * unit_height, 'ASB Effect Size, log₂ : {:.2f}'.format(ef), size=str(p_value_text_h*unit_height), anchor='middle', color='#000000de', font='PT Sans, Arial')
    # txt_fdr = transform.TextElement(((m-1)/2 + add_letters + 0.7 + concordance_indent) * unit_width, fdr_text_y * unit_height, 'ASB FDR: ' + get_scientific_text(fdr), size=str(p_value_text_h*unit_height), anchor='middle', color='#000000de', font='PT Sans, Arial')
    # fig.append([txt_ef, txt_fdr])
    #
    # bracket_thick = 10
    # conc_label_space = 0.1
    #
    # place_letter_on_svg(fig, os.path.expanduser('~/PARAMETERS/letters/rect2.svg'), conc_label_space*2*unit_width+concordance_w, (1/2 + bracket_thick/600/2) * unit_height, (1/2 + full_gap + text_h + indent + hill_gap + strands_h/2 + ref_height - bracket_thick/600)*unit_height, bracket_thick/300*unit_width)
    # place_letter_on_svg(fig, os.path.expanduser('~/PARAMETERS/letters/rect2.svg'), conc_label_space*2*unit_width+concordance_w, (1/2 - bracket_thick/600/2) * unit_height, bracket_thick/600*unit_height, (add_letters + concordance_indent - 3*conc_label_space)*unit_width - concordance_w)
    # place_letter_on_svg(fig, os.path.expanduser('~/PARAMETERS/letters/rect2.svg'), conc_label_space*2*unit_width+concordance_w, (1 + full_gap + text_h + indent + hill_gap + strands_h/2 + ref_height - bracket_thick/600/2) * unit_height, bracket_thick/600*unit_height, (add_letters + concordance_indent - 3*conc_label_space)*unit_width - concordance_w)

    fig.save("output.svg")
