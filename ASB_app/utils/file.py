def pack(line):
    return '\t'.join(map(str, line)) + '\n'


def process_row(row, param, header):
    row_dict = dict(zip(header, row))
    row_dict['PEAK_CALLS'] = 'None' if row_dict['PEAK_CALLS'] == 0 else 'Single' if row_dict[
                                                                                        'PEAK_CALLS'] == 1 else 'Multiple'
    row_dict['RS_ID'] = 'rs' + str(row_dict['RS_ID'])
    if param == 'TF':
        row_dict['MOTIF_ORIENTATION'] = '+' if row_dict['MOTIF_ORIENTATION'] else '-' if row_dict[
                                                                                             'MOTIF_ORIENTATION'] == 0 else 'None'
        row_dict['MOTIF_CONCORDANCE'] = {
            'Concordant': 'Concordant',
            'Discordant': 'Discordant',
            'Weak Concordant': 'Weak concordant',
            'Weak Discordant': 'Weak discordant',
            'No Hit': 'No hit',
            None: 'None',
        }[row_dict['MOTIF_CONCORDANCE']]
    return tuple(row_dict[name] for name in header)
