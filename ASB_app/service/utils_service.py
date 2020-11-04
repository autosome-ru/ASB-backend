import os
import tempfile

from flask import send_file
import pandas as pd


def annotate_with_context(release, filename, field):
    SNP = release.SNP
    t = pd.read_table(os.path.expanduser(filename))
    rs_ids = [int(x[2:]) for x in t[field].tolist()]
    ann = ['rs' + str(snp.rs_id) + ' ' + snp.context[:24] + '[' + snp.context[24] + '/' + snp.alt + ']' + snp.context[25:] for snp in SNP.query.filter(SNP.rs_id.in_(rs_ids))]
    file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
    file.write('\n'.join(ann))
    file.flush()

    return send_file(
        file.name,
        cache_timeout=0,
        mimetype="text/tsv",
        as_attachment=True
    )