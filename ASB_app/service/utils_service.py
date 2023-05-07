import os
import tempfile

from flask import send_file
import pandas as pd

from ASB_app.models import abstract_models


class UtilsService:
    def __init__(self, release):
        self.release = release
        for model in abstract_models:
            setattr(self, model.__name__, getattr(release, model.__name__))

    def annotate_with_context(self, filename, field):
        SNP = self.release.SNP
        try:
            t = pd.read_table(os.path.expanduser(filename))
            rs_ids = [int(x[2:]) for x in t[field].tolist()]
        except:
            return False, {}
        ann = ['rs' + str(snp.rs_id) + ' ' + snp.context[:24] + '[' + snp.context[24] + '/' + snp.alt + ']' + snp.context[25:] for snp in SNP.query.filter(SNP.rs_id.in_(rs_ids))]
        file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
        file.write('\n'.join(ann))
        file.flush()

        return True, send_file(
            file.name,
            cache_timeout=0,
            mimetype="text/tsv",
            as_attachment=True
        )
