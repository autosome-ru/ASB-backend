from ASB_app.releases import current_release
import sys
import os
import pandas as pd

SNP = current_release.SNP
t = pd.read_table(os.path.expanduser(sys.argv[1]), header=None, comment='#')
rs_ids = [int(x[2:]) for x in t[field].tolist()]
out = pd.DataFrame()

ann = []
snps = SNP.query.filter(SNP.rs_id.in_(rs_ids)).all()
for snp in snps:
    ann.append('rs' + str(snp.rs_id) + ' ' + snp.context[:24] + '[' + snp.context[24] + '/' + snp.alt + ']' + snp.context[25:])

out['context'] = ann

out.to_csv(os.path.expanduser(sys.argv[2]), sep='\t', index=False, header=False)
