from ASB_app import *
from ASB_app.releases import current_release

S = current_release.SNP
session = current_release.session

with open('D:\\Sashok\\snps_{}.bed'.format(current_release.name), 'w', encoding='utf-8') as f:
    for s in S.query:
        assert s.fdr_class and s.fdr_class != '1'
        f.write('\t'.join(map(str, [s.chromosome, s.position - 1, s.position, 'rs' + str(s.rs_id)])) + '\n')
