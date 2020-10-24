from .placeholders import *
from .adastra_models import *
from .ananastra_models import *
from ASB_app.releases import Release, current_release


for release in Release.__subclasses__():
    release.db.create_all(release.name)
    release.session.commit()

current_release.db.create_all('candidates')

SNP = current_release.SNP
session = current_release.session

possible_tf_asbs = SNP.query.filter(
    SNP.tf_aggregated_snps.any()
).group_by(SNP.rs_id).count()

possible_cl_asbs = SNP.query.filter(
    SNP.tf_aggregated_snps.any()
).group_by(SNP.rs_id).count()

possible_all_asbs = session.query(SNP.rs_id.distinct()).count()

possible_tf_candidates = CandidateSNP.query.filter_by(ag_level='TF').group_by(CandidateSNP.rs_id).count()

possible_cl_candidates = CandidateSNP.query.filter_by(ag_level='CL').group_by(CandidateSNP.rs_id).count()

possible_all_candidates = session.query(CandidateSNP.rs_id.distinct()).count()
