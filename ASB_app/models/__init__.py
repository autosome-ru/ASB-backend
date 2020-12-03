from .placeholders import *
from .adastra_models import *
from .ananastra_models import *
from ASB_app.releases import Release, current_release


for release in Release.__subclasses__():
    release.db.create_all(release.name)
    release.session.commit()

current_release.db.create_all('candidates')
current_release.db.create_all('tickets')
