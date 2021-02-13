from ASB_app.exceptions import ReleaseNotFound


class Release:
    db = None
    session = None
    api = None
    name = 'abstract'
    version = '0'
    full_version = '0.0'


class ReleaseSoos(Release):
    name = 'soos'
    version = '1'
    full_version = '1.6'


class ReleaseFord(Release):
    name = 'ford'
    version = '2'
    full_version = '2.1'


class ReleaseDan(Release):
    name = 'dan'
    version = '3'
    full_version = '3.0'


def get_release_by_version(version):
    for release in Release.__subclasses__():
        if release.version == version:
            return release
    raise ReleaseNotFound('No release: v{}'.format(version))


current_release = ReleaseDan
