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


# class ReleaseSusan(Release):
#     name = 'susan'
#     version = '3'
#     full_version = '3.5'
#
#
# class ReleaseZanthar(Release):
#     name = 'zanthar'
#     version = '4'
#     full_version = '4.0'


class ReleaseBillCipher(Release):
    name = 'billcipher'
    version = '5'
    full_version = '5.1'


def get_release_by_version(version):
    for release in Release.__subclasses__():
        if release.version == version:
            return release
    raise ReleaseNotFound('No release: v{}'.format(version))


current_release = ReleaseBillCipher
