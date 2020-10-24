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
    full_version = '1.2'


class ReleaseTest(Release):
    name = 'test'
    version = '2'
    full_version = '2.0'


current_release = ReleaseSoos
