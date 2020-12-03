import os
import pathlib
dir = pathlib.Path(__file__).parent.absolute()


class Config(object):
    # SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(dir, 'AD_ASTRA.db')
    SQLALCHEMY_BINDS = {
        'soos': 'sqlite:///' + os.path.join(dir, 'AD_ASTRA_soos.db'),
        'ford': 'sqlite:///' + os.path.join(dir, 'AD_ASTRA_ford.db'),
        'candidates': 'sqlite:///' + os.path.join(dir, 'AD_ASTRA_candidates.db'),
        'tickets': 'sqlite:///' + os.path.join(dir, 'AD_ASTRA_tickets.db'),
    }
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    RESTPLUS_VALIDATE = True
    LOGGER_LEVEL = 'DEBUG'
    MAX_CONTENT_LENGTH = 2 * 1024 * 1024
    EXECUTOR_MAX_WORKERS = 1
