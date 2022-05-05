import os
import pathlib
dir = pathlib.Path(__file__).parent.absolute()


class Config(object):
    SQLALCHEMY_BINDS = {
        'soos': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_soos?charset=utf8',
        'susan': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_susan?charset=utf8',
        'zanthar': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_zanthar?charset=utf8',
        'candidates_zanthar': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_candidates_zanthar?charset=utf8',
        'billcipher': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_billcipher?charset=utf8',
        'tickets': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_tickets?charset=utf8',
    }
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    RESTPLUS_VALIDATE = True
    LOGGER_LEVEL = 'DEBUG'
    MAX_CONTENT_LENGTH = 2 * 1024 * 1024
    EXECUTOR_MAX_WORKERS = 1
    RAW_MYSQL_USERNAME = 'adastra'
    RAW_MYSQL_PASSWORD = 'pass=ADASTRA880'
    RAW_MYSQL_HOSTNAME = 'localhost'
