import os
import pathlib
dir = pathlib.Path(__file__).parent.absolute()


class Config(object):
    # SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(dir, 'AD_ASTRA.db')
    # 'mysql://local_user:local_password9@localhost:3306/lct?charset=utf8'
    SQLALCHEMY_BINDS = {
        'soos': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_soos?charset=utf8',
        'ford': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_ford?charset=utf8',
        'dan': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_dan?charset=utf8',
        # 'dnase': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_dnase?charset=utf8',
        'candidates': 'mysql://adastra:pass=ADASTRA880@localhost:3306/adastra_candidates?charset=utf8',
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
