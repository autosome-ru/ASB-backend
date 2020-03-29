import os
import pathlib
dir = pathlib.Path(__file__).parent.absolute()


class Config(object):
    SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(dir, 'AD_ASTRA.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    RESTPLUS_VALIDATE = True
    LOGGER_LEVEL = 'DEBUG'
