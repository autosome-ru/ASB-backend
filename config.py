import os


class Config(object):
    SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.expanduser('~/AD_ASTRA_database/AD_ASTRA.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    RESTPLUS_VALIDATE = True
    LOGGER_LEVEL = 'DEBUG'
