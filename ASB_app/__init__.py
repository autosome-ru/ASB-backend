from flask import Flask, Blueprint
from flask_migrate import Migrate
from flask_restplus import Api
from flask_sqlalchemy import SQLAlchemy
from jsonschema import FormatChecker
from sqlalchemy import MetaData

from config import Config

# from flask_sslify import SSLify


app = Flask(__name__)
app.config.from_object(Config)
logger = app.logger
logger.setLevel(app.config['LOGGER_LEVEL'])

naming_convention = {
    "ix": 'ix_%(column_0_label)s',
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(column_0_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s"
}
db = SQLAlchemy(app, metadata=MetaData(naming_convention=naming_convention))
session = db.session


blueprint = Blueprint('api', __name__, url_prefix='/api/v1')
api = Api(blueprint, version='1.0', title='AD ASTRA API', description='AD ASTRA API', format_checker=FormatChecker())
app.register_blueprint(blueprint)


migrate = Migrate(app, db)


@app.after_request
def apply_caching(response):
    response.headers['Access-Control-Allow-Origin'] = '*'
    return response


import ASB_app.models as models
from . import routes
