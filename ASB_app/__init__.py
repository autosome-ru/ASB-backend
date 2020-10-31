from flask import Flask, Blueprint
from flask_migrate import Migrate
from flask_restplus import Api
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import MetaData
from flask_apscheduler import APScheduler
from flask_executor import Executor
from flask_cors import CORS
from ASB_app.releases import Release

from config import Config

app = Flask(__name__)
app.config.from_object(Config)
CORS(app)
logger = app.logger
logger.setLevel(app.config['LOGGER_LEVEL'])

naming_convention = {
    "ix": 'ix_%(column_0_label)s',
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(column_0_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s"
}


for release in Release.__subclasses__():
    blueprint = Blueprint('api_{}'.format(release.name), __name__, url_prefix='/api/v{}'.format(release.version))
    new_api = Api(blueprint, version=release.full_version, title='ADASTRA API', description='ADASTRA API (release {})'.format(release.name))
    setattr(release, 'api', new_api)
    app.register_blueprint(blueprint)
    new_db = SQLAlchemy(app, metadata=MetaData(naming_convention=naming_convention))
    setattr(release, 'db', new_db)
    setattr(release, 'session', new_db.session)
    setattr(release, 'migrate', Migrate(app, new_db))

scheduler = APScheduler(app=app)
executor = Executor(app)
scheduler.start()


# @app.after_request
# def apply_caching(response):
#     response.headers['Access-Control-Allow-Origin'] = '*'
#     return response


from . import models
from . import routes
from . import utils
