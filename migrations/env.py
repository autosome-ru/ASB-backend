from __future__ import with_statement
from alembic import context
from sqlalchemy import engine_from_config, pool
from logging.config import fileConfig
import logging

# this is the Alembic Config object, which provides
# access to the values within the .ini file in use.
config = context.config

# Interpret the config file for Python logging.
# This line sets up loggers basically.
fileConfig(config.config_file_name)
logger = logging.getLogger('alembic.env')

# add your model's MetaData object here
# for 'autogenerate' support
# from myapp import mymodel
# target_metadata = mymodel.Base.metadata
from flask import current_app
config.set_main_option('sqlalchemy.url',
                       current_app.config.get('SQLALCHEMY_DATABASE_URI'))
target_metadata = current_app.extensions['migrate'].db.metadata

# other values from the config, defined by the needs of env.py,
# can be acquired:
# my_important_option = config.get_main_option("my_important_option")
# ... etc.


def my_compare_type(context, inspected_column, metadata_column, inspected_type, metadata_type):
    """
    To automatically detect changes in string length (String(50) -> String(60))
     one can provide compare_type parameter to context (disabled by default).
    Default algorithm tries to convert mysql type TINYINT(1) to Boolean.
    Since mysql has no 'fair' Boolean type and uses TINYINT(1) instead, it happens each time.
    This function overrides default comparator to avoid this problem.

    Links:
    https://github.com/miguelgrinberg/Flask-Migrate/issues/143
    https://alembic.sqlalchemy.org/en/latest/autogenerate.html#comparing-types
    """
    # This patch is required only for mysql database, otherwise use default algorithm
    if context.dialect.name != 'mysql':
        return None

    # tinyint(1) -> bool
    if metadata_type.python_type is bool and inspected_type.python_type is int and getattr(inspected_type, 'display_width') == 1:
        return False
    else:
        # default algorithm
        return None


def run_migrations_offline():
    """Run migrations in 'offline' mode.

    This configures the context with just a URL
    and not an Engine, though an Engine is acceptable
    here as well.  By skipping the Engine creation
    we don't even need a DBAPI to be available.

    Calls to context.execute() here emit the given string to the
    script output.

    """
    url = config.get_main_option("sqlalchemy.url")
    context.configure(url=url,
                      include_object=include_object,
                      compare_type=my_compare_type,
                      compare_server_default=True)

    with context.begin_transaction():
        context.run_migrations()


def include_object(object, name, type_, reflected, compare_to):
    if type_ == 'table' and name in ('GAMAUNSUBS', 'VKADS', 'VKGROUPS', 'audio_messages', 'audio_sends', 'lectparameters', 'processed_rooms', 'reviews', 'webinars_stats', 'webinars_temp', 'vk_teachers_audio'):
        return False
    return True


def run_migrations_online():
    """Run migrations in 'online' mode.

    In this scenario we need to create an Engine
    and associate a connection with the context.

    """

    # this callback is used to prevent an auto-migration from being generated
    # when there are no changes to the schema
    # reference: http://alembic.zzzcomputing.com/en/latest/cookbook.html
    def process_revision_directives(context, revision, directives):
        if getattr(config.cmd_opts, 'autogenerate', False):
            script = directives[0]
            if script.upgrade_ops.is_empty():
                directives[:] = []
                logger.info('No changes in schema detected.')

    engine = engine_from_config(config.get_section(config.config_ini_section),
                                prefix='sqlalchemy.',
                                poolclass=pool.NullPool)

    connection = engine.connect()
    context.configure(connection=connection,
                      target_metadata=target_metadata,
                      process_revision_directives=process_revision_directives,
                      include_object=include_object,
                      compare_type=my_compare_type,
                      # compare_server_default=True,  # Raises error since migration passes None to re.sub
                      render_as_batch=True,
                      **current_app.extensions['migrate'].configure_args)

    try:
        with context.begin_transaction():
            context.run_migrations()
    finally:
        connection.close()

if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
