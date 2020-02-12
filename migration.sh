export FLASK_APP=./run.py
DB=AD_ASTRA.db

mysqldump $DB >"$HOME/lct_backup_$(date +%d%m).sql"

. venv/bin/activate

if [ -f ./pre-migration.sql ]; then
    mysql $DB <pre-migration.sql;
fi

flask db migrate
flask db upgrade

if [ -f ./post-migration.sql ]; then
    mysql $DB <post-migration.sql;
fi
