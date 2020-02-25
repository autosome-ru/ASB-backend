export FLASK_APP=./run.py
DB=~/AD_ASTRA_database/AD_ASTRA.db

cp $DB "$HOME/ad_astra_backup_$(date +%d%m).db"

. venv/bin/activate

if [ -f ./pre-migration.sql ]; then
    sqlite3 $DB <pre-migration.sql;
fi

flask db migrate
flask db upgrade

if [ -f ./post-migration.sql ]; then
    sqlite3 $DB <post-migration.sql;
fi
