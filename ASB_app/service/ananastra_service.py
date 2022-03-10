import csv
import os
import shutil
import tempfile
from datetime import datetime, timedelta
import pandas as pd
from hashlib import md5

from flask import send_file

from ASB_app import executor, logger
from ASB_app.exceptions import FileNotProcessed
from ASB_app.models import Ticket
from ASB_app.releases import current_release
from ASB_app.utils.aggregates import TsvDialect

session = current_release.session


def get_ticket_id_from_path(file_path):
    return os.path.basename(file_path).replace('.tsv', '')


def get_ticket(ticket_id):
    return Ticket.query.get_or_404(ticket_id)


def create_processed_path(ticket_id):
    ticket_dir = get_path_by_ticket_id(ticket_id, path_type='dir', ext='')
    if not os.path.isdir(ticket_dir):
        os.mkdir(ticket_dir)


def get_tickets_dir(suffix='',):
    return os.path.join(os.path.expanduser('~/'), 'adastra',
                        'tickets', suffix)


def get_path_by_ticket_id(ticket_id, path_type='input', ext='.tsv'):
    if path_type == 'dir':
        return get_tickets_dir(ticket_id)
    else:
        return os.path.join(
            get_tickets_dir('accepted' if path_type == 'input' else get_tickets_dir(ticket_id)),
            '{}{}{}'.format(ticket_id, '' if path_type == 'input' else '.' + path_type, ext)
        )


def create_ticket(ticket_id, user_id):
    ticket = Ticket(
        ticket_id=ticket_id,
        status='Created',
        user_id=user_id,
        date_created=datetime.now(),
        fdr=None,
        expiration_date=datetime.now() + timedelta(days=3)
    )
    session.add(ticket)
    session.commit()
    return ticket


def update_ticket_status(ticket_id, status):
    ticket = get_ticket(ticket_id)
    ticket.status = status
    session.commit()


def log_hash(output):
    with open(os.path.join(get_tickets_dir('logs'), 'hash_sums.log'), 'a') as out:
        out.write(output + '\n')


def start_processing_ticket(ticket_id):
    calc_hash = md5()
    update_ticket_status(ticket_id, 'Processing')
    with open(get_path_by_ticket_id(ticket_id), 'rb') as r:
        calc_hash.update(r.read())
    output = calc_hash.hexdigest()
    log_hash(output)


def delete_ticket(ticket_id):
    ticket = get_ticket(ticket_id)
    if ticket.status == 'Processing':
        logger.info('Deleting a ticket in process: {}'.format(ticket_id))
    ticket_dir = get_path_by_ticket_id(ticket_id, path_type='dir')
    if os.path.isdir(ticket_dir):
        shutil.rmtree(ticket_dir)
    session.delete(ticket)
    session.commit()
    return True


def modify_null(field):
    if field == 'None':
        return None
    if pd.isna(field):
        return None
    return field


def get_result(ticket_id, param, limit, format):
    ticket = get_ticket(ticket_id)
    if ticket.status != 'Processed':
        raise FileNotProcessed
    out_file = get_path_by_ticket_id(ticket_id, path_type=param)
    if format == 'json':
        out = pd.read_table(out_file, na_values=['None', 'NaN', 'nan', '', 'NULL'],
                            nrows=None if limit == 0 else limit)
        new_header = {x: x.lower() for x in out.columns}
        out.rename(columns=new_header, inplace=True)
        return out.to_json(orient='records')
    elif format == 'tsv':
        file = tempfile.NamedTemporaryFile('wt', suffix='.tsv')
        csv_writer = csv.writer(file, dialect=TsvDialect)
        with open(out_file) as out:
            for number, line in enumerate(out):
                if limit != 0 and number == limit + 1:
                    break
                if param != 'not_found':
                    if number == 0:
                        csv_writer.writerow([x.lower() for x in line.strip('\n').split('\t')])
                        continue
                csv_writer.writerow(line.strip('\n').split('\t'))
        file.flush()
        return send_file(
            file.name,
            cache_timeout=0,
            mimetype="text/tsv",
            as_attachment=True
        )


def delete_all_tickets():
    for ticket in Ticket.query:
        ticket_id = ticket.ticket_id
        if ticket.status == 'Processing':
            continue
        ticket_dir = get_path_by_ticket_id(ticket_id, path_type='dir')
        if os.path.isdir(ticket_dir):
            shutil.rmtree(ticket_dir)
        session.delete(ticket)
        session.commit()


def get_ticket_position_in_queue(ticket_id):
    count = 0
    for id, future in executor.futures._futures.items():
        if id == ticket_id:
            return count
        if getattr(future, '__proxied_object')._state in ('RUNNING', 'PENDING'):
            count += 1
