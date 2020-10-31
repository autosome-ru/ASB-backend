import os
from ASB_app.models import Ticket
from ASB_app.releases import current_release
import pandas as pd

session = current_release.session


def get_ticket_id_from_path(file_path):
    return os.path.basename(file_path).replace('.tsv', '')


def get_ticket(ticket_id):
    return Ticket.query.get_or_404(ticket_id)


def create_processed_path(ticket_id, *subdirs):
    ticket_dir = get_path_by_ticket_id(ticket_id, path_type='dir', ext='')
    if not os.path.isdir(ticket_dir):
        os.mkdir(ticket_dir)
    for subdir in subdirs:
        path = os.path.join(ticket_dir, subdir)
        if not os.path.isdir(path):
            os.mkdir(path)


def get_tickets_dir(suffix=''):
    return os.path.join(os.path.expanduser('~/'), 'adastra', 'tickets', suffix)


def get_path_by_ticket_id(ticket_id, path_type='input', ext='.tsv'):
    if path_type == 'dir':
        ext = ''
    return os.path.join(
        get_tickets_dir(),
        *{
            'input': ['accepted'],
            'dir': ['processed'],
            'tf': ['processed', '{}'.format(ticket_id), 'tf'],
            'cl': ['processed', '{}'.format(ticket_id), 'cl'],
            'tf_sum': ['processed', '{}'.format(ticket_id), 'tf_sum'],
            'cl_sum': ['processed', '{}'.format(ticket_id), 'cl_sum'],
        }[path_type],
        '{}{}'.format(ticket_id, ext)
    )


def create_ticket(ticket_id):
    ticket = Ticket(
        ticket_id=ticket_id,
        status='Created',
    )
    session.add(ticket)
    session.commit()
    return ticket


def update_ticket_status(ticket_id, status):
    ticket = get_ticket(ticket_id)
    ticket.status = status
    session.commit()


def delete_ticket(ticket_id):
    ticket = get_ticket(ticket_id)
    if ticket.status == 'Processing':
        return False
    input_file = get_path_by_ticket_id(ticket_id)
    if os.path.isfile(input_file):
        os.remove(input_file)
    ticket_dir = get_path_by_ticket_id(ticket_id, path_type='dir')
    if os.path.isdir(ticket_dir):
        for path_type in ('tf', 'cl', 'tf_sum', 'cl_sum'):
            report = get_path_by_ticket_id(ticket_id, path_type=path_type)
            if os.path.isfile(report):
                os.remove(report)
        for dir in os.listdir(ticket_dir):
            os.rmdir(os.path.join(ticket_dir, dir))
        os.rmdir(ticket_dir)
    session.delete(ticket)
    session.commit()
    return True


def get_result(ticket_id, param):
    ticket = get_ticket(ticket_id)
    if ticket.status != 'Processed':
        return False, {}
    result = []
    out_file = get_path_by_ticket_id(ticket_id, path_type=param)
    with open(out_file) as out:
        header = []
        for number, line in enumerate(out):
            if number == 1001:
                break
            if number == 0:
                header = [x.lower() for x in line.strip('\n').split('\t')]
                continue
            result.append(dict(zip(header, line.strip('\n').split('\t'))))
    return True, result


def delete_all_tickets():
    for ticket in Ticket.query:
        ticket_id = ticket.ticket_id
        if ticket.status == 'Processing':
            continue
        input_file = get_path_by_ticket_id(ticket_id)
        if os.path.isfile(input_file):
            os.remove(input_file)
        ticket_dir = get_path_by_ticket_id(ticket_id, path_type='dir')
        if os.path.isdir(ticket_dir):
            for path_type in ('tf', 'cl', 'tf_sum', 'cl_sum'):
                report = get_path_by_ticket_id(ticket_id, path_type=path_type)
                if os.path.isfile(report):
                    os.remove(report)
            for dir in os.listdir(ticket_dir):
                os.rmdir(os.path.join(ticket_dir, dir))
            os.rmdir(ticket_dir)
        session.delete(ticket)
        session.commit()


def annotate_with_context(filename, field):
    SNP = current_release.SNP
    t = pd.read_table(os.path.expanduser(filename))
    rs_ids = [int(x[2:]) for x in t[field].tolist()]
    ann = ['rs' + str(snp.rs_id) + ' ' + snp.context[:24] + '[' + snp.context[24] + '/' + snp.alt + ']' + snp.context[25:] for snp in SNP.query.filter(SNP.rs_id.in_(rs_ids))]
    return '\n'.join(ann)
