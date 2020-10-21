import os
from ASB_app.models import Ticket
from ASB_app import session


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
        tf_report = get_path_by_ticket_id(ticket_id, path_type='tf')
        cl_report = get_path_by_ticket_id(ticket_id, path_type='cl')
        if os.path.isfile(tf_report):
            os.remove(tf_report)
        if os.path.isfile(cl_report):
            os.remove(cl_report)
        for dir in os.listdir(ticket_dir):
            os.rmdir(dir)
        os.remove(ticket_dir)
    session.delete(ticket)
    session.commit()
    return True


# def get_result(ticket_id, param):
#     ticket = get_ticket(ticket_id)
#     if ticket.status != ''