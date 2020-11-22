from ASB_app import scheduler, logger
from ASB_app.models import Ticket
from ASB_app.service import ananastra_service
from datetime import datetime, timedelta

from ASB_app.releases import current_release
session = current_release.session


@scheduler.task('cron', id='TicketScheduller', hour='*', max_instances=1)
def check_outdated_tickets():
    logger.info('Started tickets check')
    for (ticket_id, ) in session.query(Ticket.ticket_id).filter(datetime.now() > Ticket.expiration_date):
        ok = ananastra_service.delete_ticket(ticket_id)
        if not ok:
            logger.info('Could not delete ticket: {}'.format(ticket_id))
    for (ticket_id, ) in session.query(Ticket.ticket_id).filter(Ticket.status == 'Failed'):
        ok = ananastra_service.delete_ticket(ticket_id)
        if not ok:
            logger.info('Could not delete ticket: {}'.format(ticket_id))
    logger.info('Finished tickets check')
