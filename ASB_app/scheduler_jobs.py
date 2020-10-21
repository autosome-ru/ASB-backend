from ASB_app import scheduler, session, logger
from ASB_app.models import Ticket
from ASB_app.service import ananastra_service
from datetime import datetime, timedelta


@scheduler.task('cron', id='TicketScheduller', minute='0,15,30,45', max_instances=1)
def check_outdated_tickets():
    logger.info('Started tickets check')
    for ticket_id in session.query(Ticket.ticket_id).filter(datetime.now() - Ticket.date_created > timedelta(hours=1)):
        ok = ananastra_service.delete_ticket(ticket_id)
        if not ok:
            logger.info('Could not delete ticket: {}'.format(ticket_id))
