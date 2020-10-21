from ASB_app import db
from datetime import datetime


class Ticket(db.Model):
    __tablename__ = 'tickets'

    ticket_id = db.Column(db.String(50), primary_key=True)
    status = db.Column(db.Enum('Created', 'Processing', 'Processed', 'Failed'), server_default='Created', nullable=False)
    date_created = db.Column(db.DATETIME, default=datetime.now(), nullable=False)
    meta_info = db.Column(db.JSON, server_default='{}')

    def __repr__(self):
        return '<AnanastraTicket {0.ticket_id}, created {0.date_created}, {0.status}>'.format(self)
