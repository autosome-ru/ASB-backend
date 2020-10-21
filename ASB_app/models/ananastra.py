from ASB_app import db


class Ticket(db.Model):
    __tablename__ = 'tickets'

    ticket_id = db.Column(db.String(50), primary_key=True)
    status = db.Column(db.Enum('Created', 'Processing', 'Processed', 'Failed'))
    meta_info = db.Column(db.JSON)

    def __repr__(self):
        return '<AnanastraTicket {0.ticket_id}, {0.status}>'.format(self)
