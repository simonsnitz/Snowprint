from api import db

class Cluster(db.Model):
    __searchable__ = ['body']
    __tablename__ = "cluster"

    id = db.Column(db.Integer, primary_key=True)

    cluster_percent_identity = db.Column(db.Integer, index=True)
    family = db.Column(db.String(16), index=True)



    def __repr__(self):
        return '<Cluster {}>'.format(self.id)