from api import db


class Alignment(db.Model):
    __searchable__ = ['body']
    __tablename__ = "alignment"

    id = db.Column(db.Integer, primary_key=True, unique=True)

        # Protein accession ID of the BLASTed regulator.
    query_id = db.Column(db.String(16), index=True, unique=True)

        # Metadata of homologs
        # [ {"acc": accession_ID, "identity": percent_identity, "qcoverage": query_coverage}, {...}, ...]
    homologs = db.Column(db.Text(4096))     # This will be a stringified JSON object

    def __repr__(self):
        return '<Alignment {}>'.format(self.id)


class Regulator(db.Model):
    __searchable__ = ['body']
    __tablename__ = "regulator"

    id = db.Column(db.Integer, primary_key=True, unique=True)

    prot_id = db.Column(db.String(16), index=True, unique=True)
    genome_id = db.Column(db.String(16))
    organism = db.Column(db.String(128))
    start_pos = db.Column(db.Integer)
    stop_pos = db.Column(db.Integer)
    strand = db.Column(db.String(4))
    organism_id = db.Column(db.Integer)

    operons = db.relationship("Association", back_populates="regulator")

    operator = db.relationship("Operator", back_populates="regulator")
    operator_id = db.Column(db.Integer, db.ForeignKey('operator.id'))

    def __repr__(self):
        return '<Regulator {}>'.format(self.id)


class Operon(db.Model):
    __searchable__ = ['body']
    __tablename__ = "operon"

    id = db.Column(db.Integer, primary_key=True, unique=True)

        # Sequence of the entire extracted operon
    operon_seq = db.Column(db.String(4096))

        # Operon data
        # [ {"acc": string, "annotation": str, "direction": str, "start": int, "stop": int}, {...}, ... ]
    operon = db.Column(db.Text(4096))

        # New columns which help *define* an operon
    genome_id = db.Column(db.String(16))
    start_pos = db.Column(db.Integer)
    stop_pos = db.Column(db.Integer)

    regulators = db.relationship("Association", back_populates="operon")


    def __repr__(self):
        return '<Operon {}>'.format(self.id)


class Association(db.Model):
    __searchable__ = ['body']
    __tablename__ = "association"

    id = db.Column(db.Integer, primary_key=True, unique=True, autoincrement=True)

    regulator_id = db.Column(db.ForeignKey(Regulator.id))
    operon_id = db.Column(db.ForeignKey(Operon.id))

    reg_index = db.Column(db.Integer)
    reg_type = db.Column(db.Integer)
    regulated_seq = db.Column(db.String(4096))

    operon = db.relationship("Operon", back_populates="regulators")
    regulator = db.relationship("Regulator", back_populates="operons")


    def __repr__(self):
        return '<Association {}>'.format(self.id)


class Operator(db.Model):
    __searchable__ = ['body']
    __tablename__ = "operator"

    id = db.Column(db.Integer, primary_key=True, unique=True)

    number_seqs = db.Column(db.Integer)
    consensus_score = db.Column(db.Float(128))
    validated = db.Column(db.Boolean)

        # Operator motif JSON. [ {"base": str, "score": float"}, {...}, ...]
    motif = db.Column(db.Text(4096))

        # All sequences JSON. [ {"prot_id": str, "aligned_seq": str}, {...}, ...]
    aligned_seqs = db.Column(db.Text(4096))

    regulators = db.relationship("Regulator", foreign_keys="Regulator.operator", backref='operator', lazy='dynamic')


    def __repr__(self):
        return '<Operator {}>'.format(self.id)