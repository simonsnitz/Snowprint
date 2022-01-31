from api import db

'''
class Cluster(db.Model):
    __searchable__ = ['body']
    __tablename__ = "cluster"

    id = db.Column(db.Integer, primary_key=True)

        # cluster info
    family = db.Column(db.String(16), index=True)
    cluster_percent_identity = db.Column(db.Integer, index=True)

        # info on representative protein
    rep_EMBL = db.Column(db.String(16), index=True, unique=True)
    rep_genome = db.Column(db.String(16), index=True)
    rep_operon = db.Column(db.PickleType)
    reg_index = db.Column(db.Integer, index=True)
    rep_ligand_SMILES = db.Column(db.String(256), index=True)
    rep_ligand_NAME = db.Column(db.String(256), index=True)
    
    homologs = db.Column(db.PickleType)

        # essentially this is the operator motif
    dna_binding_motif = db.Column(db.PickleType)
        # compounds this cluster is known to recognize. SMILES with associated weights
    ligand_binding_motif = db.Column(db.PickleType)
        # name motif of genes in the associated operons
    gene_context_motif = db.Column(db.PickleType)


    def __repr__(self):
        return '<Cluster {}>'.format(self.id)
'''

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
    start = db.Column(db.Integer)
    stop = db.Column(db.Integer)
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
        # [ {"acc": string, "annotation": string, "direction": string, "start": int, "stop": int}, {...}, ... ]
    operon = db.Column(db.Text(4096))

    regulators = db.relationship("Association", back_populates="operon")


    def __repr__(self):
        return '<Operon {}>'.format(self.id)


class Association(db.Model):
    __tablename__ = "association"

    regulator_id = db.Column(db.ForeignKey(Regulator.id), primary_key=True)
    operon_id = db.Column(db.ForeignKey(Operon.id), primary_key=True)

    reg_index = db.Column(db.Integer)
    reg_type = db.Column(db.Integer)
    regulated_seq = db.Column(db.String(4096))

    operon = db.relationship("Operon", back_populates="regulators")
    regulator = db.relationship("Regulator", back_populates="operons")


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