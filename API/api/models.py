from api import db

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