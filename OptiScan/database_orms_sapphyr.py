import sqlalchemy as sq
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

Base = declarative_base()


class RunLineage(Base):
    __tablename__ = "run"
    id = sq.Column(sq.String, primary_key=True)
    organism = sq.Column(sq.String)
    date = sq.Column(sq.String)

class FlowcellLineage(Base):
    __tablename__ = "flowcell"
    id = sq.Column(sq.String, primary_key=True)
    run_id = sq.Column(sq.String, sq.ForeignKey("run.id"))
    fc = sq.Column(sq.String)


class ScanLinage(Base):
    __tablename__ = "scan"
    id = sq.Column(sq.String, primary_key=True)
    run_id = sq.Column(sq.Integer, sq.ForeignKey("run.id"))
    fc_id = sq.Column(sq.Integer, sq.ForeignKey("fc.id"))
    scan = sq.Column(sq.String)
    mem_path = sq.Column(sq.String)


class BankLineage(Base):
    __tablename__ = "bank_array"
    id = sq.Column(sq.String, primary_key=True)
    scan_id = sq.Column(sq.Integer, sq.ForeignKey("scan.id"))
    scan = relationship(ScanLinage)
    col_id = sq.Column(sq.Integer)
    mean_molecule_length = sq.Column(sq.Float)
    molecule_count = sq.Column(sq.Integer)
    minimum_molecule_length = sq.Column(sq.Float)
    abstract_threshold = sq.Column(sq.Float)
    memmap_start = sq.Column(sq.Integer)
    memmap_end = sq.Column(sq.Integer)
    longest_shape = sq.Column(sq.Integer)


class Molecule(Base):
    __tablename__ = "molecule"
    id = sq.Column(sq.String, primary_key=True)
    col_id = sq.Column(sq.Integer, sq.ForeignKey("col_array.id"))
    col = relationship(ColArray) #fix for next db
    index_in_col = sq.Column(sq.Integer)
    # memmap_start = sq.Column(sq.Integer)
    # memmap_end = sq.Column(sq.Integer)
