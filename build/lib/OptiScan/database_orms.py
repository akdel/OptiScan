import sqlalchemy as sq
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker


Base = declarative_base()


class RunLineage(Base):
    __tablename__ = "run"
    id = sq.Column(sq.String, primary_key=True)
    organism = sq.Column(sq.String)
    date = sq.Column(sq.String)


class ScanLinage(Base):
    __tablename__ = "scan"
    id = sq.Column(sq.Integer, primary_key=True)
    run_id = sq.Column(sq.Integer, sq.ForeignKey("run.id"))
    scan = sq.Column(sq.Integer)
    bnx_path = sq.Column(sq.String)
    mol_path = sq.Column(sq.String)
    stitch_path = sq.Column(sq.String)
    lab_path = sq.Column(sq.String)
    tiff_path = sq.Column(sq.String)
    mem_path = sq.Column(sq.String)


class ColArray(Base):
    __tablename__ = "col_array"
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
