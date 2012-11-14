# from http://www.tylerlesmann.com/2009/apr/27/copying-databases-across-platforms-sqlalchemy/

from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.mysql import LONGBLOB, ENUM
from sqlalchemy.types import VARCHAR
import sys
import os


def make_session(connection_string):
    if os.path.exists(connection_string.split("///")[1]) and \
        connection_string.startswith("sqlite"):
            print >>sys.stderr, "attempting to add to existing sqlite database"
    engine = create_engine(connection_string, echo=False, convert_unicode=True)
    Session = sessionmaker(bind=engine)
    engine.connect()
    return Session(), engine

def quick_mapper(table):
    Base = declarative_base()
    class GenericMapper(Base):
        __table__ = table
    return GenericMapper

def mirror(genome, tables, connection_string):
    #source, sengine = make_session(from_db)
    #smeta = MetaData(bind=sengine)
    destination, dengine = make_session(connection_string)

    for table_name in tables:
        table = genome.table(table_name)
        print 'Processing', table_name
        table = Table(table_name, genome.Base.metadata, autoload=True)

        # need to prefix the indexes with the table name to avoid collisions
        for idx in table.indexes:
            idx.name = table_name + "." + idx.name

        for col in table.columns:
            # convert mysql-specific types to varchar
            if isinstance(col.type, (LONGBLOB, ENUM)):
                col.type = VARCHAR()

        table.metadata.create_all(dengine)
        destination.commit()
        #[c.index for c in table.columns]
        NewRecord = quick_mapper(table)
        columns = table.columns.keys()
        for record in genome.session.query(table):
            data = dict(
                (str(column), getattr(record, column)) for column in columns
            )
            destination.merge(NewRecord(**data))
        destination.commit()
    destination.close()

if __name__ == "__main__":
    from cruzdb import Genome
    g = Genome('hg18', host="localhost", user="brentp")
    #print g.chromInfo
    #print g.table('chromInfo')

    mirror(g, ['chromInfo', 'cpgIslandExt', 'refGene'], 'sqlite:////tmp/u.db')
