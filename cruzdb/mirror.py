# from http://www.tylerlesmann.com/2009/apr/27/copying-databases-across-platforms-sqlalchemy/

from sqlalchemy import create_engine, MetaData, Table, select
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.mysql import LONGBLOB, ENUM
from sqlalchemy.types import VARCHAR
import sqlalchemy
import sys
import os


def make_session(connection_string):
    if "///" in connection_string and \
        os.path.exists(connection_string.split("///")[1]) and \
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

def page_query(q, limit=30000):
    offset = 0
    while True:
        r = False
        for elem in q.limit(limit).offset(offset):
           r = True
           yield elem
        offset += limit
        if not r:
            break

def mirror(genome, tables, connection_string):
    #source, sengine = make_session(from_db)
    #smeta = MetaData(bind=sengine)
    destination, dengine = make_session(connection_string)

    for table_name in tables:
        # cause it ot be mapped
        table = genome.table(table_name)
        print 'Processing', table_name
        table = Table(table_name, genome.Base.metadata, autoload=True,
                autoload_with=genome.engine)

        # need to prefix the indexes with the table name to avoid collisions
        for idx in table.indexes:
            idx.name = table_name + "." + idx.name

        for col in table.columns:
            # convert mysql-specific types to varchar
            if isinstance(col.type, (LONGBLOB, ENUM)) \
                    and not 'mysql' in connection_string:
                col.type = VARCHAR()
            elif str(col.type) == "VARCHAR" \
                    and "mysql" in connection_string:
                if col.type.length is None:
                    col.type.length = 95
            if not "mysql" in connection_string:
                if str(col.type).lower().startswith("set("):
                    col.type = VARCHAR(15)

        try:
            table.create(dengine)
        except sqlalchemy.exc.OperationalError:
            pass

        destination.commit()
        ins = table.insert()

        columns = table.columns.keys()
        records = []
        for i, record in enumerate(page_query(getattr(genome, table_name))):
            #print i, "before"
            data = dict(
                (str(column), getattr(record, column)) for column in columns
            )
            records.append(data)
            if 0 == i % 30000:
                destination.execute(ins, records)
                records = []
                print >>sys.stderr, "processing record %i" % i
                destination.commit()
        destination.commit()

    from . import Genome
    return Genome(engine=dengine)

if __name__ == "__main__":
    if False:
        from cruzdb import Genome
        g = Genome('hg18', host="localhost", user="brentp")
        #print g.chromInfo
        #print g.table('chromInfo')

        mirror(g, ['chromInfo', 'cpgIslandExt', 'refGene'], 'sqlite:////tmp/u.db')
