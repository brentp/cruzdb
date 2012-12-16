# from http://www.tylerlesmann.com/2009/apr/27/copying-databases-across-platforms-sqlalchemy/

from sqlalchemy import create_engine, Table, select, Enum, Column, MetaData
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.mysql import LONGBLOB, ENUM
from sqlalchemy.dialects.postgresql import ENUM as PG_ENUM
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

def page_query(q, limit=70000):
    offset = 0
    while True:
        elem = None
        for elem in q.limit(limit).offset(offset):
           yield elem
        offset += limit
        if elem is None:
            break

def set_table(genome, table, table_name, connection_string, metadata):
    """
    alter the table to work between different
    dialects
    """
    table = Table(table_name, genome.Base.metadata, autoload=True,
                    autoload_with=genome.engine, extend_existing=True)

    print "\t".join([c.name for c in table.columns])
    # need to prefix the indexes with the table name to avoid collisions
    for i, idx in enumerate(table.indexes):
        idx.name = table_name + "." + idx.name + "_ix" + str(i)

    cols = []
    for i, col in enumerate(table.columns):
        # convert mysql-specific types to varchar
        print col.name, col.type, isinstance(col.type, ENUM)
        if isinstance(col.type, (LONGBLOB, ENUM)):

            if 'sqlite' in connection_string:
                col.type = VARCHAR()
            elif 'postgres' in connection_string:
                if isinstance(col.type, ENUM):
                    print dir(col)
                    col.type = PG_ENUM(*col.type.enums, name=col.name,
                        create_type=True)
                else:
                    col.type = VARCHAR()
        elif str(col.type) == "VARCHAR" \
                and ("mysql" in connection_string \
                or "postgres" in connection_string):
            if col.type.length is None:
                col.type.length = 48 if col.name != "description" else None
        if not "mysql" in connection_string:
            if str(col.type).lower().startswith("set("):
                col.type = VARCHAR(15)
        cols.append(col)
    table = Table(table_name, genome.Base.metadata, *cols,
            autoload_replace=True, extend_existing=True)

    return table

def mirror(genome, tables, connection_string):
    destination, dengine = make_session(connection_string)
    dmeta = MetaData(bind=dengine)

    for table_name in tables:
        # cause it ot be mapped
        table = genome.table(table_name)
        print >>sys.stderr, 'Processing', table_name

        table = set_table(genome, table, table_name,
                connection_string, dmeta)
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
        destination.execute(ins, records)
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
