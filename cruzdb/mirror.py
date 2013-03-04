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
    Session = sessionmaker(bind=engine, autoflush=False, expire_on_commit=False,
            autocommit=False)
    engine.connect()
    return Session(), engine

def page_query(q, session, limit=8000):
    #if q.count() < 80000:
    #    for rn q.all()
    offset = 0
    while True:
        elem = None
        for elem in session.execute(q.offset(offset).limit(limit)):
            yield elem
        offset += limit
        if elem is None:
            break

def set_table(genome, table, table_name, connection_string, metadata):
    """
    alter the table to work between different
    dialects
    """
    table = Table(table_name, genome._metadata, autoload=True,
                    autoload_with=genome.bind, extend_existing=True)

    #print "\t".join([c.name for c in table.columns])
    # need to prefix the indexes with the table name to avoid collisions
    for i, idx in enumerate(table.indexes):
        idx.name = table_name + "." + idx.name + "_ix" + str(i)

    cols = []
    for i, col in enumerate(table.columns):
        # convert mysql-specific types to varchar
        #print col.name, col.type, isinstance(col.type, ENUM)
        if isinstance(col.type, (LONGBLOB, ENUM)):

            if 'sqlite' in connection_string:
                col.type = VARCHAR()
            elif 'postgres' in connection_string:
                if isinstance(col.type, ENUM):
                    #print dir(col)
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

    table = Table(table_name, genome._metadata, *cols,
            autoload_replace=True, extend_existing=True)

    return table

def mirror(genome, tables, connection_string):
    destination, dengine = make_session(connection_string)
    dmeta = MetaData(bind=dengine)

    orig_counts = []
    for table_name in tables:
        # cause it ot be mapped
        table = getattr(genome, table_name)._table
        print >>sys.stderr, 'Mirroring', table_name

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
        table_obj = getattr(genome, table_name)._table
        t = getattr(genome, table_name)
        for ii, record in enumerate(page_query(table_obj.select(), t.session)):
            data = dict(
                (str(column), getattr(record, column)) for column in columns
            )
            records.append(data)
            if ii % 20000 == 0 and ii > 0:
                destination.execute(ins, records)
                print >>sys.stderr, "processing record %i" % ii
                destination.commit()
                records = []
        destination.execute(ins, records)
        destination.commit()
        orig_counts.append(getattr(genome, table_name).count())

    destination, dengine = make_session(connection_string)
    from . import Genome
    newg = Genome(connection_string)
    new_counts = [getattr(newg, table_name).count() for table_name in tables]
    for tbl, oc, nc in zip(tables, orig_counts, new_counts):
        if oc != nc: print >>sys.stderr, "ERROR: mirrored table '%s' has %i \
            rows while the original had %i" % (tbl, nc, oc)
    return newg

if __name__ == "__main__":
    if True:
        from cruzdb import Genome
        g = Genome('hg18')

        mirror(g, ['chromInfo'], 'sqlite:////tmp/u.db')
