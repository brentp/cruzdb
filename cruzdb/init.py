from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker


Session = None

def initialize_sql(engine):
    global Session
    Session = sessionmaker(autoflush=True)
    #if Session.registry.has(): return Session
    #Session = scoped_session(sessionmaker(bind=engine))
    session = Session(bind=engine)
    Base = declarative_base()
    Base.metadata.bind = engine
    Base.metadata.create_all(engine)
    return session, Base
