from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker

#Session = sessionmaker(autoflush=True, expire_on_commit=False,
#       autocommit=False)

#Session = scoped_session(sessionmaker()) #bind=engine))

# http://flask.pocoo.org/snippets/22/

def initialize_sql(engine, Session):
    #Session = sessionmaker(autoflush=True)
    Base = declarative_base(bind=engine)
    Base.metadata.bind = engine
    #Base.metadata.create_all(engine)
    Session.configure(bind=engine)
    session = Session()
    return session, Base
