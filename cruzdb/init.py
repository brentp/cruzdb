from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker


#Session = scoped_session(sessionmaker()) #bind=engine))
Session = sessionmaker(autoflush=False, expire_on_commit=False)

def initialize_sql(engine):
    global Session
    #Session = sessionmaker(autoflush=True)
    session = Session(bind=engine)
    Base = declarative_base()
    Base.metadata.bind = engine
    Base.metadata.create_all(engine)
    return session, Base
