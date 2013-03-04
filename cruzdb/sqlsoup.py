"""

"""

from sqlalchemy import Table, MetaData, join
from sqlalchemy import schema, sql, util
from sqlalchemy.engine.base import Engine
from sqlalchemy.orm import scoped_session, sessionmaker, mapper, \
                            class_mapper, relationship, session,\
                            object_session, attributes
from sqlalchemy.orm.interfaces import MapperExtension, EXT_CONTINUE
from sqlalchemy.sql import expression

__version__ = '0.9.0'
__all__ = ['SQLSoupError', 'SQLSoup', 'SelectableClassType', 'TableClassType', 'Session']

Session = scoped_session(sessionmaker())
"""SQLSoup's default session registry.

This is an instance of :class:`sqlalchemy.orm.scoping.ScopedSession`,
and provides a new :class:`sqlalchemy.orm.session.Session`
object for each application thread which refers to it.

"""

class AutoAdd(MapperExtension):
    def __init__(self, scoped_session):
        self.scoped_session = scoped_session

    def instrument_class(self, mapper, class_):
        class_.__init__ = self._default__init__(mapper)

    def _default__init__(ext, mapper):
        def __init__(self, **kwargs):
            for key, value in kwargs.iteritems():
                setattr(self, key, value)
        return __init__

    def init_instance(self, mapper, class_, oldinit, instance, args, kwargs):
        session = self.scoped_session()
        state = attributes.instance_state(instance)
        session._save_impl(state)
        return EXT_CONTINUE

    def init_failed(self, mapper, class_, oldinit, instance, args, kwargs):
        sess = object_session(instance)
        if sess:
            sess.expunge(instance)
        return EXT_CONTINUE

class SQLSoupError(Exception):
    pass

class ArgumentError(SQLSoupError):
    pass

# metaclass is necessary to expose class methods with getattr, e.g.
# we want to pass db.users.select through to users._mapper.select
class SelectableClassType(type):
    """Represent a SQLSoup mapping to a :class:`sqlalchemy.sql.expression.Selectable`
    construct, such as a table or SELECT statement.
    
    """

    def insert(cls, **kwargs):
        raise SQLSoupError(
            'SQLSoup can only modify mapped Tables (found: %s)' \
              % cls._table.__class__.__name__
        )

    def __clause_element__(cls):
        return cls._table

    def __getattr__(cls, attr):
        if attr == '_query':
            # called during mapper init
            raise AttributeError()
        return getattr(cls._query, attr)

    def __getitem__(cls, key):
        return cls._query[key]

class TableClassType(SelectableClassType):
    """Represent a SQLSoup mapping to a :class:`sqlalchemy.schema.Table`
    construct.
    
    This object is produced automatically when a table-name
    attribute is accessed from a :class:`.SQLSoup` instance.
    
    """
    def insert(cls, **kwargs):
        o = cls()
        o.__dict__.update(kwargs)
        return o

    def relate(cls, propname, *args, **kwargs):
        """Produce a relationship between this mapped table and another
        one. 
        
        This makes usage of SQLAlchemy's :func:`sqlalchemy.orm.relationship`
        construct.
        
        """
        class_mapper(cls)._configure_property(propname, relationship(*args, **kwargs))

def _is_outer_join(selectable):
    if not isinstance(selectable, sql.Join):
        return False
    if selectable.isouter:
        return True
    return _is_outer_join(selectable.left) or _is_outer_join(selectable.right)

def _selectable_name(selectable):
    if isinstance(selectable, sql.Alias):
        return _selectable_name(selectable.element)
    elif isinstance(selectable, sql.Select):
        return ''.join(_selectable_name(s) for s in selectable.froms)
    elif isinstance(selectable, schema.Table):
        return selectable.name
    else:
        x = selectable.__class__.__name__
        if x[0] == '_':
            x = x[1:]
        return x

def _class_for_table(session, engine, selectable, base_cls, mapper_kwargs):
    selectable = expression._clause_element_as_expr(selectable)
    mapname = _selectable_name(selectable)
    # Py2K
    if isinstance(mapname, unicode): 
        engine_encoding = engine.dialect.encoding 
        mapname = mapname.encode(engine_encoding)
    # end Py2K

    if isinstance(selectable, Table):
        klass = TableClassType(mapname, (base_cls,), {})
    else:
        klass = SelectableClassType(mapname, (base_cls,), {})

    def _compare(self, o):
        L = list(self.__class__.c.keys())
        L.sort()
        t1 = [getattr(self, k) for k in L]
        try:
            t2 = [getattr(o, k) for k in L]
        except AttributeError:
            raise TypeError('unable to compare with %s' % o.__class__)
        return t1, t2

    # python2/python3 compatible system of 
    # __cmp__ - __lt__ + __eq__

    def __lt__(self, o):
        t1, t2 = _compare(self, o)
        return t1 < t2

    def __eq__(self, o):
        t1, t2 = _compare(self, o)
        return t1 == t2

    for m in ['__eq__', '__lt__']:
        setattr(klass, m, eval(m))
    klass._table = selectable
    klass.c = expression.ColumnCollection()
    mappr = mapper(klass,
                   selectable,
                   extension=AutoAdd(session),
                   **mapper_kwargs)

    for k in mappr.iterate_properties:
        klass.c[k.key] = k.columns[0]

    klass._query = session.query_property()
    return klass

class SQLSoup(object):
    """Represent an ORM-wrapped database resource."""

    def __init__(self, engine_or_metadata, base=object, session=None):
        """Initialize a new :class:`.SQLSoup`.

        :param engine_or_metadata: a string database URL, :class:`.Engine` 
          or :class:`.MetaData` object to associate with. If the
          argument is a :class:`.MetaData`, it should be *bound*
          to an :class:`.Engine`.
        :param base: a class which will serve as the default class for 
          returned mapped classes.  Defaults to ``object``.
        :param session: a :class:`.ScopedSession` or :class:`.Session` with
          which to associate ORM operations for this :class:`.SQLSoup` instance.
          If ``None``, a :class:`.ScopedSession` that's local to this 
          module is used.

        """

        self.session = session or Session
        self.base=base

        if isinstance(engine_or_metadata, MetaData):
            self._metadata = engine_or_metadata
        elif isinstance(engine_or_metadata, (basestring, Engine)):
            self._metadata = MetaData(engine_or_metadata)
        else:
            raise ArgumentError("invalid engine or metadata argument %r" % 
                                engine_or_metadata)

        self._cache = {}
        self.schema = None

    @property
    def bind(self):
        """The :class:`sqlalchemy.engine.base.Engine` associated with this :class:`.SQLSoup`."""
        return self._metadata.bind

    engine = bind

    def delete(self, instance):
        """Mark an instance as deleted."""

        self.session.delete(instance)

    def execute(self, stmt, **params):
        """Execute a SQL statement.

        The statement may be a string SQL string,
        an :func:`sqlalchemy.sql.expression.select` construct, or a 
        :func:`sqlalchemy.sql.expression.text` 
        construct.

        """
        return self.session.execute(sql.text(stmt, bind=self.bind), **params)

    @property
    def _underlying_session(self):
        if isinstance(self.session, session.Session):
            return self.session
        else:
            return self.session()

    def connection(self):
        """Return the current :class:`sqlalchemy.engine.base.Connection` in use by the current transaction."""

        return self._underlying_session._connection_for_bind(self.bind)

    def flush(self):
        """Flush pending changes to the database.

        See :meth:`sqlalchemy.orm.session.Session.flush`.

        """
        self.session.flush()

    def rollback(self):
        """Rollback the current transction.

        See :meth:`sqlalchemy.orm.session.Session.rollback`.

        """
        self.session.rollback()

    def commit(self):
        """Commit the current transaction.

        See :meth:`sqlalchemy.orm.session.Session.commit`.

        """
        self.session.commit()

    def expunge(self, instance):
        """Remove an instance from the :class:`.Session`.

        See :meth:`sqlalchemy.orm.session.Session.expunge`.

        """
        self.session.expunge(instance)

    def expunge_all(self):
        """Clear all objects from the current :class:`.Session`.

        See :meth:`.Session.expunge_all`.

        """
        self.session.expunge_all()

    def map_to(self, attrname, tablename=None, selectable=None, 
                    schema=None, base=None, mapper_args=util.immutabledict()):
        """Configure a mapping to the given attrname.

        This is the "master" method that can be used to create any 
        configuration.

        :param attrname: String attribute name which will be
          established as an attribute on this :class:.`.SQLSoup`
          instance.
        :param base: a Python class which will be used as the
          base for the mapped class. If ``None``, the "base"
          argument specified by this :class:`.SQLSoup`
          instance's constructor will be used, which defaults to
          ``object``.
        :param mapper_args: Dictionary of arguments which will
          be passed directly to :func:`.orm.mapper`.
        :param tablename: String name of a :class:`.Table` to be
          reflected. If a :class:`.Table` is already available,
          use the ``selectable`` argument. This argument is
          mutually exclusive versus the ``selectable`` argument.
        :param selectable: a :class:`.Table`, :class:`.Join`, or
          :class:`.Select` object which will be mapped. This
          argument is mutually exclusive versus the ``tablename``
          argument.
        :param schema: String schema name to use if the
          ``tablename`` argument is present.


        """
        if attrname in self._cache:
            raise SQLSoupError(
                "Attribute '%s' is already mapped to '%s'" % (
                attrname,
                class_mapper(self._cache[attrname]).mapped_table
            ))

        if tablename is not None:
            if not isinstance(tablename, basestring):
                raise ArgumentError("'tablename' argument must be a string."
                                    )
            if selectable is not None:
                raise ArgumentError("'tablename' and 'selectable' "
                                    "arguments are mutually exclusive")

            selectable = Table(tablename, 
                                        self._metadata, 
                                        autoload=True, 
                                        autoload_with=self.bind, 
                                        schema=schema or self.schema)
        elif schema:
            raise ArgumentError("'tablename' argument is required when "
                                "using 'schema'.")
        elif selectable is not None:
            if not isinstance(selectable, expression.FromClause):
                raise ArgumentError("'selectable' argument must be a "
                                    "table, select, join, or other "
                                    "selectable construct.")
        else:
            raise ArgumentError("'tablename' or 'selectable' argument is "
                                    "required.")

        if not selectable.primary_key.columns and not  \
                             'primary_key' in mapper_args:
            if tablename:
                raise SQLSoupError(
                            "table '%s' does not have a primary "
                            "key defined" % tablename)
            else:
                raise SQLSoupError(
                            "selectable '%s' does not have a primary "
                            "key defined" % selectable)

        mapped_cls = _class_for_table(
            self.session,
            self.engine,
            selectable,
            base or self.base,
            mapper_args
        )
        self._cache[attrname] = mapped_cls
        return mapped_cls


    def map(self, selectable, base=None, **mapper_args):
        """Map a selectable directly.

        The class and its mapping are not cached and will
        be discarded once dereferenced (as of 0.6.6).

        :param selectable: an :func:`.expression.select` construct.
        :param base: a Python class which will be used as the
          base for the mapped class. If ``None``, the "base"
          argument specified by this :class:`.SQLSoup`
          instance's constructor will be used, which defaults to
          ``object``.
        :param mapper_args: Dictionary of arguments which will
          be passed directly to :func:`.orm.mapper`.

        """

        return _class_for_table(
            self.session,
            self.engine,
            selectable,
            base or self.base,
            mapper_args
        )

    def with_labels(self, selectable, base=None, **mapper_args):
        """Map a selectable directly, wrapping the 
        selectable in a subquery with labels.

        The class and its mapping are not cached and will
        be discarded once dereferenced (as of 0.6.6).

        :param selectable: an :func:`.expression.select` construct.
        :param base: a Python class which will be used as the
          base for the mapped class. If ``None``, the "base"
          argument specified by this :class:`.SQLSoup`
          instance's constructor will be used, which defaults to
          ``object``.
        :param mapper_args: Dictionary of arguments which will
          be passed directly to :func:`.orm.mapper`.

        """

        # TODO give meaningful aliases
        return self.map(
                    expression._clause_element_as_expr(selectable).
                            select(use_labels=True).
                            alias('foo'), base=base, **mapper_args)

    def join(self, left, right, onclause=None, isouter=False, 
                base=None, **mapper_args):
        """Create an :func:`.expression.join` and map to it.

        The class and its mapping are not cached and will
        be discarded once dereferenced (as of 0.6.6).

        :param left: a mapped class or table object.
        :param right: a mapped class or table object.
        :param onclause: optional "ON" clause construct..
        :param isouter: if True, the join will be an OUTER join.
        :param base: a Python class which will be used as the
          base for the mapped class. If ``None``, the "base"
          argument specified by this :class:`.SQLSoup`
          instance's constructor will be used, which defaults to
          ``object``.
        :param mapper_args: Dictionary of arguments which will
          be passed directly to :func:`.orm.mapper`.

        """

        j = join(left, right, onclause=onclause, isouter=isouter)
        return self.map(j, base=base, **mapper_args)

    def entity(self, attr, schema=None):
        """Return the named entity from this :class:`.SQLSoup`, or 
        create if not present.

        For more generalized mapping, see :meth:`.map_to`.

        """
        try:
            return self._cache[attr]
        except KeyError, ke:
            return self.map_to(attr, tablename=attr, schema=schema)

    def __getattr__(self, attr):
        return self.entity(attr)

    def __repr__(self):
        return 'SQLSoup(%r)' % self._metadata

