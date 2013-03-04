import sqlsoup
from sqlalchemy import Table, util


class Genome(sqlsoup.SQLSoup):

    def map_to(self, attrname, tablename=None, selectable=None,
                    schema=None, base=None, mapper_args=util.immutabledict()):
        tbl = Table(tablename, self._metadata, autoload=True,
                 autoload_with=self.bind, schema=schema or self.schema)

        # make a fake primary key
        pids = [x for x in tbl.columns if x.name in ('chrom', 'chromStart', 'name',
                                                          'txStart', 'kgID')
                                       or x.primary_key or x.unique]
        models = __import__("cruzdb.models", globals(), locals(), [], -1).models
        try:
            base = getattr(models, tablename)
        except AttributeError:
            base = models.Feature

        mapper_args = dict(mapper_args)
        mapper_args['primary_key'] = pids
        return sqlsoup.SQLSoup.map_to(self, attrname, tablename, selectable,
                                       schema, base=base, mapper_args=mapper_args)

if __name__ == "__main__":
    db = Genome('mysql://genome@genome-mysql.cse.ucsc.edu/hg19')

    print db.cpgIslandExt.first()
    print db.refGene.first()

