from sqlalchemy.sql import expression
import sqlalchemy
from sqlalchemy.ext import compiler


class group_concat_distinct_sep(expression.FunctionElement):
    name = "group_concat_distinct_sep"


@compiler.compiles(group_concat_distinct_sep, 'sqlite')
def _group_concat_sep_sqlite(element, compiler, **kw):
    """
    A work-around for an SQLite bug:
    http://sqlite.1065341.n5.nabble.com/GROUP-CONCAT-with-DISTINCT-bug-td99696.html
    """

    if len(element.clauses) == 2:
        separator = compiler.process(element.clauses.clauses[1])
    else:
        separator = ','

    return "REPLACE(REPLACE(GROUP_CONCAT(DISTINCT REPLACE({}, ',', '💩')), ',', {}), '💩', ',')".format(
        compiler.process(element.clauses.clauses[0]),
        separator,
    )
