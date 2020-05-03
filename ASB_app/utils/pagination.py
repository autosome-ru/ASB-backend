__all__ = ['PaginationMixin', 'is_field', 'is_relationship', 'get_entity_from_relationship']

import re

from sqlalchemy.sql.elements import True_
from sqlalchemy.inspection import inspect
from sqlalchemy import desc
from flask_sqlalchemy import Model

from ASB_app import exceptions, session
from ._parser import FilterParser
from ..models import TranscriptionFactor, CellLine, SNP, TranscriptionFactorSNP, CellLineSNP


def is_field(entity_cls, attr):
    """
    Shows if attribute attr represents field to which filtering can be applied.
    Raises an exception if arguments has wrong type.
    :param type entity_cls: entity class
    :param str attr: attribute name
    :return bool:
    """
    if attr == '__mapper__':
        return False
    mapper = inspect(entity_cls)
    return attr in mapper.all_orm_descriptors and attr not in mapper.relationships


def is_relationship(entity_cls, attr):
    """
    Shows if attribute attr represents relationship to any other entity.
    Raises an exception if arguments has wrong type.
    :param type entity_cls: entity class
    :param str attr: relationship name
    :return bool:
    """
    return attr in inspect(entity_cls).relationships


def get_entity_from_relationship(entity_cls, attr):
    """
    Return entity on which relationship named attr shows.
    :param type entity_cls: entity class
    :param str attr: relationship name
    :return type: another entity
    """
    if not is_relationship(entity_cls, attr):
        raise exceptions.ParserError('"{}" is not a relationship of {}'.format(attr, entity_cls))
    return getattr(entity_cls, attr).mapper.class_


class PaginationMixin(FilterParser):
    BaseEntity = None

    def __init__(self, base_entity_cls=None):
        """
        Создает объект класса PaginationMixin с установленным значением базовой сущности.
        Созданный таким образом объект можно использовать для доступа к функции paginate без создания подклассов.
        """
        if base_entity_cls is not None:
            self.BaseEntity = base_entity_cls

    @staticmethod
    def _mapping_const(const):
        # FIXME: Временное решение для парсинга фильтров ' == "NULL"' и ' != "NULL"'
        return None if const == "NULL" else const

    def _check_entity_type(self):
        if not issubclass(self.BaseEntity, Model):
            raise ValueError('PaginationMixin requires set BaseEntity')

    def _evaluate_tree(self, syntax_tree):
        if isinstance(syntax_tree, FilterParser.NodeAnd):
            return self._evaluate_tree(syntax_tree.args[0]) & self._evaluate_tree(syntax_tree.args[1])
        elif isinstance(syntax_tree, FilterParser.NodeOr):
            return self._evaluate_tree(syntax_tree.args[0]) | self._evaluate_tree(syntax_tree.args[1])
        elif isinstance(syntax_tree, FilterParser.NodeSimpleFilter):
            return self._construct_simple_filter(*syntax_tree.args)
        elif isinstance(syntax_tree, FilterParser.NodeAggregationFilter):
            return self._construct_aggregation_filter(*syntax_tree.args)
        else:
            raise TypeError('Syntax tree must be instance of subclass of SyntaxNode, got {}'.format(syntax_tree))

    def _construct_simple_filter(self, field, op, const):
        if not is_field(self.BaseEntity, field):
            raise exceptions.ParserError('"{}" is not a field of {}'.format(field, self.BaseEntity))
        return self.mapping_operators[op](getattr(self.BaseEntity, field), self._mapping_const(const))

    def _construct_aggregation_filter(self, aggr, relation, expr, op, const):
        related_entity = get_entity_from_relationship(self.BaseEntity, relation)
        sub_expr = PaginationMixin(related_entity)._evaluate_tree(expr)
        op = self.mapping_operators[op]
        relation = getattr(self.BaseEntity, relation)
        # One has to use HAS when related entity is a 'parent' and ANY/ALL when 'child'
        # Otherwise raises sqlalchemy.exc.InvalidRequestError
        if aggr == 'ANY':
            return op(relation.any(sub_expr), self._mapping_const(const))
        elif aggr == 'ALL':
            return ~op(relation.any(~sub_expr), self._mapping_const(const))
        elif aggr == 'HAS':
            return op(relation.has(sub_expr), self._mapping_const(const))
        else:
            raise NotImplementedError

    def parse_filters(self, filters_str):
        """
        Преобразовывает строку с условиями, разделенными запятой, в список условий sqlalchemy.
        Выражения, не имеющие вид "column operator value" без запятых или с неизвестными значениями
         для column или operator, игнорируются.
        :param filters_str:
        :return sqlalchemy.sql.elements.ClauseElement: SQLAlchemy object, representing all condition in filters_str.
        """
        if not filters_str:
            return True_()
        syntax_tree = self.parse(self.tokens_stream(filters_str))
        return self._evaluate_tree(syntax_tree)

    def parse_order_clauses(self, sorting_str):
        """
        Преобразовывает строку с критериями сортировки, разделенными ::, в список критериев sqlalchemy.
        Выражения, не имеющие вид "[-]column" неизвестными значениями для column, игнорируются.
        :param sorting_str:
        :return:
        """
        if not sorting_str:
            return [], [], []
        order_clauses = []
        additional_columns = []
        additional_join_tuples = []
        for criterion in sorting_str.strip().split('::'):
            descending = False
            if criterion.startswith('-'):
                descending = True
                criterion = criterion[1:]

            if is_field(self.BaseEntity, criterion):
                variable = getattr(self.BaseEntity, criterion)
                order_clauses.append(variable.desc() if descending else variable)
            elif re.search('^(CL|TF)@log_p_value_(ref|alt)@[^@]+$', criterion) is not None:
                if self.BaseEntity != SNP:
                    continue
                what_for, field, name = criterion.split('@')
                aggregation_class = {'TF': TranscriptionFactor, 'CL': CellLine}[what_for]
                aggregation_entity = aggregation_class.query.filter_by(name=name).one_or_none()
                if not aggregation_entity:
                    continue
                aggregation_id = getattr(aggregation_entity, {'TF': 'tf_id', 'CL': 'cl_id'}[what_for])

                aggregated_snp_class = {'TF': TranscriptionFactorSNP, 'CL': CellLineSNP}[what_for]
                additional_join_tuples.append(
                    (
                        aggregated_snp_class,
                        (getattr(
                            aggregated_snp_class,
                            {'TF': 'tf_id', 'CL': 'cl_id'}[what_for]
                        ) == aggregation_id) &
                        (aggregated_snp_class.chromosome == SNP.chromosome) &
                        (aggregated_snp_class.position == SNP.position) &
                        (aggregated_snp_class.alt == SNP.alt)
                    )
                )
                label = '_'.join(map(str, [what_for, field, aggregation_id]))
                additional_columns.append(getattr(aggregated_snp_class, field).label(label))
                order_clauses.append(desc(label) if descending else label)

        return order_clauses, additional_columns, additional_join_tuples

    def paginate(self, args, extra_filters=()):
        """
        Returns list of self.BaseEntity objects taking into account the parameters passed in args.
        :param args: dictionary with the following keys:
                     size - number of objects should be returned. Interpreted as capacity of one page.
                            No limitations if size is not provided or equals to zero.
                     offset - number of objects should be skipped before pagination.
                     page - number of page to return. Ignored, if size is not positive integer.
                     filter - string specifying filtering rules. See `parse_filters` method.
                     order_by - string specifying order of objects in the selection. See `parse_order_clauses` method.
        :param extra_filters: list of additional filters in sqlalchemy format, like 'User.id == 4'.
                              Use this parameter to restrict access to objects without changing filtering string.
        :return: list of BaseEntity objects.
        """
        self._check_entity_type()

        page = args.get('page')
        size = args.get('size')
        offset = args.get('offset')
        filters_str = args['filter']
        sorting_str = args.get('order_by')

        parsed_filters = self.parse_filters(filters_str)
        order_clauses, additional_columns, additional_join_tuples = self.parse_order_clauses(sorting_str)

        query = self.BaseEntity.query

        for model, join_criterion in additional_join_tuples:
            query = query.join(model, join_criterion)

        query = query.filter(parsed_filters, *extra_filters)

        for column in additional_columns:
            query = query.add_column(column)

        query = query.order_by(*order_clauses)

        if size:
            start = offset + size * (page - 1)
            query = query.slice(start, start + size)
        elif offset:
            query = query.offset(offset)

        if additional_columns:
            return [item[0] for item in query]
        else:
            return query.all()

    def items_count(self, filters_str="", extra_filters=()):
        """
        Return amount of self.BaseEntity objects in the database satisfying given filters.
        Supports both filters_str (see `parse_filters` method) and prepared sqlalchemy conditions.
        :param str filters_str: string specifying filtering rules. See `parse_filters` method.
        :param extra_filters: list of additional filters in sqlalchemy format, like 'User.id == 4'.
                              Use this parameter to restrict access to objects without changing filtering string.
        :return int: number of BaseEntity objects.
        """
        self._check_entity_type()
        parsed_filters = self.parse_filters(filters_str)
        return self.BaseEntity.query.filter(parsed_filters, *extra_filters).count()
