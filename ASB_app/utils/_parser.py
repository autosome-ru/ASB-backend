from collections import OrderedDict
from more_itertools import peekable
import operator
import re

from ASB_app import exceptions


class Token:
    __slots__ = ['type_', 'value', 'state']

    def __init__(self, type_, value, state=0):
        self.type_ = type_
        self.value = value
        self.state = state

    def __repr__(self):
        return 'Token<{0.type_}; "{0.value}">'.format(self)


class Rule:
    __slots__ = ['token_type', 'token_len', 'func']

    def __init__(self, token_type, token_len, func):
        self.token_type = token_type
        self.token_len = token_len
        self.func = func

    def apply(self, stack):
        args = [stack.pop() for _ in range(self.token_len)]
        args.reverse()
        prev_state = args[0].state
        stack.append(Token(self.token_type, self.func(*self._values(args)), prev_state))
        return prev_state

    @staticmethod
    def _values(args):
        return [arg.value for arg in args]


class SyntaxNode:
    """
    Base element to store abstract syntax tree.
    """
    def __init__(self, *args):
        self.args = args


class Parser:
    lr_actions = {}
    lr_goto = {}
    rules = []
    token_patterns = {}

    def parse(self, tokens_iter):
        try:
            return self._parse(tokens_iter)
        except StopIteration:
            raise exceptions.ParserError('Unexpected end of input. Use Parser.tokens_stream generator')
        except KeyError as e:
            raise exceptions.ParserError('Unexpected token: {}'.format(e))
        except IndexError:
            raise exceptions.ParserError('LR tables are corrupted')

    def tokens_stream(self, line):
        while line:
            for regexp in self.token_patterns:
                match = regexp.match(line)
                if match:
                    type_ = self.token_patterns[regexp]
                    if type_:
                        yield Token(type_, match.group(1))
                    line = line[match.end():]
                    break
            else:
                # No one of patterns matched remain part of line.
                raise exceptions.ParserError('Unknown pattern near {}'.format(line))
        yield Token('$end', None)

    def _parse(self, tokens_iter):
        self.state = 0
        self.stack = []
        tokens = peekable(tokens_iter)
        while True:
            token = next(tokens)
            token.state = self.state
            next_state = self.lr_actions[self.state][token.type_]
            if next_state > 0:
                # Just push new token into the stack, change state and read next token
                self.stack.append(token)
            elif next_state < 0:
                # Apply one of rules to existing tokens
                tokens.prepend(token)  # put token back into the stream
                rule = self.rules[-next_state]
                prev_state = rule.apply(self.stack)
                next_state = self.lr_goto[prev_state][rule.token_type]
            else:
                assert token.type_ == '$end', 'LR tables are corrupted: only $end token may terminate parsing'
                break
            self.state = next_state

        assert len(self.stack) == 1
        return self.stack[0].value


class FilterParser(Parser):
    # Do not use this dictionary to compare python values.
    # Some of its values may create sqlalchemy filters even for two scalars, like int.
    mapping_operators = {
        'EQ': operator.eq, '=': operator.eq, '==': operator.eq,
        'GT': operator.gt, '>': operator.gt,
        'GE': operator.ge, '>=': operator.ge,
        'LT': operator.lt, '<': operator.lt,
        'LE': operator.le, '<=': operator.le,
        'NE': operator.ne, '!=': operator.ne}

    class NodeAnd(SyntaxNode):
        def __init__(self, left, _and_token, right):
            super().__init__(left, right)

    class NodeOr(SyntaxNode):
        def __init__(self, left, _or_token, right):
            super().__init__(left, right)

    class NodeSimpleFilter(SyntaxNode):
        pass

    class NodeAggregationFilter(SyntaxNode):
        def __init__(self, aggr, _lb, relation, _sc, expr, _rb, op, const):
            super().__init__(aggr, relation, expr, op, const)

    lr_actions = {
        0: {'FIELD': 2, 'LBRACKET': 3, 'AGGREGATE': 4},
        3: {'FIELD': 2, 'LBRACKET': 3, 'AGGREGATE': 4},
        5: {'FIELD': 2, 'LBRACKET': 3, 'AGGREGATE': 4},
        6: {'FIELD': 2, 'LBRACKET': 3, 'AGGREGATE': 4},
        9: {'FIELD': 14},
        15: {'FIELD': 2, 'LBRACKET': 3, 'AGGREGATE': 4},
        4: {'LBRACKET': 9},
        1: {'$end': 0, 'AND': 5, 'OR': 6},
        10: {'$end': -4, 'AND': -4, 'OR': -4, 'RBRACKET': -4},
        11: {'$end': -5, 'AND': 5, 'OR': -5, 'RBRACKET': -5},
        12: {'$end': -1, 'AND': -1, 'OR': -1, 'RBRACKET': -1},
        13: {'$end': -2, 'AND': -2, 'OR': -2, 'RBRACKET': -2},
        19: {'$end': -3, 'AND': -3, 'OR': -3, 'RBRACKET': -3},
        8: {'AND': 5, 'OR': 6, 'RBRACKET': 13},
        16: {'AND': 5, 'OR': 6, 'RBRACKET': 17},
        2: {'OPERATOR': 7},
        17: {'OPERATOR': 18},
        7: {'CONST': 12},
        18: {'CONST': 19},
        14: {';': 15}
    }
    lr_goto = {
        0: {'expression': 1},
        3: {'expression': 8},
        5: {'expression': 10},
        6: {'expression': 11},
        15: {'expression': 16}
    }
    # Note: this object should be sorted since two letter operators (>=) should not be parsed as (>) + (=),
    #  while regexp (=|>|>=) does exactly this error, and *unordered* dict {'>=': 'OP', '>': 'OP'} may do it too.
    # In Python 3.7+ build-in dict already sorted, but it is important here.
    token_patterns = OrderedDict({
        # Regexp object -> Token type or None to ignore.
        # Each regexp must contain one group. The first group (not the whole match!) is used as token value.
        re.compile(r'(==|>=|<=|!=)'): 'OPERATOR',
        re.compile(r'([<=>])'): 'OPERATOR',
        # Match only whole word. It is one of ways to avoid problem with (le)ct_id
        re.compile(r'\b(EQ|GT|GE|LT|LE|NE)\b'): 'OPERATOR',  # option: re.IGNORECASE
        re.compile(r'(ANY|ALL|HAS|SUM|COUNT)'): 'AGGREGATE',
        re.compile(r'(\()'): 'LBRACKET',
        re.compile(r'(\))'): 'RBRACKET',
        re.compile(r'(;)'): ';',
        re.compile(r'"([^"]*)"'): 'CONST',
        re.compile(r'([a-z_\d]+)'): 'FIELD',
        re.compile(r'(AND|,)'): 'AND',
        re.compile(r'(OR)'): 'OR',
        re.compile(r'(\s+)'): None
    })

    rules = [
        Rule("S'", 1, None),
        Rule('expression', 3, NodeSimpleFilter),
        Rule('expression', 3, lambda lb, filter_, rb: filter_),
        Rule('expression', 8, NodeAggregationFilter),
        Rule('expression', 3, NodeAnd),
        Rule('expression', 3, NodeOr),
    ]


assert all(
    any(regexp.match(op) for (regexp, value) in FilterParser.token_patterns.items() if value == 'OPERATOR')
    for op in FilterParser.mapping_operators
), "Not all operators can be parsed by regular expressions!"