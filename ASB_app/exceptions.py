class ASBException(Exception):
    pass


class ParsingError(ASBException):
    """
    Is raised when there is a custom validation error
    """


class ParserError(ASBException):
    """
    Is raised while incorrect filters parsing
    """
