class ASBException(Exception):
    pass


class ParsingError(ASBException):
    """
    Is raised when there is a custom validation error
    """
