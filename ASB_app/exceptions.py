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


class ReleaseNotFound(ASBException):
    """
    Is raised when specified release version is not found
    """


class FileNotProcessed(ASBException):
    """
    Is raised when trying to download ANANASTra ticket ptocessing result for an unprocessed ticket
    """