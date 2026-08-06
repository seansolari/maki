class DatasetError(Exception):
    pass


class DatasetAlreadyRegisteredError(DatasetError):
    pass


class DatasetNotFoundError(DatasetError):
    pass


class DatasetMaterializationError(DatasetError):
    pass


class DatasetChecksumError(DatasetMaterializationError):
    pass


class DatasetTransformationError(DatasetMaterializationError):
    pass