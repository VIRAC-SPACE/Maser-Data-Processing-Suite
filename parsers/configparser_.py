"""
parse configure files
"""
import configparser


class Singleton(type):
    """
    singleton base class
    """
    __instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls.__instances:
            cls.__instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        else:
            cls.__instances[cls].__init__(*args, **kwargs)
        return cls.__instances[cls]


class ConfigParser(metaclass=Singleton):
    """
    singleton config parser
    """

    def __init__(self, config_file_path):
        if config_file_path is not None:
            self._config_file_path = config_file_path
            self._config = configparser.RawConfigParser()
            self._config.read(self._config_file_path)

    def get_config(self, section, key):
        """

        :param section: configuration file section
        :param key: key of section
        :return: configuration file value from section with key param key
        """
        return self._config.get(section, key)

    def get_items(self, section):
        """

        :param section: section of configuration file
        :return: all items from section of configuration file
        """
        return dict(self._config.items(section))
