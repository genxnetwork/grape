import os
import sys
import json


class ReferenceDirectory:
    @staticmethod
    def _get_valid_content_filepath():
        """
        Path to the JSON file with valid reference directory content. The file contains
        relative paths of the reference files along with their content size. This data is
        used to validate downloaded reference data instead of MD5 hash, since hashing
        takes significant amount of time.
        """

        module = sys.modules[ReferenceDirectory.__module__]
        return os.path.join(
            os.path.dirname(module.__file__),
            'reference_directory_content.json'
        )

    @staticmethod
    def _get_content(reference_directory_path):
        """
        Get a content structure of the reference directory.
        Return a dictionary of relative paths and the files content size.
        """

        content = {}
        for root, _, filenames in os.walk(reference_directory_path):
            for filename in filenames:
                filepath = os.path.join(root, filename)
                relative_path = os.path.relpath(filepath, reference_directory_path)
                content[relative_path] = os.path.getsize(filepath)

        return content

    def __init__(self, path):
        self.path = path
        self.content = self._get_content(path)

    def is_valid(self):
        if not os.path.exists(self.path):
            return False

        with open(self._get_valid_content_filepath(), 'r') as dump_file:
            content = json.load(dump_file)

        return content == self.content

    def to_json(self, filepath):
        with open(filepath, 'w') as dump_file:
            json.dump(self.content, dump_file)
