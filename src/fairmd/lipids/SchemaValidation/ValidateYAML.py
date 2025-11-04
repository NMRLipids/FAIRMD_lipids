import json
import logging
import os

import yaml
from jsonschema import Draft202012Validator, FormatChecker

"""
Tools for validating YAML files against predefined JSON Schemas.

Schemas are stored under the `Schema` subdirectory, and this module
provides helpers to validate both in-memory dicts and on-disk YAML files.
"""


logger = logging.getLogger(__name__)
default_schema_path = os.path.join(os.path.dirname(__file__), "Schema", "info_yml_schema.json")


def validate_info_dict(instance: dict, schema_path: str = default_schema_path):
    """
    Validate an info dict against a schema dict.
    Returns a list of jsonschema.ValidationError objects which is empty with valid.
    """
    with open(schema_path, encoding="utf-8") as f:
        schema = json.load(f)
    validator = Draft202012Validator(schema, format_checker=FormatChecker())
    return list(validator.iter_errors(instance))


def validate_info_file(info_file_path: str, schema_path: str = default_schema_path):
    """
    Validate an info file (YML/YAML) on disk against a JSON schema file.
    Returns a list of ValidationError objects (empty if valid).
    """
    with open(info_file_path, encoding="utf-8") as f:
        instance = yaml.safe_load(f)
    return validate_info_dict(instance, schema_path)
