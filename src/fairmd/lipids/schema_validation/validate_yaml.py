"""
Tools for validating YAML files against predefined JSON Schemas.

This module can be run either as a module or as a script.

Module usage:
-------------
.. code-block:: bash

    python -m fairmd.lipids.schema_validation.validate_yaml --schema readme README.yaml

Script usage:
-------------
.. code-block:: bash

    # Validate a README.yaml file
    python validate_yaml.py --schema readme README.yaml
    python validate_yaml.py --schema info info.yml
    # Multiple files can be validated at once
    python validate_yaml.py README.yaml other/README.yaml

"""

import argparse
import datetime
import json
import logging
import os
import sys
from typing import Literal

import yaml
from jsonschema import Draft7Validator, FormatChecker, SchemaError, ValidationError

logger = logging.getLogger(__name__)

default_info_schema_path = os.path.join(os.path.dirname(__file__), "schema", "info_yml_schema.json")
default_readme_yaml_schema_path = os.path.join(os.path.dirname(__file__), "schema", "readme_yaml_schema.json")


schema_type_options = Literal["info", "readme"]


def _filter_yaml_date_string_type_errors(errors: list[ValidationError]) -> list[ValidationError]:
    """
    PyYAML may parse unquoted ISO dates (e.g. 2021-02-23) into datetime.date.
    JSON Schema expects a string (often with format: date).
    This filters out only that specific type mismatch.
    """
    out: list[ValidationError] = []
    for e in errors:
        if e.validator == "type" and "string" in e.validator_value:
            inst = getattr(e, "instance", None)
            if isinstance(inst, (datetime.date, datetime.datetime)):
                continue
        out.append(e)
    return out


def validate_info_dict(instance: dict, schema_path: str = default_info_schema_path):
    """
    Validate an info dict against a schema dict.
    Returns a list of jsonschema.ValidationError objects which is empty with valid.
    """
    with open(schema_path, encoding="utf-8") as f:
        schema = json.load(f)
    validator = Draft7Validator(schema, format_checker=FormatChecker())
    errors = list(validator.iter_errors(instance))
    return _filter_yaml_date_string_type_errors(errors)


def validate_info_file(info_file_path: str, schema_path: str = default_info_schema_path):
    """
    Validate an info file (YML/YAML) on disk against a JSON schema file.
    Returns a list of ValidationError objects (empty if valid).
    """
    with open(info_file_path, encoding="utf-8") as f:
        instance = yaml.safe_load(f)
        if not isinstance(instance, dict) or not instance:
            raise ValueError("YAML did not contain a non-empty mapping")
    return validate_info_dict(instance, schema_path)


def validate_readme_dict(instance: dict, schema_path: str = default_readme_yaml_schema_path):
    """
    Validate a README.yaml dict against the README JSON schema.

    Returns a list of jsonschema.ValidationError objects (empty if valid).
    """
    with open(schema_path, encoding="utf-8") as f:
        schema = json.load(f)
    validator = Draft7Validator(schema, format_checker=FormatChecker())
    errors = list(validator.iter_errors(instance))
    return _filter_yaml_date_string_type_errors(errors)


def validate_readme_file(readme_file_path: str, schema_path: str = default_readme_yaml_schema_path):
    """
    Validate a README.yaml file on disk against the README JSON schema.

    Returns a list of ValidationError objects (empty if valid).
    """
    with open(readme_file_path, encoding="utf-8") as f:
        instance = yaml.safe_load(f)
        if not isinstance(instance, dict) or not instance:
            raise ValueError("YAML did not contain a non-empty mapping")
    return validate_readme_dict(instance, schema_path)


def run_file(path: str, schema_type: schema_type_options) -> None:
    """
    Validate a single YAML file against the selected schema.

    On success, logs "OK: <path>" and returns normally.
    On failure, logs detailed schema errors and raises an exception.

    Raises:
        FileNotFoundError:
            If the file does not exist or is not a regular file.

        RuntimeError:
            If the file was successfully read but failed schema validation.

        ValueError:
            If the YAML does not contain a non-empty mapping (invalid structure).

        OSError, yaml.YAMLError, json.JSONDecodeError, SchemaError:
            For I/O errors, invalid YAML, invalid JSON schema, or other runtime failures.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(path)

    if schema_type == "info":
        errors = validate_info_file(path)
    else:
        errors = validate_readme_file(path)

    if not errors:
        logger.info("OK: %s", path)
        return

    logger.error("INVALID: %s", path)
    for err in errors:
        keys = ".".join(str(p) for p in err.path) if err.path else "<root>"
        logger.error("  -> %s (at %s)", err.message, keys)

    raise RuntimeError("Schema validation failed")


def main() -> int:
    """
    Command-line entry point for YAML schema validation.

    Validates one or more YAML files against either the FAIRMD
    info.yml schema or the README.yaml schema.

    :returns: Process exit code:
              0 = all files valid
              1 = at least one file failed schema validation
              2 = at least one file was missing or not a regular file
              3 = at least one file could not be read, parsed, or validated due to YAML, JSON, or runtime errors
    """
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(description="Validate YAML files against FAIRMD info / README schemas.")
    parser.add_argument(
        "--schema",
        "-s",
        choices=["info", "readme"],
        default="readme",
        help="Schema to use: 'readme' (default) or 'info'.",
    )
    parser.add_argument(
        "files",
        nargs="+",
        help="YAML files to validate.",
    )

    args = parser.parse_args()

    exit_code = 0
    for f in args.files:
        path = os.path.normpath(f)
        try:
            run_file(path, args.schema)
        except RuntimeError:
            exit_code = max(exit_code, 1)
        except FileNotFoundError:
            logger.error("File not found (or not a file): %s", path)
            exit_code = max(exit_code, 2)
        except (ValueError, OSError, yaml.YAMLError, json.JSONDecodeError, SchemaError) as e:
            logger.error("ERROR: %s", path)
            logger.error("  -> %s: %s", type(e).__name__, e)
            exit_code = max(exit_code, 3)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
