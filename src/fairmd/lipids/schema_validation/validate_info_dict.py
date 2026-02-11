import copy

from fairmd.lipids.molecules import lipids_set, molecule_ff_set, solubles_set
from fairmd.lipids.schema_validation.engines import software_dict


class YamlBadConfigError(Exception):
    """Custom Exception class for parsing the yaml configuration"""


def parse_valid_config_settings(info_yaml: dict, logger) -> tuple[dict, list[str]]:
    """
    Parse, validate and update dict entries from yaml configuration file.

    Args:
        info_yaml (dict): info.yaml of database to add
    Raises:
        KeyError: Missing required key in info.yaml
        YamlBadConfigException: Incorrect or incompatible configuration
    Returns:
        dict: updated sim dict
        list[str]: list of filenames to download
    """
    sim = copy.deepcopy(info_yaml)  # mutable objects are called by reference in Python

    # STEP 1 - check supported simulation software
    if "SOFTWARE" not in sim:
        msg = "'SOFTWARE' Parameter missing in yaml"
        raise KeyError(msg)

    if sim["SOFTWARE"].upper() in software_dict:
        logger.info(f"Simulation uses supported software '{sim['SOFTWARE'].upper()}'")
    else:
        msg = f"Simulation uses unsupported software '{sim['SOFTWARE'].upper()}'"
        raise YamlBadConfigError(msg)

    software_sim = software_dict[sim["SOFTWARE"].upper()]  # related to dicts in this file

    # STEP 2 - check required keys defined by sim software used
    software_required_keys = [k for k, v in software_sim.items() if v["REQUIRED"]]

    # are ALL required keys are present in sim dictionary?
    missing_keys = [k for k in software_required_keys if k not in sim or sim[k] is None]
    if missing_keys:
        msg = (
            f"Required '{sim['SOFTWARE'].upper()}' sim keys missing or "
            f"not defined in conf file: {', '.join(missing_keys)}"
        )
        raise YamlBadConfigError(msg)

    logger.debug(
        f"all {len(software_required_keys)} required '{sim['SOFTWARE'].upper()}' sim keys are present",
    )

    # STEP 4 - Check that all entry keys provided for each simulation are valid
    files_tbd = []

    #   loop config entries
    for key_sim, value_sim in sim.items():
        logger.debug(f"processing entry: sim['{key_sim}'] = {value_sim!s}")

        if key_sim.upper() == "SOFTWARE":
            continue

        # STEP 4.1.
        # Anne: check if key is in molecules_dict, molecule_numbers_dict or
        # molecule_ff_dict too
        if (
            (key_sim.upper() not in software_sim)
            and (key_sim.upper() not in solubles_set)
            and (key_sim.upper() not in lipids_set)
            and (key_sim.upper() not in molecule_ff_set)
        ):
            _es = (
                f"key_sim '{key_sim}' in {sim['SOFTWARE'].lower()}_dict' : {key_sim.upper() in software_sim}",
                f"key_sim '{key_sim}' in molecules_dict : {key_sim.upper() in solubles_set}",
                f"key_sim '{key_sim}' in lipids_dict : {key_sim.upper() in lipids_set}",
                f"key_sim '{key_sim}' in molecule_ff_dict : {key_sim.upper() in molecule_ff_set}",
            )
            logger.error(_es)
            msg = (
                f"'{key_sim}' not supported: Not found in "
                f"'{sim['SOFTWARE'].lower()}_dict', 'molecules_dict',"
                f" 'lipids_dict' and 'molecule_ff_dict'"
            )
            raise YamlBadConfigError(msg)
        if key_sim.upper() not in software_sim:  # hotfix for unkown yaml keys. TODO improve check 4.1?
            logger.warning(
                f"ignoring yaml entry '{key_sim}', not found in '{sim['SOFTWARE'].lower()}_dict'",
            )
            continue

        # STEP 4.2.
        # entries with files information to contain file names in arrays
        if "TYPE" in software_sim[key_sim]:
            if "file" in software_sim[key_sim]["TYPE"]:  # entry_type
                logger.debug(
                    f"-> found '{key_sim}:{software_sim[key_sim]}' of 'TYPE' file",
                )  # DEBUG

                if value_sim is None:
                    logger.debug(f"entry '{key_sim}' has NoneType value, skipping")
                # already a list -> ok
                elif isinstance(value_sim, list):
                    logger.debug(f"value_sim '{value_sim}' is already a list, skipping")
                    files_tbd.extend(value_sim)
                else:
                    value_sim_splitted = value_sim.split(";")

                    if len(value_sim_splitted) == 0:
                        msg = f"found no file to download for entry '{key_sim}:{software_sim[key_sim]}'"
                        raise YamlBadConfigError(msg)
                    # in case there are multiple files for one entry
                    if len(value_sim_splitted) > 1:
                        files_list = [x.strip() for x in value_sim.split(";")]
                        sim[key_sim] = files_list  # replace ; separated string with list
                    else:
                        sim[key_sim] = [
                            [f.strip()] for f in value_sim_splitted
                        ]  # IMPORTANT: Needs to be list of lists for now
                    files_tbd.extend(f[0] for f in sim[key_sim])

                # STEP 4.3.
                # Batuhan: In conf file only one psf/tpr/pdb file allowed each
                # (can coexist), multiple TRJ files are ok
                # TODO true for all sim software?
                # TODO add dict entry param "unique" instead?
                if key_sim.upper() in ["PSF", "TPR", "PDB"] and len(sim[key_sim]) > 1:
                    msg = f"only one '{key_sim}' entry file allowed, but got {len(sim[key_sim])}: {sim[key_sim]}"
                    raise YamlBadConfigError(msg)

        else:
            logger.warning(
                f"skipping key '{key_sim}': Not defined in software_sim library",
            )

    logger.info("found %d resources:\n%s", len(files_tbd), "\n".join(files_tbd))

    return sim, files_tbd
