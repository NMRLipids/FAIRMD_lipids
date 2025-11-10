# General rules, code style, and formating

Regarding versions:
- We support python3.10+

Regarding conventions:
- We use `ruff` as a linter and code-style checker.
- We use `variable_name`, `ClassName`, and `method_name`. For global variables, we use `UPPER_CASE_NAMES`.
- Sometimes, this convention can be broken, but a rule should be written in the end to allow `ruff` to ignore the breach.
- We recommend using extended `ruff` checking by setting `developer/ruff-dev.toml` as a default ruff condig. This set of rules is agreed to be our guide in the code style. 

Please install pre-commit hook:
```
pre-commit install
```
This will make sure that no code can be committed to the repository until it satisfies the `ruff` check.

We are using automatic testing with `pytest` please consult (test-README)[Scripts/tests/README.md]. We use both unit and regression testing, however the coverage of code with unit-tests remains relatively low.

# Repository rules

1. Please create an issue, and then pull-request, not vice versa.
2. Please keep main branch in your fork updated so that it can be rebased (we prefer rebasing over merging).
3. Once you address an issue, please refer it in the commit message by using a phrase such a 'Partially fixes #000'.
4. You must make sure that the tests don't break before you make a PR.
5. We don't oblige contributors to squash commits. We appreciate if the commit has a considerable size.
6. We require at least one review from organisation member to accept a pull-request.
7. Check representation of yourself in `AUTHORS.md` file.
8. You are always welcome to participate in repository discussions and NMRlipids community events to develop code together.

# PATH handling in our code

As the code operates with an over-filesystem database, universal path handling is crucial. To construct paths, we use global variables of the type NMLDB\_**XXXX**\_PATH, where **XXXX** could be:

- DATA points by default to `./BilayerData` folder. Can be used to point to somewhere inside the Database but not for molecules, simulations, and experiments. Currently is used for Rankings as well.
- MOL points by default to `{FMDL_DATA_PATH}/Molecules` folder, from where molecule lists are initialized. Is used to get access to something in molecule-folders.
- SIMU points by default to `{FMDL_DATA_PATH}/Simulations` folder. Is used to get access to a certain simulation.
- EXP points by default to `{FMDL_DATA_PATH}/Experiments` folder. Is used to get access to a certain experiment.

We currently construct paths using `os.path.join(a,b)`.

#  Getting started

To help with developing start by installing the development dependencies. Our continuous
integration pipeline is based on [Tox](https://tox.readthedocs.io/en/latest/). So you
need to install `tox` first

```bash
    pip install tox
    # or
    conda install -c conda-forge tox
```

Then go to the [develop project](https://github.com/NMRLipids/FAIRMD_lipids/) page, hit the
``Fork`` button and clone your forked branch to your machine.

```bash
  git clone git@github.com:your-user-name/fairmd-lipids.git
```

Now you have a local version on your machine which you can install by

```bash
  cd fairmd-lipids
  pip install -e .
```

This install the package in development mode, making it importable globally and allowing
you to edit the code and directly use the updated version. To see a list of all
supported tox environments please use

```bash
  tox list
```

# Running the tests

The testsuite is implemented using the [pytest](https://docs.pytest.org/en/stable/)
framework and should be set-up and run in an isolated virtual environment with
[tox](https://tox.readthedocs.io/en/latest/). All tests can be run with

```bash
  tox                  # all tests
```

If you wish to test only specific functionalities, for example:

```bash
  tox -e lint          # code style
  tox -e tests-min     # unit tests of the main library
  tox -e tests-all     # regression tests
```

You can also use `tox -e format` to use tox to do actual formatting instead of just
testing it. Also, you may want to setup your editor to automatically apply the
[ruff](https://ruff.rs/docs/) code formatter when saving your files, there are plugins
to do this with all major editors.

# Contributing to the documentation

The documentation is written in reStructuredText (rst) and uses
[sphinx](https://www.sphinx-doc.org) documentation generator. In order to modify the
documentation, first create a local version on your machine as described above. Then,
build the documentation with

```bash
    tox -e docs
```

You can then visualize the local documentation with your favorite browser using the
following command (or open the :file:`docs/build/html/index.html` file manually).

```bash

    # on linux, depending on what package you have installed:
    xdg-open docs/build/html/index.html
    firefox docs/build/html/index.html

    # on macOS:
    open docs/build/html/index.html
```

# Data handling

FAIRMD Lipids separates codespace from
[the Database](https://github.com/NMRLipids/BilayerData) since June 2025 (v.1.1.0). 
Data contribution rules are moved there accordingly.

# Contributing data

If you would like to contribute simulations or experiments, please consult instructions for [Data Contribution](https://nmrlipids.github.io/FAIRMD_lipids/stable/dbcontribute.html).
