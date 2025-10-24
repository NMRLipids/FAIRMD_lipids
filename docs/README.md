To manually recreate documentation, run:

```bash
tox -e docs
```

**Note!** You should be already in the environment suitable for working with
fairmd-lipids project.

Then run:
```
make html
```
If you want to run apidoc which produces rst files from the project folder structure
under Scripts and produce html files run:

```
make all
```
This command combines running the python file: `source/run_apidoc` with the `make html`
command.

### Custom module rst files
If you want to include custom templates for individual modules, for instance AddData.py: 

Look at name of the RST file produced by `make all` within `source/auto_gen` folder,
with this example it's `Scripts.BuildDatabank.AddData.rst` then put your custom rst file
with the same name within this folder. It will not be overwritten. 

### Custom apidoc and sphinx templates used 

There are custom templates used by apidoc within the `source/_templates` folder. The
main change within `source/_templates/footer.html` template relative to the standard
template is adding links to repository and license. These links can be changed through
the python file `source/conf.py` with the dictionary called `html_context`.

The main change to the `source/_templates/package.rst.jinja` is removing differentiation
between namespace packages and normal python packages.

The original apidoc templates can be found
[here](https://github.com/sphinx-doc/sphinx/tree/master/sphinx/templates/apidoc). The
original footer template can be found
[here](https://github.com/readthedocs/sphinx_rtd_theme/blob/master/sphinx_rtd_theme/footer.html).

Deleting the custom templates within `source/_templates` would make Sphinx use default
templates. 
