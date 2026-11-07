ResSim Studio - Windows release (MinGW build)

This folder holds only the programs:
  ResSimStudio.exe   desktop application (Model, Results, Geomodel, Studies)
  ResSimCLI.exe      command-line runner (--regress, --compare, MPI runs)
  ResGeo.exe         standalone geomodelling program
  ressim.dll         C API (Python package: ../python)
  USER_GUIDE.html    user guide (Help > User guide)
  lang/              UI translations

Everything else is kept once, in the parent folder (not copied here):
  ../examples, ../decks   example and test decks (the application lists them from there)
  ../python               Python package and examples
  ../research             research reports, data and papers
  ../docs, ../*.md        user guide source, ENGINES_AND_MATH, ROADMAP, INDUSTRY_GAP, PETREL_PARITY, UI review
  ../tools                helper scripts (OPM test download, nightly regression, i18n)
The test program (ressim_tests.exe) is built by build_mingw.sh into build/mingw and is not shipped.
