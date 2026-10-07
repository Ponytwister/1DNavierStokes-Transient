# Navier desktop test package

Extract the entire ZIP to a writable folder and open `bin/NavierGui.exe`.
Keep its DLLs, `qt.conf`, and the plugins folder together. Qt and a compiler do
not need to be installed on the destination machine. This build targets x64
Windows 10 (1809+) / Windows 11 with Qt 6.10.3 and MinGW GCC 13.1.0.

Choose an existing experiment database and an output directory from the **File**
menu, then Run.
Open a `.txt`, `.csv`, or `.tsv` report with NavierGui from File Explorer to
load it directly in the **Report** tab. The Report tab's **Open report...**
button accepts the same formats. CSV and TSV delimiters are detected from the
header; existing space-delimited reports remain supported.

Open the **Model controls** tab to edit settings. Navigation, closing, and Run
remain available while editing. Run asks you to accept changed controls; unchanged
controls run without a prompt. **Update preset** saves to the loaded setup file;
**Save preset...** prompts for a new file. Database defaults are unchanged.
**Cancel** discards edits not yet accepted for a run. Starting a run clears old results.
No experiment database is shipped or selected automatically. Test with a copy
of your database: loading can create solution records. Report export, profile
saving, and fitted-input saving are explicit actions. Cancellation is cooperative;
closing waits for workers, and failed saves keep the window open for retry.

The **Report** tab lets you select experiment/profile blocks, profile types and
metadata columns. After a run choose **Prepare report...** or **Use current run**
to generate directly from results, without writing a source file first. Existing
text files can still be opened. Direct generation offers separate analytical-zero
and true bound-bead rows plus every species' model-output profile in its declared
model units. It retains full double precision until you choose a number format.
Toggle profile values, headers, units/axes and blank rows to control the layout.
**D-A / integrals** selects the original stacked D-A layout or six separate,
individually selectable metadata columns. Split columns repeat the block's
values on every selected profile row; missing source values remain `-`.
Choose original numeric text, fixed decimals or scientific notation; formatting
rounds data values for presentation without recalculating results or axes.
The preview shows exactly the selected output. **Copy for Excel** copies that
output with tabs between cells; paste using Excel's **Paste Special > Text** to
keep your sheet's formatting. **Export selected report...** writes the same data
to a UTF-8 `.tsv` or `.txt` file for tab-delimited import. No data selection disables
copy/export. **Export legacy report** retains the original 13-row text format,
including its historical analytical-zero row labeled as bound beads.
Opening a file resets selections to All and keeps formatting choices. Generating
again retains matching selections; newly available types default to selected.
Report settings do not change the source data or database. Excel interprets
imported numbers using its own precision and locale settings.

## Dependency and redistribution record

This is a local development/test package, not a completed public release.
The repository has not selected an application distribution license. Before
redistribution, resolve the combined application license, provide corresponding
source/build instructions and all applicable notices, and validate on a clean
Windows machine. No code signing or installer is provided.

- Qt 6.10.3 Core, Gui, Widgets and plugins are shared libraries. The included
  `notices/qt-sbom` records the installed Qt Base kit, including third-party
  components; it is a superset of the deployed binaries, not a package manifest.
  Obtain the matching Qt source, license texts and third-party notices from
  https://download.qt.io/archive/qt/6.10/6.10.3/ before redistribution. The SBOM
  is not a substitute for those texts or source obligations.
- ALGLIB 4.02.0 is statically linked. Its source headers declare GPL v2 or later;
  copies of GPL v2 and v3 are in `notices/alglib`. A Qt LGPL deployment alone
  does not settle the combined application's distribution requirements.
- Eigen 3.4.0 is header-only; its supplied license files are in `notices/eigen`.
- SQLite is statically linked; its source declares public-domain dedication.
  See https://sqlite.org/copyright.html.
- GCC 13.1.0 / MinGW-w64 runtime libraries are dynamically linked. Their supplied
  GCC runtime exception, GPL and MinGW/winpthreads notices are in
  `notices/compiler`.

The package intentionally contains no production databases, generated model
results, tests, compiler executables, or Qt development headers.
