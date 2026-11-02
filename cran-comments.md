## Submission

This is the first CRAN submission of this package (version 3.0.0). The
version number continues the numbering of the public GitHub releases.

## R CMD check results

0 errors | 0 warnings | 1 note

The note comes from the "CRAN incoming feasibility" check:

* This is a new submission.
* The DESCRIPTION cites Krishnamoorthy and Mathew (2009),
  <doi:10.1002/9780470473900>. The DOI resolves (doi.org and Crossref), but
  the publisher's website (Wiley) answers automated requests with HTTP 403,
  so the URL check may report it as "Forbidden".
* If the spell check flags words in the DESCRIPTION, they are the author
  names "Krishnamoorthy" and "Mathew" and the acronym "MCMC" (Markov chain
  Monte Carlo, spelled out in the text).

## Test environments

* local: Windows 11 x64, R 4.2.3 (`R CMD check --as-cran`: 0 errors,
  0 warnings; notes only for the new submission and the local environment)
* GitHub Actions (`R-CMD-check` workflow from r-lib/actions, all `Status: OK`):
  * ubuntu-latest: R-devel (2026-10-21 r90732), R 4.6.1 (release),
    R 4.5.3 (oldrel-1)
  * macos-latest (aarch64): R 4.6.1
  * windows-latest: R 4.6.1

## Notes for the reviewers

* Sanofi is listed as copyright holder and funder (`cph`, `fnd`) because the
  code was developed by Sanofi employees. The maintainer address is the
  maintainer's own mailbox and is monitored.
* The package is licensed under Apache-2.0; the NOTICE file is shipped in
  `inst/NOTICE` as required by section 4(d) of the license.
* `start_app_view()` launches an interactive Shiny application. Its example
  is wrapped in `if (interactive())`; all other examples run in well under
  one second on small synthetic data. The tests run without network access
  and write no files.
* The application code only reads files that the user uploads through the
  interface (up to 30 MB). It writes no files and makes no network requests
  itself (the Bootswatch theme supplied by 'bslib' may load a web font in
  the user's browser). The `shiny.maxRequestSize` option is set while the
  application runs and restored afterwards.
* The package contains one small synthetic data set (`data_InputFileTester2`,
  19 rows) for examples and tests. It contains no real data.
* The `URL` field lists the user guide (GitHub Pages, built from the `docs`
  folder of the repository) and the development repository. The guide uses
  simulated data only.
* The two-sided tolerance interval (Krishnamoorthy and Mathew, 2009) was
  corrected relative to the earlier GitHub versions: draws are counted when
  their own interval lies inside the candidate interval, the reported
  achieved credibility is the real achieved percentage, and the two-sided
  plots show the upper quantile on the x-axis and the lower on the y-axis.
  The tests cover this; see `NEWS.md`.
* There are no reverse dependencies (new package).
