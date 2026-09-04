# Open a DDI assessment report template

Crate a new DDI report template as R markdown file in the working
directory, and open it. The new report contains examplinib data for
reference. Requires RStudio for the editor; otherwise the file is
written and the path is returned.

## Usage

``` r
ddi_report(name = "DDI-assessment.qmd", overwrite = FALSE)
```

## Arguments

- name:

  Destination file. Defaults to `DDI-assessment.qmd` in the working
  directory.

- overwrite:

  Overwrite `name` if it already exists.

## Value

The destination filename, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
ddi_report()
} # }
```
