# Runnable MicroBiotR tutorial

[Download the R Notebook source](MicroBiotR-tutorial.Rmd) and [download the executed 16S HTML example](MicroBiotR-tutorial.html). On GitHub, open the file and use **Download raw file**. Save the `.Rmd` to an analysis folder and open it in RStudio. The HTML example opens directly in a browser after downloading; GitHub's ordinary file view is not an HTML hosting service.

The notebook uses the Flatly theme, Tango code highlighting, floating three-level contents, numbered sections, paged data frames, and inline colors/typography matching the supplied DADA2 data-integrity notebook. Its CSS is included inside the `.Rmd`, so no external stylesheet is needed.

The notebook is a standalone tutorial, not a package vignette. Install MicroBiotR and the tutorial dependencies in the [root README](../README.md). Package installation does not copy this repository file or download datasets.

## Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `workflow` | `16s` | `cytometry` runs steps 1–13, `16s` runs 14, `rna` runs 15, `all` runs all workflows. |
| `data_dir` | `data` | Folder containing the matching Zenodo example files. Relative paths resolve from the notebook folder when knitting. |
| `output_dir` | `results` | Root for separate `cytometry`, `16s`, and `rna` output folders. Repeated runs may overwrite those outputs. |
| `download_data` | `false` | Download missing required example files from Zenodo. Existing files must pass the published MD5 checksums. |
| `execute` | `true` | Execute selected analysis chunks. Use `false` to render code and prose only. |
| `reuse_som` | `false` | For cytometry, load `SOM.csv`, `cluster_information.csv`, and mapped FCS files already generated in the selected output folder from the same inputs. |

## Execute or render

Edit the YAML header and choose **Run All** in RStudio to execute chunks. Notebook **Preview** shows saved chunk outputs; it does not rerun every analysis. To execute the selected workflow and produce a standalone HTML report, choose **Knit to HTML** or run this from the folder containing the `.Rmd`:

```r
rmarkdown::render("MicroBiotR-tutorial.Rmd", output_format = "html_document")
```

For another workflow, override the parameters:

```r
rmarkdown::render("MicroBiotR-tutorial.Rmd", output_format = "html_document",
                  params = list(workflow = "rna", data_dir = "data",
                                output_dir = "results-rna"))
```

`rmarkdown`, `knitr`, and Pandoc are required for rendering. RStudio bundles Pandoc; outside RStudio, make it available on PATH or through `RSTUDIO_PANDOC`. The notebook does not install or update packages. Complete cytometry SOM training can take substantial time and memory; cached outputs must correspond to the same data and parameter settings.

The distributed HTML contains actual 16S outputs and code for all workflows. The original README images remain historical illustrations. Computational checks do not establish biological validity or independent classifier performance.

## Maintainer validation

Keep analysis outputs outside Git. Render `html_document` for 16S and RNA, and render cytometry/all modes against validated cached SOM outputs. Test `execute = FALSE` and notebook preview separately. Use the [current test report](../docs/tutorial-local-test-report.md) for tested versions, scope, and reproducible commands. Changes to README analysis calls should be reflected in the notebook; its paths and chunk execution flags intentionally differ from the linear README.
