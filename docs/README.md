# Analysis and maintenance documentation

This folder explains the workflows, statistical assumptions, input validation and current implementation limitations.

* [Micro.R](Micro.md): every public workflow function, its operations, rationale, results and examples.
* [dependent.R](dependent.md): the six low-level helpers and their data contracts.
* [Gating](Gating.md): instrument-specific standalone preprocessing and all six workflow helpers.
* [Shared validation](validation.md): sample-ID matching, matrices and independent test evaluation.
* [Statistical methods](statistical-methods.md): assumptions, what reported tests mean and primary references.
* [Examples and validation](examples.md): how to run synthetic examples and check documentation.
* [Validation results](validation-results.md): checks performed and remaining verification limits.
* [Implementation notes](implementation-notes.md): remaining limitations and requested issue fixes.

The executable scripts retain their original paths. Markdown guides are repository documentation; installed R help is in [`man/`](../man/README.md). No real FCS or sample metadata are bundled.
