# Code Availability — Draft text for manuscript

## Recommended version (Nature Communications)

Code for the tandem repeat expansion analyses is available at the GitHub
repository ASD-SZ-WGS-Japan:
https://github.com/KushimaItaru/ASD-SZ-WGS-Japan (commit `e0664d2`).
The repository contains the main entry-point wrapper scripts for the STRling
and ExpansionHunter Denovo (EHdn) pipelines, the underlying helper scripts,
configuration files, and documentation required to reproduce the primary
tandem repeat burden analyses. Wrapper scripts organize execution order and
job dependencies but do not alter the analytical logic.

## Shorter version (if space is limited)

Code for the tandem repeat expansion analyses is available at
https://github.com/KushimaItaru/ASD-SZ-WGS-Japan (commit `e0664d2`).
The repository includes the main entry-point wrapper scripts for the STRling
and EHdn pipelines, together with helper scripts, configuration files, and
documentation for reproducing the primary burden analyses.

---

*Note: Replace the commit hash above with the final release tag (e.g.,
`v1.0.0`) before manuscript submission.*

---

## GitHub About field (repository description)

Publication-ready tandem repeat expansion burden analysis pipeline for
GRIFIN-PD WGS using STRling and EHdn.

## GitHub README opening paragraph

This repository contains the publication-ready workflow for tandem repeat
expansion (TRE) burden analyses in the GRIFIN-PD whole-genome sequencing
study. It provides entry-point wrapper scripts for the STRling and EHdn
pipelines, together with the underlying helper scripts, execution order,
shared configuration, and documentation needed to reproduce the main burden
analyses.
