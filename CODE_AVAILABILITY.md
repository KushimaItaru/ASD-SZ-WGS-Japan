# Code Availability — Draft text for manuscript

## Recommended version (Nature Communications)

Code used for the tandem repeat expansion analyses is available at
https://github.com/KushimaItaru/ASD-SZ-WGS-Japan. The repository includes the
main entry-point wrapper scripts for reproducing the STRling and EHdn analyses,
together with the underlying helper scripts, configuration files, and
documentation described in the README. For STRling, the primary entry-point
scripts are `strling/01_strling_build_panel.sh`,
`strling/02_strling_casecontrol_call_and_outliers.sh`,
`strling/03_strling_casecontrol_burden.sh`, and
`crosscaller/04_tre_crosscaller_compare.sh`. For EHdn, the primary entry-point
scripts are `ehdn/01_ehdn_setup_and_sample_lists.sh`,
`ehdn/02_ehdn_depth.sh`, `ehdn/03_ehdn_profile_and_merge.sh`, and
`ehdn/04_ehdn_casecontrol_burden.sh`. These wrapper scripts organize execution
order and job dependencies but do not alter the analytical logic.

## Shorter version (if space is limited)

Code for the tandem repeat expansion analyses is available at
https://github.com/KushimaItaru/ASD-SZ-WGS-Japan. The repository contains
the main entry-point wrapper scripts for the STRling and EHdn pipelines,
the underlying helper scripts, configuration files, and a README describing
execution order and dependencies. The wrapper scripts organize workflow
execution and job dependencies without changing the analytical logic.

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
