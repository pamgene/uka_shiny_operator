# upstream_kinase_analysis_shiny_operator

##### Description

The `upstream_kinase_analysis_shiny_operator` converts a list of peptides into a list of kinases. 

##### Usage

Input projection|.
---|---
`y-axis`        | normalized data (log / VSN / Combat-corrected), single y value / cell
`row`           | peptides
`column`| observations (samples)
`color`| grouping factor (no color needed if input is single observation / foldchange.)


Input parameters|.
---|---
`Kinase_family`      | Type of kinase family, either PTK or STK, default is PTK
`Lock_kinase_family` | Either Yes or No, default is Yes



##### Details

The input is a cross-tab view with a single value per cell. A data color is needed to group the observations if there are two groups (test / control). If input is a single observation (usually foldchange), color is not needed.

##### Versions and UKA database:

Version 0.3.0 and before: UKA db 2023 (contains Array Layouts 86312, 86402, 86412, 87102)
- version 0.3.3: scoreplot Kinase name colors fixed
Version 0.4.0: UKA db 2024 (contains Array Layouts 86312, 86402, 86412, 87102, 87202)
