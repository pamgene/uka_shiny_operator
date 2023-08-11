# upstream_kinase_analysis_shiny_operator

##### Description

The `upstream_kinase_analysis_shiny_operator` converts a list of peptides into a list of kinases. 

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, single y value per cell
`row`           | peptides
`column`| observations
`color`|grouping factor for the analysis (no color needed if input is single observation / foldchange.)


Input parameters|.
---|---
`Kinase_family`      | Type of kinase family, either PTK or STK, default is PTK
`Lock_kinase_family` | Either Yes or No, default is Yes



##### Details

The input is a cross-tab view with a single value per cell. A data color is needed to group the observations if there are two groups. If input is a single observation (usually foldchange), color is not needed.

