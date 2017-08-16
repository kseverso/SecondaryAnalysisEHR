# Secondary Analysis of Electronic Healthcare Records

This repository contains code to perform secondary analysis of electronic healthcare records obtained from [MIMIC](https://mimic.physionet.org/), an open, de-identified database of healthcare records. Specifically, this project was interested in the impact of vasopressin on serum lactate levels in ICU patients with sepsis. Eligible patients were over 18 years old, were admitted to either the medical or surgical ICU, and had serum lactate measurements taken 6 hours prior to or up to 6 hours after ICU admission with a subsequent measurement taken during hours 21-27.

The analysis is done in two parts. The jupyter notebooks query the database and apply most of the inclusion criteria. The R script performs the statistical analysis using [DesignMatch](https://cran.r-project.org/web/packages/designmatch/designmatch.pdf). The data files are not available and the code is provided for reference.
