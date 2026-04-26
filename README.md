# Raleigh Crime Spatial Analysis

## Overview
This project analyzes Raleigh, North Carolina crime data using spatial point process methods. The two crime categories studied are Larceny and Simple Assault.

## Tools Used
- R
- sf
- splancs
- MASS
- maps

## Methods
- Complete spatial randomness (CSR) testing
- L-function analysis
- Monte Carlo simulation envelopes
- Log relative risk surface estimation
- Monte Carlo relabeling for statistical significance

## Analysis
The project treats Larceny and Simple Assault incidents as spatial point patterns within the Raleigh township. Each pattern is compared against complete spatial randomness using L-functions and simulation envelopes.

The two crime types are also compared using a log relative risk surface to identify areas where one crime type is relatively more concentrated than the other.

## Key Findings
- Larceny showed strong clustering across most distances.
- Simple Assault showed clustering at shorter distances, randomness at intermediate distances, and inhibition at larger distances.
- The log relative risk surface was close to zero across most of Raleigh, suggesting broadly similar spatial patterns.
- Monte Carlo relabeling identified localized areas where Simple Assault was significantly more concentrated, with a smaller region where Larceny was more concentrated.

## Goal
Demonstrate how spatial statistics and point process methods can be used to analyze geographic crime patterns and identify statistically meaningful differences between crime types.
