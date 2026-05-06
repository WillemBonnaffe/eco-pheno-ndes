# Eco-pheno NDEs

## Overview

This repository contains data and code to define and train a neural difference equation system on time series of mean body mass of fish and population size.
The aim is to quantify the effects of population density, mean body mass, winter and summer temperature, and harvesting on phenotype and population dynamics.
The approach identifies the degree of non-linearity in effects and interactions between variables that is supported by the time series data.
For more details see the associated publication ([ref]).

## Setting up

You need to install R. 
Computations were performed in R version 4.5.2 (2025-10-31).
Running the code only relies on base R packages.

## General structure of repository

All code and data can be found in the folder `src`.
All files that begin with the file indicator `f_`, e.g. `f_slp.r`, only contain function definitions supporting computations in the main scripts.
The main executable scripts are marked by the file indicator `m<script_number>_` where the script number indicates the order in which to scripts should be executed, e.g. `m1_prepare_raw_data.r` needs to be executed first.

## File descriptions

* `src/raw`
|_ * `src/raw`
