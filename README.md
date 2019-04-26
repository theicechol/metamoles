# MetaMolES
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://travis-ci.com/metamoles/metamoles.svg?branch=master)](https://travis-ci.com/metamoles/metamoles)
[![Coverage Status](https://coveralls.io/repos/github/metamoles/metamoles/badge.svg?branch=master)](https://coveralls.io/github/metamoles/metamoles?branch=master)

![logo](https://github.com/theicechol/metamoles/blob/master/figures/metamoles_logo.png)

See our poster for the quarter summary [here](https://github.com/theicechol/metamoles/blob/master/docs/MetaMolES-Poster.pdf)

#### Current Team Members: Ellie (elliej3), Ice (theicechol), Phil (philipjleung)

##### Previous Team Members: Stephen (blasks), and Yeon Mi (ymhwang414)

### Current Work

**UW-DIRECT Project on metabolite retrosynthetic analysis**
- Aim to utilize data science and software engineering intuition to find, and predict, a plausible metabolic pathway for production of a given molecule with retrosynthetic analysis approach. 
- Current focus is to find a novel promiscuous substrate for enzymatic transformation.

### Future Work
-	Add compounds with canonical SMILES string but no isomeric SMILES string
- Test inclusion of additional features, such as full chemical fingerprints and enzyme descriptors
-	Explore alternatives models such as SVMs, neural networks, decision trees/random forests, and ensemble methods
- Extend approach to include non-promiscuous enzymes
- Include simple chemical transformation for biocatalysis application

### For thorough explanation, use cases, and project background --- Check out our [Wiki](https://github.com/theicechol/metamoles/wiki)

### Dependencies
- BioPython
- PubChemPy
- RDKit
- pandas
- sklearn
- scipy
- numpy
__________

## Graphical Explanation of our Goal
![Expected_Output](https://github.com/theicechol/metamoles/blob/master/figures/Expected_RetroSynthesis_Output.png)

## A core workflow for MetaMolES
![CoreWorkflow](https://github.com/theicechol/metamoles/blob/master/figures/CoreWorkflow.png)

## How we calculate distance of desired compounds with enzyme substrate scope
![Distance](https://github.com/theicechol/metamoles/blob/master/figures/project_outline2.png)

## How we curate the data
![Functions](https://github.com/theicechol/metamoles/blob/master/figures/DataCuration_Query_Workflow.png)
