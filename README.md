# MetaMolES
![logo](https://github.com/theicechol/metamoles/blob/master/figures/metamoles_logo.png)

See our poster for the quarter summary [here](https://github.com/theicechol/metamoles/blob/master/docs/MetaMolES-Poster.pdf)

### Current Work

**UW-DIRECT Project on metabolite retrosynthetic analysis**
- Aim to utilize data science and software engineering intuition to find, and predict, a plausible metabolic pathway for production of a given molecule with retrosynthetic analysis approach. 
- Current focus is to find a novel promiscuous substrate for enzymatic transformation.

### Future Work
•	Add compounds with canonical SMILES string but no isomeric SMILES string
•	Test inclusion of additional features, such as full chemical fingerprints and enzyme descriptors
•	Explore alternatives models such as SVMs, neural networks, decision trees/random forests, and ensemble methods
•	Extend approach to include non-promiscuous enzymes
•	Include simple chemical transformation for biocatalysis application


### For thorough explanation and project background --- Check out our [Wiki](https://github.com/theicechol/metamoles/wiki)

# Team Members

#### Ellie (elliej3), Ice (theicechol), Phil (philipjleung), Stephen (blasks), and Yeon Mi (ymhwang414)

## Use Cases

**Basic**
- **Search for enzyme to create input molecule:** This use cases starts when a user inputs the PubChem ID of a molecule into the Metamoles package. This use case ends when a prediction about an enzyme that could produce the input molecule is made.
Components:
  - User input
    -ID
  - Enzyme data collection
    - KEGG promiscuous reactions info
    - SMILES of products, etc
  - Compound data collection
    - SMILES
  - Data generation for model
    - Comparison of products for each enzyme
  - Model generation
    - Train and test
  - Model output
    - Prediction
    - Quality of prediction
- **Retrieve information about model used for prediction:** This case starts when a user inputs the PubChem ID of a molecule. This use case ends when a prediction is made.
Components:
  - User input
    - ID
  - Enzyme data collection
    - KEGG promiscuous reactions info
    - SMILES of products, etc
  - Compound data collection
    - SMILES
    - Data generation for model
    - Comparison of products for each enzyme
  - Model generation
    - Train and test
  - Model output
    - Prediction
    - Quality of prediction

- **Train the regression/ML prior to deployment:** This use case begins when a developer inputs the PubChem ID of a molecule. This use case ends when the model has accessed all the training data.
Components:
  - Enzyme data collection
  - Compound data collection
  - Data generation for model
    - Comparisons
  - Regression/ML selection
  - Train the model with generated data
- **Test the regression/ML:** This use case starts when a developer uses known information to test the reliability of the model. This case ends when the model meets a level of accuracy that is appropriate for deployment.
Components:
  - Trained model
  - Withheld data

**Advanced**
- Real time search (calculate when called) versus recent update search (some data on hand)?
- Explore enzyme promiscuity (look at what’s happening in the background)
- Report experimental outcomes (positive and negative if people use our prediction)
- Offer different ML models or prediction methods (edited) 
- Included non-promismuous enzyme in the reaction search (or some non-enzymatic pathway)

## A core workflow for MetaMolES
![CoreWorkflow](https://github.com/theicechol/metamoles/blob/master/figures/CoreWorkflow.png)

## How we calculate distance of desired compounds with enzyme substrate scope
![Distance](https://github.com/theicechol/metamoles/blob/master/figures/project_outline2.png)

## How we curate the data
![Functions](https://github.com/theicechol/metamoles/blob/master/figures/DataCuration_Query_Workflow.png)
