# MetaMolES
![logo](https://github.com/theicechol/metamoles/blob/master/figures/metamoles_logo.png)

UW-DIRECT Project on metabolite retrosynthetic analysis
- Aim to utilize data science and software engineering intuition to find, and predict, a plausible metabolic pathway for production of a given molecule with retrosynthetic analysis approach. 
- Current focus is to find a novel promiscuous substrate for enzymatic transformation.

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

![Outline_1](https://github.com/theicechol/metamoles/blob/master/figures/project_outline1.png)
![Outline_2](https://github.com/theicechol/metamoles/blob/master/figures/project_outline2.png)
![function_outline](https://github.com/theicechol/metamoles/blob/master/figures/function_flowchart.png)

funciton outline will be edited with correct function names and better outlook! 
