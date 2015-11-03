# NeuroGLM

This will be a MATLAB toolbox for the analysis of neural and behavioural data under a decision making task. The toolbox will incorporate models to quantify to what extent recorded neural spiking activity relates to external and internal variables of decision.

Classical approaches rely on Signal Detection Theory to quantify latent features of perception (sensitivity, bias). A new GLM approach has been developed which offers greater flexibility for modeling the effect of internal and external variables of interest on a task.

## How to use
### Loading GLM object
The GLM object is created using a data block provided by *expRef*
```
>> g = GLM(expRef)
```

The GLM object contains the following fields:
```
>> g = GLM('2015-05-28_1_Laveran');
>> g
g =
  GLM with properties:

             expRef: '2015-05-28_1_Laveran'
        modelString: []
    parameterLabels: []
      parameterFits: []
    parameterBounds: []
     parameterStart: []
                 ZL: []
                 ZR: []
               data: [1x1 struct]
              p_hat: []
```
`g.data` is a struct containing the relevant behavioural data for that block. 

### Assigning a model
Models are assigned using the setModel method and a *modelString* argument
```
>> g = g.setModel(modelString)
```
Note that since MATLAB classes are **value classes**, functions applied to instances of the class must be reassigned back to the class. i.e. the syntax is `g = g.fcn` not  `g.fcn`.

The modelStrings currently included are:

* Offset
* ifC
* fullContrasts
* CL+CR
* C^N
* C^NL^NR
* C50
* Supersaturation

As well as subset models of these (append '-subset' to modelString).

```
>> g = g.setModel('C^N-subset');
>> g
g = 
  GLM with properties:

             expRef: '2015-05-28_1_Laveran'
        modelString: 'C^N-subset'
    parameterLabels: {'Offset_L'  'ScaleL_L'  'Offset_R'  'ScaleR_R'  'N'}
      parameterFits: []
    parameterBounds: [2x5 double]
     parameterStart: [0 0 0 0 0]
                 ZL: @(P,CL,CR)(P(1)+P(2).*CL.^P(5))
                 ZR: @(P,CL,CR)(P(3)+P(4).*CR.^P(5))
               data: [1x1 struct]
              p_hat: []
```

The file **GLM.m** outlines the configuration for each model. Your own model can be specified by defining the following properties:

* modelString
* Parameter labels
* Estimation bounds for the parameters
* ZL function (3 input arguments: parameter vector, contrast left scalar/vector, contrast right scalar/vector) 
* ZR function (3 input arguments: same as above)
* (optional) starting parameter values for optimisation

### Model fitting
The model can then be fit with the following
```
>> g = g.fit;
>> g
g = 
  GLM with properties:

             expRef: '2015-05-28_1_Laveran'
        modelString: 'C^N-subset'
    parameterLabels: {'Offset_L'  'ScaleL_L'  'Offset_R'  'ScaleR_R'  'N'}
      parameterFits: [-0.2790 5.8255 -1.3476 6.1142 0.5169]
    parameterBounds: [2x5 double]
     parameterStart: [0 0 0 0 0]
                 ZL: @(P,CL,CR)(P(1)+P(2).*CL.^P(5))
                 ZR: @(P,CL,CR)(P(3)+P(4).*CR.^P(5))
               data: [1x1 struct]
              p_hat: []
```

Leave-1-out crossvalidation can be called with `g.fitCV` and the p_hat field will be populated with predicted probabilities for held-out observation.

### Plotting
GLM object data can be plotted standalone with `g.plotData`. 

After running a non-crossvalidated fitting, the GLM fit itself can be plotted using `g.plotFit`

The parameters of a fit can be plotted using `g.plotParams`.

### Special case: optogenetic inactivation sessions 
The subclass **laserGLM** defines extra behaviour for experiments involving optogenetic inactivation. Sessions involving laser inactivation can be loaded into the GLM framework using `g = laserGLM(expRef)` and the same fitting/plotting functionality applies. Currently the code only removes laser trials but this subclass can be expanded to incorporate laser vs nolaser comparisons. 

### Special case: model simulations
The model fitting process can be validated by simulating from a specified model type with given parameter values. This is facilitated by the subclass **simulateGLM**. Simulations are based on contrast stimuli found in a real experimental session, and therefore simulation GLM objects are created using the following syntax
```
g = simulateGLM(expRef,trueModelString,trueModelParams)
```
Where expRef is the real experimental session from which contrast stimuli are taken, trueModelString & trueModelParams are model configurations used for generating the simulated data. The simulateGLM object can be fitted and plotted as usual.

## How to contribute
The class **GLM** is a superclass containing general methods for model fitting and plotting. Extra functionality can be added by defining a new subclass using the code framework:
```
classdef subclass < GLM
    % Subclass information
    
    properties
    end
    
    methods
        function obj = laserGLM(expRef)
            obj@GLM(expRef);
            % EXTRA CODE HERE FOR THE CONSTRUCTOR
        end
        
        function obj = fit(obj)
            % EXAMPLE FUNCTION OVERLOADING
        end    
    end
end
```
