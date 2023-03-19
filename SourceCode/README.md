# Source code of the implemented algorithm and closed-loop simulation test for artificial pancreas(AP):


Since the Matlab and C++ versions of the MPC algorithm are implemented using the CasADi framework, the CasADi file needs to be called to run the code. CasADi files can be downloaded directly from the official website. 

The input path of the "addpath" function in the Matlab version of the code needs to be changed to the path where the CasADi file is located.
	***/Matlab_version/LMPC/LMPC_single/CasADiInit_LMPC_single.m***
	***/Matlab_version/LMPC/LMPC_single/CasADiCompute_LMPC_single.m***
	***/Matlab_version/LMPC/LMPC_multiple/CasADiInit_LMPC_multiple.m***
	***/Matlab_version/LMPC/LMPC_multiple/CasADiCompute_LMPC_multiple.m***
	***/Matlab_version/NMPC/NMPC_single/CasADiInit_NMPC_single.m***
	***/Matlab_version/NMPC/NMPC_multiple/CasADiInit_NMPC_multiple.m***	


## The structure of the folder and the information of the main code file
### Matlab version of Linear model predictive control(LMPC) for single-hormone and closed-loop simulation test
Main files:
- ***DriverAP_LMPC_single.m:*** 		Driver file(main function).
- ***ParamEstimate_dataset.m:*** 		Generate data set used for parameter estimation.
- ***ParamEstimateMLE.m:*** 			The MLE method is used for parameter estimation.
- ***CasADiInit_LMPC_single.m:*** 		Initialization of OCP solver using CasADi.
- ***ClosedLoopSimulation_LMPC_single.m:*** 	Closed-loop simulation.
- ***KF_single.m:*** 				The discrete dynamic kalman filters.
- ***CasADiCompute_LMPC_single.m:*** 		Obtain the solution by OCP solver using CasADi.
- ***EulerMaruyamaSDE_Hovorka_single.m:*** 	Use the Euler-Maruyama method to simulate extended Hovorka model for single-hormone
- ***HovorkaModel_single.m:***		Extended Hovorka model for single-hormone used in closed-loop simulation.
  - ***HovarkaOutput.m***
  - ***HovarkaSensor.m***
  - ***HovorkaModelSteadyState_single.m:***	Extended Hovorka model for single-hormone used in steady state calculation.
  - ***HovorkaModely2x_single.m***
- ***MVPStateJacobian_single.m:***		Extended MVP model for single-hormone used in parameter estimation.
  - ***MVPStateJacobian_single_noDiffusion.m***
  - ***MVPOutputJacobian_single.m***		
  - ***MVPModelSteadyState_single.m:***	Extended MVP model for single-hormone used in steady state calculation.
  - ***MVPModely2x_single.m***		
---	
### Matlab version of Linear model predictive control(LMPC) for multiple-hormone and closed-loop simulation test
Main files:
- ***DriverAP_LMPC_multiple.m:***   		Driver file(main function).
- ***ParamEstimate_dataset.m:***		Generate data set used for parameter estimation.
- ***ParamEstimateMLE.m:*** 			The MLE method is used for parameter estimation.
- ***CasADiInit_LMPC_multiple.m:*** 		Initialization of OCP solver using CasADi.
- ***ClosedLoopSimulation_LMPC_multiple.m:*** 	Closed-loop simulation.
- ***KF_multiple.m:*** 			The discrete dynamic kalman filters.
- ***CasADiCompute_LMPC_multiple.m:*** 	Obtain the solution by OCP solver using CasADi.
-  ***EulerMaruyamaSDE_Hovorka_multiple.m:*** 	Use the Euler-Maruyama method to simulate extended Hovorka model for multiple-hormone.
- ***HovorkaModel_multiple.m:***		Extended Hovorka model for single-hormone used in closed-loop simulation.
  - ***HovarkaOutput.m***
  - ***HovarkaSensor.m***
  - ***HovorkaModelSteadyState_multiple.m:***	Extended Hovorka model for multiple-hormone used in steady state calculation.
  - ***HovorkaModely2x_multiple.m***
- ***MVPStateJacobian_single.m:***		Extended MVP model for single-hormone used in parameter estimation.
  - ***MVPStateJacobian_single_noDiffusion.m***
  - ***MVPOutputJacobian_single.m***		
  - ***MVPStateJacobian_multiple.m:***		Extended MVP model for multiple-hormone used in CDEKF.
  - ***MVPStateJacobian_single_noDiffusion.m***
  - ***MVPOutputJacobian_multiple.m***	
  - ***MVPModelSteadyState_multiple.m:***	Extended MVP model for multiple-hormone used in steady state calculation.
  -   ***MVPModely2x_multiple.m***	
  ---
### Matlab version of Nonlinear model predictive control(NMPC) for single-hormone and closed-loop simulation test
Main files:
-  ***DriverAP_NMPC_single.m:*** 		Driver file(main function).
-  ***ParamEstimate_dataset.m:*** 		Generate data set used for parameter estimation.
-  ***ParamEstimateMLE.m:*** 			The MLE method is used for parameter estimation.
-  ***CasADiInit_NMPC_single.m:*** 		Initialization of OCP solver using CasADi.
-  ***ClosedLoopSimulation_NMPC_single.m:*** 	Closed-loop simulation.
-  ***CDEKF_single.m:*** 			The continuous-discrete extended Kalman filter(CDEKF).
-  ***RK4_EKF_predict_single.m:*** 		Use the explicit Runge–Kutta 4(ERK4) method to predict the states and its covariance.
-  ***EulerMaruyamaSDE_Hovorka_single.m:*** 	Use the Euler-Maruyama method to simulate extended Hovorka model for single-hormone.
- ***HovorkaModel_single.m:***		Extended Hovorka model for single-hormone used in closed-loop simulation.
  - ***HovarkaOutput.m***
  - ***HovarkaSensor.m***
  - ***HovorkaModelSteadyState_single.m:***	Extended Hovorka model for single-hormone used in steady state calculation.
  - ***HovorkaModely2x_single.m***
- ***MVPStateJacobian_single.m:***	Extended MVP model for single-hormone used in parameter estimation.
  - ***MVPStateJacobian_single_noDiffusion.m***
  - ***MVPOutputJacobian_single.m***		
  - ***MVPModelSteadyState_single.m:***	Extended MVP model for single-hormone used in steady state calculation.
  - ***MVPModely2x_single.m***
---
### Matlab version of Nonlinear model predictive control(NMPC) for multiple-hormone and closed-loop simulation test
Main files:
-  ***DriverAP_NMPC_multiple.m:*** 		Driver file(main function).
-  ***ParamEstimate_dataset.m:*** 		Generate data set used for parameter estimation.
-  ***ParamEstimateMLE.m:*** 			The MLE method is used for parameter estimation.
-  ***CasADiInit_NMPC_multiple.m:*** 		Initialization of OCP solver using CasADi.
-  ***ClosedLoopSimulation_NMPC_multiple.m:*** 	Closed-loop simulation.
-  ***CDEKF_multiple.m:*** 			The continuous-discrete extended Kalman filter(CDEKF).
-  ***RK4_EKF_predict_multiple.m:*** 		Use the explicit Runge–Kutta 4(ERK4) method to predict the states and its covariance.
-  ***EulerMaruyamaSDE_Hovorka_multiple.m:*** 	Use the Euler-Maruyama method to simulate extended Hovorka model for multiple-hormone.
- ***HovorkaModel_multiple.m:***		Extended Hovorka model for single-hormone used in closed-loop simulation.
  - ***HovarkaOutput.m***
  - ***HovarkaSensor.m***
  - ***HovorkaModelSteadyState_multiple.m:***	Extended Hovorka model for multiple-hormone used in steady state calculation.
  - ***HovorkaModely2x_multiple.m***
- ***MVPStateJacobian_single.m	:***	Extended MVP model for single-hormone used in parameter estimation.
  - ***MVPStateJacobian_single_noDiffusion.m***
  - ***MVPOutputJacobian_single.m***		
  - ***MVPStateJacobian_multiple.m:***		Extended MVP model for multiple-hormone used in CDEKF.
  - ***MVPStateJacobian_single_noDiffusion.m***
  - ***MVPOutputJacobian_multiple.m***		
  - ***MVPModelSteadyState_multiple.m:***	Extended MVP model for multiple-hormone used in steady state calculation.
  - ***MVPModely2x_multiple.m***
---
### C++ version of Nonlinear model predictive control(NMPC) for single-hormone and closed-loop simulation test
Main files:
- ***main.cpp:*** 				Driver file to run Closed-loop simulation.
- ***CasADiSolver.h & CasADiSolver.cpp:*** 	Construct OCP solver by using CasADi.
- ***explicitRK4.h & explicitRK4.cpp:***		Use the explicit Runge–Kutta 4(ERK4) method to predict the states and its covariance.
- ***functionUsed.h & functionUsed.cpp:***	Function used in closed-loop simulation, such as the process noise and sensor noise.
- ***HovorkaModel.h & HovorkaModel.cpp:***	Extended Hovorka Model for single-hormone.
- ***MVPmodel.h & MVPmodel.cpp:***  		Extended MVP Model for single-hormone.
---
### C++ version of Nonlinear model predictive control(NMPC) for multiple-hormone and closed-loop simulation test 
Main files:
- ***main.cpp:*** 				Driver file to run Closed-loop simulation.
- ***CasADiSolver.h & CasADiSolver.cpp:*** 	Construct OCP solver by using CasADi.
- ***explicitRK4.h & explicitRK4.cpp:***		Use the explicit Runge–Kutta 4(ERK4) method to predict the states and its covariance.
- ***functionUsed.h & functionUsed.cpp:***	Function used in closed-loop simulation, such as the process noise and sensor noise.
- ***HovorkaModel.h & HovorkaModel.cpp:*** 	Extended Hovorka Model for multiple-hormone.
- ***MVPmodel.h & MVPmodel.cpp:***  		Extended MVP Model for multiple-hormone.
