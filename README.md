# Blood-Circulatory-System

University Group Project in partnership with the High Performance Computing Centre Stuttgart (HLRS).

Project Title: Modelling and Numerical Solution of the Blood Circulatory System. \
Creators: Berlinski Tomer, Costa Beatrice, Ehrlich Manuel and Moral Sánchez Elena. \
Matlab simulation (all the code in this repository) by Ehrlich Manuel and Moral Sánchez Elena.\
Supervised by Prof. Rainer Callies and Dr. Tobias Koeppl


The model is a compartmental model with a focus on the heart valves resulting in a system of ODEs. The system was solved numerically on Matlab using appropriately modified adaptive Range-Kutta methods. One of the challenges was to deal with non-differentiable points in the right hand side of the ODE.

The presentation and the poster attached were presented towards (but before) the end of our work. The report contains all our work. The Matlab code in the repository allows one to reproduce a simulation of a blood circulatory system model. 



------------------------------------------------------------------------------------------------------------------------------------------------------



In this project we have modelled and simulated the circulatory system using 14 compartments: [large arteries, arterioles, capillaries, venules, large veins, left atrium, left ventricle, large arteries, arterioles capillaries, venules, large veins, right atrium, right ventricle].

You can choose between four numerical methods. You can choose between 2 different implementations of the adaptive time step, with a tradeoff between speed and accuracy. We have non-differentiable points in the model. To address that you can choose whether to use a slightly modified model, or use a root-finding method.

METHOD_num has options @RKF45, @RK_Bogacki_Shampine, @RK_Cash_Karp, @RK_Dormand_Prince.\
METHOD_adapative has options @adaptive_deprecated (faster but less accurate) and @adaptive (slower but more accurate).\
METHOD_diff has options @F (with root_finding) and @F_modified (using modified RHS).\
