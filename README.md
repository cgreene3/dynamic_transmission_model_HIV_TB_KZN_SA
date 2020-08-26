# Impact of Isoniazid Preventive Therapy for HIV-Infected Adults in South Africa: A Compartmental Epidemiological Model

This  compartmental  model  evaluates  the  impacts  of  different  strategies  to  provide  isoniazid preventive therapy (IPT) to prevent tuberculosis (TB) among people living with HIV (PLHIV). The model was developed in the context of rural KwaZulu-Natal,  South Africa for the Delivery Optimization for Antiretroviral Therapy study (DO ART). This model was adapted from Dowdy, David W., et al. (2014) who developed a similar model for application in Rio de Janeiro, Brazil. The appendix is stored in the documentation folder.

The model code is in the file called ode_epi_HIV_TB_code.R
The calibration code is in the file called calibration_epi_HIV_TB_code.R

The model and calibration codes were run with R version 3.5.2. Before running the code you will need to install the following packages
- 'dplyr'
- 'deSolve'
- 'readxl'
- 'stringr'
- 'reshape2'
- 'ggplot2'
- 'varhandle'
- 'here'
- 'paletteer'

If you would like to use your own parameter files, you can simply update the data input files in the param_files folder:
- If you want to update the calibration parameters you can update Calibration_parameters.xlsx
- If you want to update the epi model parameter you can update Epi_model_parameters.xlsx
