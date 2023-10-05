# ROS_evaluation
Selected scripts relating to an evaluation of Australian headfire rate of spread models. 

empirical_models.py is a python translation of Vesta Mk1, Vesta Mk2, 10% windspeed rule of thumb, and McArthur Mk5 rate of spread models. Requires weather and fuel input parameters that are not able to be distributed.  

ros_analysis.R generates plots and RF models from dataframe containing ros predictions, ros observations, and weather parameters. 
