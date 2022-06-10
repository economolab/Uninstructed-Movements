import numpy as np
import torch

########## R-squared (R2) ##########

def get_R2(y,y_pred):

    """
    Function to get R2
    Parameters
    ----------
    y - the true outputs (a matrix of size number of examples x number of outputs)
    y_pred - the predicted outputs (a matrix of size number of examples x number of outputs)
    Returns
    -------
    R2_array: An array of R2s for each output

    adapted from : https://github.com/KordingLab/Neural_Decoding/blob/master/Neural_Decoding/metrics.py
    """

    R2_list=np.zeros((y.shape[1],1)) #Initialize a list that will contain the R2s for all the outputs
    for i in range(y.shape[1]): #Loop through outputs
        #Compute R2 for each output
        y_mean=torch.mean(y[:,i])
        R2=1-torch.sum( (y_pred[:,i]-y[:,i])**2) / (torch.sum((y[:,i]-y_mean)**2) )
        R2_list[i] = R2.detach().numpy() #Append R2 of this output to the list
    return R2_list #Return an array of R2s
