import numpy as np

import torch
import torch.nn as nn
import torch.nn.functional as F

### Dense Neural Network ###
class DenseNN(nn.Module):

    '''
    Class for the dense (fully-connected) neural network decoder
    Parameters
    ----------
    units: integer or vector of integers, optional, default 100
        This is the number of hidden units in each layer
        If you want a single layer, input an integer (e.g. units=400 will give you a single hidden layer with 400 units)
        If you want multiple layers, input a vector (e.g. units=[400,200]) will give you 2 hidden layers with 400 and 200 units, repsectively.
        The vector can either be a list or an array

    adapted from https://github.com/KordingLab/Neural_Decoding from tensorflow to pytorch since that's what i know
    '''
    
    def __init__(self, in_features, out_features, units=100, act=nn.ReLU, dropout=0):
        super(DenseNN, self).__init__()
        self.out_features = out_features
        self.in_features = in_features
        self.act = act()
        self.dropout = dropout

        # if 'units' is an integer, put it in the form of a vector
        try:
            units[0]
        except:
            units = [units]
        self.units = units

        # input layer
        self.fc1 = nn.Linear(self.in_features , self.units[0])

        # number of hidden layers (determined by length of units)
        self.num_layers = len(self.units)
        self.hidden = nn.ModuleList()
        for k in range(self.num_layers-1):
            self.hidden.append(nn.Linear(self.units[k], self.units[k+1]))

        # output layer
        self.out = nn.Linear(self.units[-1] , self.out_features)


    def forward(self, x): # x represents our data
        x = self.act(self.fc1(x))
        if self.dropout != 0:
            drop = nn.Dropout(p = self.dropout)
            x = drop(x)
        for hidden_layer in self.hidden:
            x = self.act(hidden_layer(x))
            if self.dropout != 0:
                drop = nn.Dropout(p = self.dropout)
                x = drop(x)   
        # x = torch.sigmoid(self.out(x))
        x = self.out(x)

        return x


    def fit(self, loss_func, optimizer, X, Y, epochs=100):
        for i in range(epochs):
            preds = self.forward(X) ## Make Predictions by forward pass through network

            loss = loss_func(preds, Y) ## Calculate Loss

            optimizer.zero_grad() ## Zero weights before calculating gradients
            loss.backward() ## Calculate Gradients
            optimizer.step() ## Update Weights

            if i % 500 == 0: ## Print MSE every 500 epochs
                print(f"Epoch: {i}  ||  MSE: {loss}")


    def predict(self, x):
        return self.forward(x)


