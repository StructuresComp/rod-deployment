import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from planner.utils import float_tensor, get_numpy, ON_GPU


class PoseInterpolator(nn.Module):
    def __init__(self, hidden_layer_units=(392,)*4):
        super(PoseInterpolator, self).__init__()
        self.input_size = 3
        self.output_size = 6
        self.layers = nn.ModuleList()
        self.layers.append(nn.Linear(self.input_size, hidden_layer_units[0]))
        for i in range(1, len(hidden_layer_units)):
            self.layers.append(nn.Linear(hidden_layer_units[i-1], hidden_layer_units[i]))
        self.layers.append(nn.Linear(hidden_layer_units[-1], self.output_size))

        self.stds = np.array([[1.80860715e+01, 5.33914771e+03, 1.07518554e+00]], dtype="float32").reshape((1, 3))
        self.means = np.array([[1.79961462e+01, 3.31403698e+03, 1.57740136e+00]], dtype="float32").reshape((1, 3))

        # Load pretrained weights
        filename = "392x4/model_(392, 392, 392, 392)_bs_2048_mse_0.00002_mae_0.00141"

        if torch.cuda.is_available():
            self.load_state_dict(torch.load("deployment_models/" + filename + ".pth"))
        else:
            self.load_state_dict(torch.load("deployment_models/" + filename + ".pth", map_location=torch.device('cpu')))
        self.eval()

    def get_data_params(self):
        return self.means, self.stds

    def preprocess_input(self, orig_data):
        data = orig_data.copy()
        data[:] -= self.means
        data[:] /= self.stds
        return data, float_tensor(data)

    def revert_data(self, data):
        data[:] *= self.stds
        data[:] += self.means

    def forward(self, state):
        x_c = state
        for i in range(len(self.layers)-1):
            x_c = F.relu(self.layers[i](x_c))
        output = self.layers[-1](x_c)
        return output
