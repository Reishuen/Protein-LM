
import torch 
import torch.nn as nn
import torch.optim as optim

import numpy as np
import matplotlib.pyplot as plt

X = np.random.rand(100, 1)
y = 2 * X + np.random.randn(100,1)

model = LinearRegression()
model.fit(X,y)

predictions = model.Predict(X)

plt.scatter(X, y, color='blue', label='Data points')
plt.plot(X, predictions, color='red', label='Regression line')
plt.xlabel('X')
plt.ylabel('y')
plt.title('Linear Regression Fit')
plt.legend()
plt.show()