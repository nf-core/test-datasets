import torch
import torch.nn as nn
from typing import Callable, Optional, Tuple

class ModelSimple(torch.nn.Module):
    """
    A simple model example.
    It takes as input a 1D tensor of any size,
    apply some convolutional layer and
    outputs a single value using a maxpooling layer and a softmax function.

    All functions `forward`, `compute_loss` and `batch` need to be implemented for any new model.
    """
    def __init__(self, kernel_size: int = 3, pool_size: int = 2):
        super(ModelSimple, self).__init__()
        self.conv1 = nn.Conv1d(in_channels=4, out_channels=1, kernel_size=kernel_size)
        self.pool = nn.MaxPool1d(pool_size, pool_size)
        self.softmax = nn.Softmax(dim=1)
        # had to change to 6 because dna sequence is shoprter
        self.linear = nn.Linear(6, 1)

    def forward(self, hello: torch.Tensor) -> dict:
        """
        Forward pass of the model.
        It should return the output as a dictionary, with the same keys as `y`.
        """
        x = hello.permute(0, 2, 1).to(torch.float32)  # permute the two last dimensions of hello
        x = self.conv1(x)
        x = self.pool(x)
        x = self.softmax(x)
        x = self.linear(x)
        x = x.squeeze()
        return x

    def compute_loss(self, output: torch.Tensor, hola: torch.Tensor, loss_fn: Callable) -> torch.Tensor:
        """
        Compute the loss.
        `output` is the output tensor of the forward pass.
        `hola` is the target tensor -> label column name.
        `loss_fn` is the loss function to be used.
        """
        return loss_fn(output, hola.to(torch.float32))

    def batch(self, x: dict, y: dict, loss_fn1: Callable, loss_fn2: Callable, optimizer: Optional[Callable] = None) -> Tuple[torch.Tensor, dict]:
        """
        Perform one batch step.
        `x` is a dictionary with the input tensors.
        `y` is a dictionary with the target tensors.
        `loss_fn1` and `loss_fn2` are the loss function to be used.

        If `optimizer` is passed, it will perform the optimization step -> training step
        Otherwise, only return the forward pass output and loss -> evaluation step

        TODO currently only returning loss1, but we could potentially summarize loss1 and loss2 in some way.
        However, note that both loss1 and loss2 are participating in the backward propagation, one after another.
        """
        output = self(**x)
        loss1 = self.compute_loss(output, **y, loss_fn=loss_fn1)
        loss2 = self.compute_loss(output, **y, loss_fn=loss_fn2)
        if optimizer is not None:
            optimizer.zero_grad()
            loss1.backward(retain_graph=True)
            loss2.backward(retain_graph=True)
            optimizer.step()
        return loss1, output
