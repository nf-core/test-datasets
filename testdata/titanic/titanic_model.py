"""Titanic dataset model implementation."""

from typing import Callable, Optional

import torch
from torch import nn
from torch.optim import Optimizer


class ModelTitanic(torch.nn.Module):
    """A simple model for Titanic dataset."""

    def __init__(
        self,
        nb_neurons_intermediate_layer: int = 7,
        nb_intermediate_layers: int = 3,
        nb_classes: int = 2,
    ) -> None:
        """Initialize model layers.

        Args:
            nb_neurons_intermediate_layer: Number of neurons in intermediate layers
            nb_intermediate_layers: Number of intermediate layers
            nb_classes: Number of output classes
        """
        super().__init__()
        self.input_layer = nn.Linear(7, nb_neurons_intermediate_layer)
        self.intermediate = nn.modules.ModuleList(
            [
                nn.Linear(nb_neurons_intermediate_layer, nb_neurons_intermediate_layer)
                for _ in range(nb_intermediate_layers)
            ],
        )
        self.output_layer = nn.Linear(nb_neurons_intermediate_layer, nb_classes)
        self.relu = nn.ReLU()
        self.softmax = nn.Softmax(dim=1)

    def forward(
        self,
        pclass: torch.Tensor,
        sex: torch.Tensor,
        age: torch.Tensor,
        sibsp: torch.Tensor,
        parch: torch.Tensor,
        fare: torch.Tensor,
        embarked: torch.Tensor,
    ) -> torch.Tensor:
        """Forward pass of the model.

        Args:
            pclass: Tensor of shape [batch_size, 1]
            sex: Tensor of shape [batch_size, 1]
            ...etc

        Returns:
            Tensor of shape [batch_size, nb_classes] containing class probabilities
        """
        # Stack features and remove the extra dimension
        x = torch.stack((pclass, sex, age, sibsp, parch, fare, embarked), dim=1).float()  # [batch_size, 7, 1]
        x = x.squeeze(-1)  # [batch_size, 7]

        # Pass through layers
        x = self.relu(self.input_layer(x))
        for layer in self.intermediate:
            x = self.relu(layer(x))
        return self.softmax(self.output_layer(x))

    def compute_loss(self, output: torch.Tensor, survived: torch.Tensor, loss_fn: Callable) -> torch.Tensor:
        """Compute the loss.

        Args:
            output: Model output tensor of shape [batch_size, nb_classes]
            survived: Target tensor of shape [batch_size, 1]
            loss_fn: Loss function (CrossEntropyLoss)

        Returns:
            Loss value
        """
        # Squeeze the extra dimension from the target tensor and ensure long dtype
        target = survived.squeeze(-1).long()
        return loss_fn(output, target)

    def batch(
        self,
        x: dict[str, torch.Tensor],
        y: dict[str, torch.Tensor],
        loss_fn: Callable,
        optimizer: Optional[Optimizer] = None,
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        """Perform one batch step.

        `x` is a dictionary with the input tensors.
        `y` is a dictionary with the target tensors.
        `loss_fn` is the loss function to be used.

        If `optimizer` is passed, it will perform the optimization step -> training step
        Otherwise, only return the forward pass output and loss -> evaluation step
        """
        output = self.forward(**x)
        loss = self.compute_loss(output, **y, loss_fn=loss_fn)

        if optimizer is not None:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        return loss, {"output": output}
