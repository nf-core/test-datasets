# ruff: noqa: PGH004
# ruff: noqa
# mypy: ignore-errors
"""This file contains the PyTorch model for the performance test of the Titanic dataset."""

from typing import Callable, Optional

import torch


class ModelTitanicPerformance(torch.nn.Module):
    def __init__(self, in_features: int = 16, out_features: int = 1, lyrs: list = [8]):
        super(ModelTitanicPerformance, self).__init__()
        activation_linear = torch.nn.Identity()
        activation_output = torch.nn.Identity()
        layers = []
        for i in range(len(lyrs)):
            layers.append(torch.nn.Linear(in_features, lyrs[i]))
            layers.append(activation_linear)
            in_features = lyrs[i]
        layers.append(torch.nn.Linear(in_features, out_features))
        layers.append(activation_output)
        self.layers = torch.nn.ModuleList(layers)

    def forward(
        self,
        Age: torch.Tensor,
        Family_Size: torch.Tensor,
        Fare: torch.Tensor,
        Parch: torch.Tensor,
        Pclass: torch.Tensor,
        Sex: torch.Tensor,
        SibSp: torch.Tensor,
        Embarked_C: torch.Tensor,
        Embarked_Q: torch.Tensor,
        Embarked_S: torch.Tensor,
        Title_Dr: torch.Tensor,
        Title_Master: torch.Tensor,
        Title_Miss: torch.Tensor,
        Title_Mr: torch.Tensor,
        Title_Mrs: torch.Tensor,
        Title_Rev: torch.Tensor,
        **kwargs: torch.Tensor,  # noqa: ARG002
    ):
        x = torch.stack(
            [
                Age,
                Family_Size,
                Fare,
                Parch,
                Pclass,
                Sex,
                SibSp,
                Embarked_C,
                Embarked_Q,
                Embarked_S,
                Title_Dr,
                Title_Master,
                Title_Miss,
                Title_Mr,
                Title_Mrs,
                Title_Rev,
            ],
            dim=1,
        )
        x = x.squeeze(-1)
        for layer in self.layers:
            x = layer(x)
        return x

    def compute_loss(self, output: torch.Tensor, Survived: torch.Tensor, loss_fn: Callable) -> torch.Tensor:
        """Compute the loss.

        Args:
            output: Model output tensor of shape [batch_size, nb_classes]
            survived: Target tensor of shape [batch_size, 1]
            loss_fn: Loss function (BCEWithLogitsLoss)

        Returns:
            Loss value
        """
        # Ensure both tensors have compatible shapes
        target = Survived.view(-1).float()  # Flatten and ensure float type
        output = output.view(-1).float()            # Flatten output as well
        
        # Ensure they have the same size
        if output.size(0) != target.size(0):
            raise RuntimeError(f"Size mismatch: output has {output.size(0)} elements, target has {target.size(0)} elements")
        
        try:
            return loss_fn(output, target)
        except (ValueError, RuntimeError) as e:
            raise RuntimeError(f"Error computing loss: {e}, output shape: {output.shape}, target shape: {target.shape}") from e

    def compute_accuracy(self, output: torch.Tensor, Survived: torch.Tensor) -> torch.Tensor:
        """Compute the accuracy.

        Args:
            output: Model output tensor of shape [batch_size, nb_classes]
            survived: Target tensor of shape [batch_size, 1]
        """
        # Squeeze the extra dimension from the target tensor and ensure long dtype
        target = Survived.view(-1)
        # Compute the accuracy
        accuracy = ((output > 0) == target).float().mean()
        return accuracy

    def batch(
        self,
        batch: dict[str, torch.Tensor],
        loss_fn: Callable = torch.nn.BCEWithLogitsLoss(),
        optimizer: Optional[torch.optim.Optimizer] = None,
    ) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
        """Perform one batch step.

        `batch` is a dictionary with the input and label tensors.
        `loss_fn` is the loss function to be used.

        If `optimizer` is passed, it will perform the optimization step -> training step
        Otherwise, only return the forward pass output and loss -> evaluation step
        """
        # Forward pass
        output = self.forward(**batch).squeeze(-1)

        # Compute loss
        loss = self.compute_loss(output, batch["Survived"], loss_fn=loss_fn)

        # Backward pass and optimization
        if optimizer is not None:
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        accuracy = self.compute_accuracy(output, batch["Survived"])
        return loss, {"accuracy": accuracy, "predictions": output}
