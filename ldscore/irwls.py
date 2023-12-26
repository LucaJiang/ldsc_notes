# Python Class for Iterative Weighted Least Squares

import numpy as np
import jackknife as jk

# TODO: add jackknife


class IRLS:
    """
    Python Class for Iterative Weighted Least Squares
    Usage:
    ```
    from irwls import IRLS
    irwls = IRLS(X, y)
    irwls.regression()
    reg_intercept = irwls.get_intercept()
    reg_coefficients = irwls.get_coefficients()
    ```
    Algorithm:
      1. Calculate the residuals:
        $$\mathbf{r}_k = \mathbf{Y} - \mathbf{X}\hat{\boldsymbol{\beta}}_k$$
      2. Calculate the weight matrix:
        $$\begin{aligned}
        \mathbf{w}_k &= \| \mathbf{r}_k \|^{p-2}\\
        \mathbf{W}_k &= \text{diag}(\mathbf{w}_k/\sum \mathbf{w}_k)
        \end{aligned}$$
      3. Update the regression coefficient:
        $$\hat{\boldsymbol{\beta}}_{k+1} = (\mathbf{X}^T\mathbf{W}_k\mathbf{X})^{-1}\mathbf{X}^T\mathbf{W}_k\mathbf{Y}$$
      4. Check convergence:
        $$\|\mathbf{r}_k - \mathbf{r}_{k+1}\| < \epsilon$$
    """

    def __init__(self, X, y, weights=None, max_iter=10, tol=1e-2):
        """
        :param X: numpy array of shape (n_samples, n_features)
        :param y: numpy array of shape (n_samples, )
        :param weights: numpy array of shape (n_samples, )
        :param max_iter: int, maximum number of iterations
        :param tol: float, tolerance for stopping criteria
        """
        self.X = X
        self.y = y
        self.weights = weights
        self.max_iter = max_iter
        self.tol = tol
        self.p = 1

    def _get_initial_beta(self):
        """
        Initialize the regression coefficient with pseudo-inverse
        """
        return np.linalg.pinv(self.X) @ self.y

    def _get_residuals(self, X, y, beta):
        return y - X @ beta

    def _get_weights(self, residuals):
        """
        Calculate the weight matrix
        """
        weights = np.abs(residuals) ** (self.p - 2)
        assert np.all(weights >= 0)  # weights should be non-negative
        weights /= weights.sum()
        return np.diag(weights)

    def _get_beta(self, X, y, weights):
        """
        Update the regression coefficient
        """
        return np.linalg.inv(X.T @ weights @ X) @ X.T @ weights @ y

    def _check_convergence(self, residuals, new_residuals):
        return np.linalg.norm(residuals - new_residuals) < self.tol

    def regression(self):
        """
        Fit the model and return the regression coefficient
        """
        # Initialize the regression coefficient and residuals
        beta = self._get_initial_beta()
        old_residuals = np.zeros_like(self.y) + np.inf
        convergence = False
        # Iterate
        for i in range(self.max_iter):
            residuals = self._get_residuals(self.X, self.y, beta)
            if self._check_convergence(residuals, old_residuals):
                convergence = True
                break
            old_residuals = residuals
            # Calculate the weight matrix
            weights = self._get_weights(residuals)
            # Update the regression coefficient
            beta = self._get_beta(self.X, self.y, weights)
        # Save results
        self.beta = beta
        self.residuals = residuals
        self.weights = weights
        self.n_iter = i + 1
        self.convergence = convergence

    def get_intercept(self):
        """
        Return the intercept
        """
        return self.beta[0]

    def get_coefficients(self):
        """
        Return the coefficients
        """
        return self.beta[1]
