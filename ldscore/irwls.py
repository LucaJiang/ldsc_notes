# Python Class for Iterative Weighted Least Squares

import numpy as np

# import jackknife as jk

# TODO: add jackknife


class IRLS:
    """
    Python Class for Iterative Weighted Least Squares
    Parameters:
    :param X: numpy array of shape (n_samples, 1)
    :param y: numpy array of shape (n_samples, 1)
    :param fix_intercept: bool, whether to fix the intercept to 0, default False
    :param weights: numpy array of shape (n_samples, )

    Usage:
    ```
    from irwls import IRLS
    irwls = IRLS(X, y)
    irwls.regression()
    reg_intercept = irwls.get_intercept()
    reg_coefficients = irwls.get_coefficients()
    ```
    """

    def __init__(
        self,
        X,
        y,
        fix_intercept=False,
        weights=None,
        max_iter=1,
    ):
        """
        :param X: numpy array of shape (n_samples, n_features)
        :param y: numpy array of shape (n_samples, )
        :param fix_intercept: bool, whether to fix the intercept to 0, default False
        :param weights: numpy array of shape (n_samples, )
        :param max_iter: int, maximum number of iterations
        :param tol: float, tolerance for stopping criteria
        """
        if not fix_intercept:
            self.X = np.concatenate([np.ones_like(X), X], axis=1)
        else:
            self.X = X
        self.y = y
        self.fix_intercept = fix_intercept
        if weights is not None:
            self.weights = weights
        else:
            self.weights = np.ones_like(y)
        self.max_iter = max_iter
        self.p = 1

    def _get_initial_beta(self):
        """
        Initialize the regression coefficient with one
        """
        # return np.ones_like(self.X[0])
        return np.linalg.pinv(self.X) @ self.y

    def _get_residuals(self, X, y, beta):
        return y - X @ beta

    def _get_weights_r(self, residuals):
        """
        Calculate the weight vector with residuals
        ref: initial weighted method of irwls
        """
        weights = np.abs(residuals) ** (self.p / 2 - 1)
        weights /= weights.sum()
        return weights

    def _get_weights_yhat(self, pred):
        """
        Calculate the weight vector with residuals
        ref: ldsc/ldscore/regression.py-line: 498
        In ldsc.r: weight = 1 / (2 * pred ** 2)
        """
        weights = pred ** (-1 / 2) / 2 / self.weights
        return weights

    def predict(self, X):
        return X @ self.beta

    def _get_beta(self, X, y, weights):
        """
        Update the regression coefficient with weighted least squares
        return: numpy array of shape (n_features, 1)
        """
        return np.linalg.lstsq(
            np.multiply(X, weights), np.multiply(y, weights), rcond=None
        )[0].reshape(-1, 1)

    def _check_convergence(self, residuals, new_residuals):
        return np.linalg.norm(residuals - new_residuals) < self.tol

    def regression(self):
        """
        Fit the model and return the regression coefficient
        """
        # Initialize the regression coefficient and residuals
        self.beta = self._get_initial_beta()
        # Iterate
        for _ in range(self.max_iter):
            # residuals = self._get_residuals(self.X, self.y, beta)
            # Calculate the weight vector
            weights = self._get_weights_yhat(self.predict(self.X))
            # Update the regression coefficient with weighted least squares
            self.beta = self._get_beta(self.X, self.y, weights)
        # Save results
        # self.residuals = residuals

    def get_intercept(self):
        """
        Return the intercept
        """
        return self.beta[0].item()

    def get_coefficients(self):
        """
        Return the coefficients
        """
        if self.fix_intercept:
            return self.beta.item()
        return self.beta[1].item()


#! TEST
if __name__ == "__main__":
    import pandas as pd

    coef = pd.read_csv("./results/test_coef.txt", sep="\t")
    # x = coef["x"].values.reshape(-1, 1)
    y = coef["Z^2"].values.reshape(-1, 1)
    l2 = coef["L2"].values.reshape(-1, 1)
    # x = l2 * 61220 / 853604  # l2 * N / M
    x = l2 * 61220 / 1173569  # l2 * N / M
    irwls = IRLS(x, y, weights=l2)
    irwls.regression()
    print("Intercept = %.4f" % irwls.get_intercept())
    print("h_g^2 = %.4f" % irwls.get_coefficients())

    # irwls_fix = IRLS(x, y - irwls.get_intercept(), fix_intercept=True)
    # irwls_fix.regression()
    # print("h_g^2 fix = %.4f" % irwls_fix.get_coefficients())
