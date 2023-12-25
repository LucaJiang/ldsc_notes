# Python Class for Iterative Weighted Least Squares

import numpy as np


class IRLS:
    """
    Python Class for Iterative Weighted Least Squares
    Algorithm:
    1. Compute an initial regression estimate $\hat{\boldsymbol{\beta}}_0$.
    2. For $k=0,1, \ldots$, compute the residuals $r_{i, k}\left(\hat{\boldsymbol{\beta}}_k\right)=y_i-\boldsymbol{x}_i \hat{\boldsymbol{\beta}}_k$ and weight matrix $\boldsymbol{W}_k=\operatorname{diag}\left(w_{1, k}, \ldots, w_{n, k}\right)$, where $w_{i, k}=w\left(r_{i, k}\left(\hat{\boldsymbol{\beta}}_k\right)\right)$ with some weight function $w$. Then update $\hat{\boldsymbol{\beta}}_{k+1}$ using WLS as in (2) with weight matrix $\boldsymbol{W}_{k^{\prime}}$ that is,
        $$
        \hat{\boldsymbol{\beta}}_{k+1}=\left(\boldsymbol{X}^{\top} \boldsymbol{W}_k \boldsymbol{X}\right)^{-1} \boldsymbol{X}^{\top} \boldsymbol{W}_k \boldsymbol{y}
        $$
    3. Stop when $\max _i\left|r_{i, k}-r_{i, k+1}\right|<\epsilon$, where $\epsilon$ may be fixed or related to the residual scale.
    ref: https://link.springer.com/referenceworkentry/10.1007/978-3-030-26050-7_169-1#Equ2
    """

    def __init__(self, X, y, weights=None, max_iter=1000, tol=1e-6):
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

    def fit(self):
        """
        Fit the model
        :return: self
        """
        ## initialize
        self._initialize()
        ## iterate
        for i in range(self.max_iter):
            ## update beta
            self._update_beta()
            ## update weights
            self._update_weights()
            ## check convergence
            if self._check_convergence():
                break
        ## return
        return self

    def predict(self, X):
        """
        Predict y
        :param X: numpy array of shape (n_samples, n_features)
        :return: numpy array of shape (n_samples, )
        """
        return np.dot(X, self.beta)

    def _initialize(self):
        """
        Initialize beta and weights
        :return: None
        """
        ## initialize beta
        self.beta = np.zeros(self.X.shape[1])
        ## initialize weights
        if self.weights is None:
            self.weights = np.ones(self.X.shape[0])

    def _update_beta(self):
        """
        Update beta
        :return: None
        """
        ## calculate z
        z = self._calculate_z()
        ## calculate w
        w = self._calculate_w(z)
        ## calculate beta
        self.beta = self._calculate_beta(w, z)

    def _calculate_z(self):
        """
        Calculate z
        :return: numpy array of shape (n_samples, )
        """
        return self.y + np.dot(self.X, self.beta) - self._link_function(self.y)

    def _calculate_w(self, z):
        """
        Calculate w
        :param z: numpy array of shape (n_samples, )
        :return: numpy array of shape (n_samples, )
        """
        return self
